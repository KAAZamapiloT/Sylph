/**
 * @file azimuthal_fft.hpp
 * @brief FFT-based Azimuthal Projection & Synthesis for Spherical Harmonics.
 *
 * Replaces the O(N_theta * N_phi * M) direct trigonometric loops in
 * SHGridTransform (Eq. 20-25 forward synthesis and Eq. 37 inverse projection)
 * with O(N_theta * N_phi * log N_phi) Fast Fourier Transforms:
 *
 *  - Power-of-2 N_phi: iterative Radix-2 Cooley-Tukey FFT.
 *  - Arbitrary  N_phi: Bluestein's Chirp-Z Transform (exact for any N_phi >= 1,
 *                      including primes, via power-of-2 cyclic convolution).
 *  - Paired-row packing: transforms two real theta-rows simultaneously inside
 *                        a single complex FFT, halving transform work.
 */
#pragma once

#include "sylph/spherical_grid.hpp"
#include "types/constants.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace sylph {

class AzimuthalFFT {
public:
    using Complex = std::complex<double>;

    struct Plan {
        int n = 0;
        int m = 0;            // FFT size (== n if is_pow2, else next_pow2(2*n - 1))
        bool is_pow2 = false;

        std::vector<int> bitrev;
        // Packed stage twiddles: for stage half = 1, 2, 4, ..., m/2,
        // entries [half - 1 .. 2*half - 2] hold exp(-/+ i * pi * k / half).
        std::vector<Complex> twiddles_fwd;
        std::vector<Complex> twiddles_inv;

        // Bluestein precomputed tables (only populated when !is_pow2)
        std::vector<Complex> chirp;       // exp(-i * pi * k^2 / n)
        std::vector<Complex> kernel_fft;  // FFT_m(b), where b[k] = exp(+i * pi * k^2 / n)
    };

    static constexpr bool is_power_of_two(int n) noexcept {
        return n > 0 && (n & (n - 1)) == 0;
    }

    static constexpr int next_pow2(int n) noexcept {
        if (n <= 1) return 1;
        int p = 1;
        while (p < n) {
            p <<= 1;
        }
        return p;
    }

    static const Plan& get_plan(int n) {
        if (n <= 0) {
            throw std::invalid_argument("AzimuthalFFT length must be positive");
        }

        thread_local std::unordered_map<int, Plan> cache;
        auto it = cache.find(n);
        if (it != cache.end()) {
            return it->second;
        }

        auto [inserted, _] = cache.emplace(n, build_plan(n));
        return inserted->second;
    }

    /**
     * @brief In-place iterative Radix-2 Cooley-Tukey FFT of length plan.m.
     *
     * Forward (inverse = false): X[r] = sum_{j=0}^{m-1} x[j] * exp(-i * 2*pi * r * j / m)
     * Inverse (inverse = true):  X[r] = sum_{j=0}^{m-1} x[j] * exp(+i * 2*pi * r * j / m)
     * (Unnormalized in both directions.)
     */
    static void fft_pow2_inplace(
        Complex* data,
        const Plan& plan,
        bool inverse
    ) noexcept {
        const int m = plan.m;
        if (m <= 1) return;

        const int* bitrev = plan.bitrev.data();
        for (int i = 0; i < m; ++i) {
            const int j = bitrev[i];
            if (i < j) {
                std::swap(data[i], data[j]);
            }
        }

        const Complex* twiddles = inverse
            ? plan.twiddles_inv.data()
            : plan.twiddles_fwd.data();

        for (int half = 1; half < m; half <<= 1) {
            const int step = half << 1;
            const Complex* w_stage = twiddles + (half - 1);

            for (int start = 0; start < m; start += step) {
                for (int k = 0; k < half; ++k) {
                    const Complex u = data[start + k];
                    const Complex v = data[start + k + half] * w_stage[k];
                    data[start + k]        = u + v;
                    data[start + k + half] = u - v;
                }
            }
        }
    }

    /**
     * @brief Arbitrary-length unnormalized complex DFT / IDFT of length plan.n.
     *
     * Uses Radix-2 directly when plan.is_pow2, otherwise Bluestein's Chirp-Z
     * transform via power-of-2 cyclic convolution of size plan.m.
     */
    static void transform(
        const Complex* in,
        Complex* out,
        const Plan& plan,
        bool inverse,
        std::vector<Complex>& scratch
    ) {
        const int n = plan.n;
        if (plan.is_pow2) {
            if (in != out) {
                for (int i = 0; i < n; ++i) {
                    out[i] = in[i];
                }
            }
            fft_pow2_inplace(out, plan, inverse);
            return;
        }

        // Bluestein's Chirp-Z Transform:
        // Forward: X[k] = chirp[k] * sum_j (x[j] * chirp[j]) * conj(chirp[|k-j|])
        // Inverse: IDFT(x) = conj(DFT(conj(x)))
        const int m = plan.m;
        scratch.assign(static_cast<std::size_t>(m), Complex(0.0, 0.0));

        const Complex* chirp = plan.chirp.data();
        if (!inverse) {
            for (int j = 0; j < n; ++j) {
                scratch[static_cast<std::size_t>(j)] = in[j] * chirp[j];
            }
        } else {
            for (int j = 0; j < n; ++j) {
                scratch[static_cast<std::size_t>(j)] = std::conj(in[j]) * chirp[j];
            }
        }

        fft_pow2_inplace(scratch.data(), plan, false);

        const Complex* kernel = plan.kernel_fft.data();
        for (int r = 0; r < m; ++r) {
            scratch[static_cast<std::size_t>(r)] *= kernel[r];
        }

        fft_pow2_inplace(scratch.data(), plan, true);

        const double inv_m = 1.0 / static_cast<double>(m);
        if (!inverse) {
            for (int k = 0; k < n; ++k) {
                out[k] = scratch[static_cast<std::size_t>(k)] * (chirp[k] * inv_m);
            }
        } else {
            for (int k = 0; k < n; ++k) {
                out[k] = std::conj(scratch[static_cast<std::size_t>(k)] * (chirp[k] * inv_m));
            }
        }
    }

    /**
     * @brief Eq. 37 Azimuthal Projection (Inverse SH Transform) via FFT.
     *
     * Computes, for each i in [0, nTheta) and m in [-(output_order-1), output_order-1]:
     *   s[i, m] = (2*pi / nPhi) * sum_{j=0}^{nPhi-1} grid(i, j) * R(m, phi_j)
     *
     * Output buffer `s` is sized nTheta * (2 * output_order - 1) with
     * mIndex(m) = m + output_order - 1.
     */
    static void project_rows(
        const Spherical_grid& grid,
        int output_order,
        std::vector<double>& s
    ) {
        const int nTheta = grid.theta_size();
        const int nPhi   = grid.phi_size();
        const int mCount = 2 * output_order - 1;
        const int mZero  = output_order - 1;

        s.assign(
            static_cast<std::size_t>(nTheta) * static_cast<std::size_t>(mCount),
            0.0
        );

        const double phiScale = 2.0 * pi / static_cast<double>(nPhi);

        // Fast path for DC-only projection (e.g. product_integral where output_order == 1)
        if (output_order == 1) {
            for (int i = 0; i < nTheta; ++i) {
                double sum = 0.0;
                for (int j = 0; j < nPhi; ++j) {
                    sum += grid(i, j);
                }
                s[static_cast<std::size_t>(i)] = phiScale * sum;
            }
            return;
        }

        const Plan& plan = get_plan(nPhi);
        const double sqrt2Scale = phiScale * std::numbers::sqrt2_v<double>;

        std::vector<Complex> buf(static_cast<std::size_t>(nPhi));
        std::vector<Complex> scratch;
        if (!plan.is_pow2) {
            scratch.resize(static_cast<std::size_t>(plan.m));
        }

        // Pair two real rows (i0, i1) into z[j] = g_0[j] + i * g_1[j]
        for (int i0 = 0; i0 < nTheta; i0 += 2) {
            const int i1 = i0 + 1;
            const bool has_pair = (i1 < nTheta);

            if (has_pair) {
                for (int j = 0; j < nPhi; ++j) {
                    buf[static_cast<std::size_t>(j)] = Complex(grid(i0, j), grid(i1, j));
                }
            } else {
                for (int j = 0; j < nPhi; ++j) {
                    buf[static_cast<std::size_t>(j)] = Complex(grid(i0, j), 0.0);
                }
            }

            transform(buf.data(), buf.data(), plan, false, scratch);

            double* row0 = s.data() + static_cast<std::size_t>(i0) * static_cast<std::size_t>(mCount);
            double* row1 = has_pair
                ? (s.data() + static_cast<std::size_t>(i1) * static_cast<std::size_t>(mCount))
                : nullptr;

            // m = 0 (DC bin: Z[0] = G0[0] + i * G1[0])
            row0[mZero] = phiScale * buf[0].real();
            if (has_pair) {
                row1[mZero] = phiScale * buf[0].imag();
            }

            // m = +/- k for k = 1 .. output_order - 1
            // Using real-signal conjugate symmetry:
            //   G0[r] = 0.5 * (Z[r] + conj(Z[r']))
            //   G1[r] = 0.5 / i * (Z[r] - conj(Z[r']))
            // where r = k mod nPhi, r' = (nPhi - r) mod nPhi.
            // Then:
            //   s[i, +k] = (2*pi / nPhi) * sqrt(2) *  Re(G[r])
            //   s[i, -k] = (2*pi / nPhi) * sqrt(2) * (-Im(G[r]))
            for (int k = 1; k < output_order; ++k) {
                const int r = k % nPhi;
                const int r_conj = (r == 0) ? 0 : (nPhi - r);

                const Complex& zr = buf[static_cast<std::size_t>(r)];
                const Complex& zc = buf[static_cast<std::size_t>(r_conj)];

                const double re0 = 0.5 * (zr.real() + zc.real());
                const double im0 = 0.5 * (zr.imag() - zc.imag());

                row0[mZero + k] =  sqrt2Scale * re0;
                row0[mZero - k] = -sqrt2Scale * im0;

                if (has_pair) {
                    const double re1 = 0.5 * (zr.imag() + zc.imag());
                    const double im1 = 0.5 * (zc.real() - zr.real());

                    row1[mZero + k] =  sqrt2Scale * re1;
                    row1[mZero - k] = -sqrt2Scale * im1;
                }
            }
        }
    }

    /**
     * @brief Eq. 20-25 Azimuthal Synthesis (Forward SH Transform) via IFFT.
     *
     * Computes, for each i in [0, nTheta) and j in [0, nPhi):
     *   grid(i, j) = sum_{m = -(order-1)}^{order-1} s[i, m] * R(m, phi_j)
     *
     * Input buffer `s` is sized nTheta * (2 * order - 1) with
     * mIndex(m) = m + order - 1.
     */
    static void synthesize_rows(
        const std::vector<double>& s,
        int order,
        Spherical_grid& grid
    ) {
        const int nTheta = grid.theta_size();
        const int nPhi   = grid.phi_size();
        const int mCount = 2 * order - 1;
        const int mZero  = order - 1;

        if (order == 1) {
            for (int i = 0; i < nTheta; ++i) {
                const double dc = s[static_cast<std::size_t>(i)];
                for (int j = 0; j < nPhi; ++j) {
                    grid(i, j) = dc;
                }
            }
            return;
        }

        const Plan& plan = get_plan(nPhi);
        constexpr double inv_sqrt2 = 1.0 / std::numbers::sqrt2_v<double>;

        std::vector<Complex> spec(static_cast<std::size_t>(nPhi));
        std::vector<Complex> scratch;
        if (!plan.is_pow2) {
            scratch.resize(static_cast<std::size_t>(plan.m));
        }

        // Pack two real syntheses into one complex spectrum:
        //   S = S^(i0) + i * S^(i1)
        // Since IDFT(S^(i0)) and IDFT(S^(i1)) are purely real:
        //   Re(IDFT(S)) = grid(i0, .),  Im(IDFT(S)) = grid(i1, .)
        for (int i0 = 0; i0 < nTheta; i0 += 2) {
            const int i1 = i0 + 1;
            const bool has_pair = (i1 < nTheta);

            spec.assign(static_cast<std::size_t>(nPhi), Complex(0.0, 0.0));

            const double* row0 = s.data() + static_cast<std::size_t>(i0) * static_cast<std::size_t>(mCount);
            spec[0] += Complex(row0[mZero], 0.0);

            for (int k = 1; k < order; ++k) {
                const int r_pos = k % nPhi;
                const int r_neg = (r_pos == 0) ? 0 : (nPhi - r_pos);

                const double c0 = inv_sqrt2 * row0[mZero + k];
                const double d0 = inv_sqrt2 * row0[mZero - k];

                spec[static_cast<std::size_t>(r_pos)] += Complex(c0, -d0);
                spec[static_cast<std::size_t>(r_neg)] += Complex(c0,  d0);
            }

            if (has_pair) {
                const double* row1 = s.data() + static_cast<std::size_t>(i1) * static_cast<std::size_t>(mCount);
                spec[0] += Complex(0.0, row1[mZero]);

                for (int k = 1; k < order; ++k) {
                    const int r_pos = k % nPhi;
                    const int r_neg = (r_pos == 0) ? 0 : (nPhi - r_pos);

                    const double c1 = inv_sqrt2 * row1[mZero + k];
                    const double d1 = inv_sqrt2 * row1[mZero - k];

                    // Multiply (c1 - i*d1) by i -> (d1 + i*c1)
                    // Multiply (c1 + i*d1) by i -> (-d1 + i*c1)
                    spec[static_cast<std::size_t>(r_pos)] += Complex( d1, c1);
                    spec[static_cast<std::size_t>(r_neg)] += Complex(-d1, c1);
                }
            }

            transform(spec.data(), spec.data(), plan, true, scratch);

            for (int j = 0; j < nPhi; ++j) {
                grid(i0, j) = spec[static_cast<std::size_t>(j)].real();
            }
            if (has_pair) {
                for (int j = 0; j < nPhi; ++j) {
                    grid(i1, j) = spec[static_cast<std::size_t>(j)].imag();
                }
            }
        }
    }

private:
    static Plan build_plan(int n) {
        Plan plan;
        plan.n = n;
        plan.is_pow2 = is_power_of_two(n);
        plan.m = plan.is_pow2 ? n : next_pow2(2 * n - 1);

        const int m = plan.m;

        // 1. Bit-reversal table for size m
        plan.bitrev.resize(static_cast<std::size_t>(m), 0);
        int log_m = 0;
        while ((1 << log_m) < m) {
            ++log_m;
        }
        for (int i = 0; i < m; ++i) {
            int rev = 0;
            for (int b = 0; b < log_m; ++b) {
                if (i & (1 << b)) {
                    rev |= (1 << (log_m - 1 - b));
                }
            }
            plan.bitrev[static_cast<std::size_t>(i)] = rev;
        }

        // 2. Stage twiddle factors for size m
        if (m > 1) {
            plan.twiddles_fwd.resize(static_cast<std::size_t>(m - 1));
            plan.twiddles_inv.resize(static_cast<std::size_t>(m - 1));

            for (int half = 1; half < m; half <<= 1) {
                const double base_angle = -pi / static_cast<double>(half);
                for (int k = 0; k < half; ++k) {
                    const double angle = base_angle * static_cast<double>(k);
                    const double c = std::cos(angle);
                    const double s = std::sin(angle);
                    const std::size_t idx = static_cast<std::size_t>(half - 1 + k);
                    plan.twiddles_fwd[idx] = Complex(c,  s);
                    plan.twiddles_inv[idx] = Complex(c, -s);
                }
            }
        }

        // 3. Bluestein chirp & pre-transformed kernel (only when n is not pow2)
        if (!plan.is_pow2) {
            plan.chirp.resize(static_cast<std::size_t>(n));
            const std::uint64_t two_n = 2ULL * static_cast<std::uint64_t>(n);

            for (int k = 0; k < n; ++k) {
                const std::uint64_t kk =
                    (static_cast<std::uint64_t>(k) * static_cast<std::uint64_t>(k)) % two_n;
                const double angle = -pi * static_cast<double>(kk) / static_cast<double>(n);
                plan.chirp[static_cast<std::size_t>(k)] = Complex(std::cos(angle), std::sin(angle));
            }

            plan.kernel_fft.assign(static_cast<std::size_t>(m), Complex(0.0, 0.0));
            plan.kernel_fft[0] = Complex(1.0, 0.0);
            for (int k = 1; k < n; ++k) {
                const Complex bk = std::conj(plan.chirp[static_cast<std::size_t>(k)]);
                plan.kernel_fft[static_cast<std::size_t>(k)]     = bk;
                plan.kernel_fft[static_cast<std::size_t>(m - k)] = bk;
            }
            fft_pow2_inplace(plan.kernel_fft.data(), plan, false);
        }

        return plan;
    }
};

} // namespace sylph
