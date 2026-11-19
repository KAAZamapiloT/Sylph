#include "verification_example.hpp"

#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/gauss_legendre.hpp"
#include "sylph/azimuthal_fft.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <types/constants.hpp>

namespace examples
{

static sylph::SHCoefficients naive_product(
    const sylph::SHCoefficients& a,
    const sylph::SHCoefficients& b,
    int output_order)
{
    int full_order = a.order() + b.order() - 1;
    int tt = full_order + output_order - 1;
    int nTheta = (tt + 1) / 2;
    int nPhi = tt;

    sylph::gauss_legendre gl(nTheta);
    std::vector<double> weights = gl.weights();
    std::vector<double> nodes = gl.nodes();

    sylph::SHCoefficients result(output_order);
    sylph::SphericalHarmonics sh;

    const double dPhi = 2.0 * sylph::pi / nPhi;

    for (int i = 0; i < nTheta; ++i) {
        double theta = std::acos(nodes[i]);
        for (int j = 0; j < nPhi; ++j) {
            double phi = j * dPhi;
            double valA = sh.reconstruct(a, theta, phi);
            double valB = sh.reconstruct(b, theta, phi);

            for (int l = 0; l < output_order; ++l) {
                for (int m = -l; m <= l; ++m) {
                    double y = sh.evaluate(l, m, theta, phi);
                    result(l, m) += valA * valB * y * weights[i] * dPhi;
                }
            }
        }
    }
    return result;
}

static void direct_azimuthal_project(
    const sylph::Spherical_grid& grid,
    int output_order,
    std::vector<double>& s)
{
    const int nTheta = grid.theta_size();
    const int nPhi   = grid.phi_size();
    const int mCount = 2 * output_order - 1;
    const int mZero  = output_order - 1;

    s.assign(static_cast<std::size_t>(nTheta) * static_cast<std::size_t>(mCount), 0.0);
    const double phiScale = 2.0 * sylph::pi / static_cast<double>(nPhi);
    const double sqrt2 = std::numbers::sqrt2_v<double>;

    for (int i = 0; i < nTheta; ++i) {
        for (int m = -output_order + 1; m <= output_order - 1; ++m) {
            double sum = 0.0;
            for (int j = 0; j < nPhi; ++j) {
                const double phi = grid.phi(j);
                double r = 1.0;
                if (m > 0) {
                    r = sqrt2 * std::cos(m * phi);
                } else if (m < 0) {
                    r = sqrt2 * std::sin(std::abs(m) * phi);
                }
                sum += grid(i, j) * r;
            }
            s[static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount)
              + static_cast<std::size_t>(m + mZero)] = phiScale * sum;
        }
    }
}

VerificationExample::VerificationExample(int width, int height)
{
    resize(width, height);
}

void VerificationExample::process_event(const SDL_Event&, SDL_Window*)
{
}

void VerificationExample::update(float)
{
    if (!run_once_) {
        run_verification();
        run_once_ = true;
    }
}

void VerificationExample::resize(int width, int height)
{
    if (height <= 0) return;
    camera_.set_aspect_ratio(static_cast<float>(width) / static_cast<float>(height));
}

void VerificationExample::render(render::Renderer&)
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void VerificationExample::run_verification()
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    sylph::SHGridTransform grid_transform;

    // ========================================================================
    // Table 1: End-to-End SH Product (Naive vs Orig Grid vs #4 FFT vs #4+#2)
    // ========================================================================
    std::cout << "\n================================================================================================================\n";
    std::cout << "         SH Product Mathematical Verification & Speedup Analysis (Optimizations #4 + #2)\n";
    std::cout << "================================================================================================================\n";
    std::cout << "Orig Grid = tgamma + Direct DFT | #4 = tgamma + FFT | #4+#2 = Normalized Recurrence + Parity + FFT\n\n";

    std::cout << std::left
              << std::setw(7)  << "Order"
              << std::setw(13) << "Naive O(N^5)"
              << std::setw(13) << "Orig Grid"
              << std::setw(13) << "#4 FFT Only"
              << std::setw(13) << "#4+#2 Fast"
              << std::setw(11) << "vs Orig"
              << std::setw(11) << "vs #4"
              << std::setw(11) << "vs Naive"
              << std::setw(15) << "Err(vs Naive)"
              << "\n";
    std::cout << "----------------------------------------------------------------------------------------------------------------\n" << std::flush;

    std::vector<int> orders = {4, 6, 8, 10, 12, 16};

    for (int order : orders) {
        sylph::SHCoefficients a(order);
        sylph::SHCoefficients b(order);
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                a(l, m) = dist(gen);
                b(l, m) = dist(gen);
            }
        }

        // Warm up FFT & Legendre recurrence caches
        (void)grid_transform.product(a, b, order);

        auto t0 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_naive = naive_product(a, b, order);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms_naive = std::chrono::duration<double, std::milli>(t1 - t0).count();

        auto t2 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_orig = grid_transform.product_direct(a, b, order);
        auto t3 = std::chrono::high_resolution_clock::now();
        double ms_orig = std::chrono::duration<double, std::milli>(t3 - t2).count();
        (void)res_orig;

        auto t4 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_fft_only = grid_transform.product_fft_legacy_polar(a, b, order);
        auto t5 = std::chrono::high_resolution_clock::now();
        double ms_fft_only = std::chrono::duration<double, std::milli>(t5 - t4).count();
        (void)res_fft_only;

        const int iters_fast = 20;
        auto t6 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_fast(order);
        for (int it = 0; it < iters_fast; ++it) {
            res_fast = grid_transform.product(a, b, order);
        }
        auto t7 = std::chrono::high_resolution_clock::now();
        double ms_fast = std::chrono::duration<double, std::milli>(t7 - t6).count() / iters_fast;

        double max_err_naive = 0.0;
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                max_err_naive = std::max(max_err_naive, std::abs(res_naive(l, m) - res_fast(l, m)));
            }
        }

        std::cout << std::left
                  << std::setw(7) << order
                  << std::fixed << std::setprecision(4)
                  << std::setw(10) << ms_naive    << "ms "
                  << std::setw(10) << ms_orig     << "ms "
                  << std::setw(10) << ms_fft_only << "ms "
                  << std::setw(10) << ms_fast     << "ms "
                  << std::setprecision(1)
                  << std::setw(8)  << (ms_orig / ms_fast)     << "x  "
                  << std::setw(8)  << (ms_fft_only / ms_fast) << "x  "
                  << std::setw(8)  << (ms_naive / ms_fast)    << "x  "
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << max_err_naive
                  << "\n" << std::flush;
    }

    // ========================================================================
    // Table 2: Isolated Azimuthal Projection (Eq. 37) - Direct vs FFT Scaling
    // ========================================================================
    std::cout << "\n================================================================================================================\n";
    std::cout << "      Optimization #4: Isolated Azimuthal Projection (Eq. 37) - Direct vs FFT Scaling\n";
    std::cout << "================================================================================================================\n";
    std::cout << "Tests both exact nPhi = 3N-2 (Bluestein) and padded nPhi = next_pow2(3N-2) (Radix-2)\n\n";

    std::cout << std::left
              << std::setw(7)  << "Order"
              << std::setw(12) << "nPhi(exact)"
              << std::setw(12) << "nPhi(pow2)"
              << std::setw(14) << "Direct"
              << std::setw(14) << "FFT(Bluest.)"
              << std::setw(14) << "FFT(Radix2)"
              << std::setw(13) << "Speedup(R2)"
              << std::setw(14) << "Max Error"
              << "\n";
    std::cout << "----------------------------------------------------------------------------------------------------------------\n" << std::flush;

    std::vector<int> az_orders = {4, 8, 12, 16, 24, 32, 48, 64};

    for (int order : az_orders) {
        const int full_order = 2 * order - 1;
        const int tt = full_order + order - 1;
        const int nTheta = (tt + 1) / 2;
        const int nPhi_exact = tt;
        const int nPhi_pow2  = sylph::AzimuthalFFT::next_pow2(tt);

        sylph::gauss_legendre gl(nTheta);
        sylph::Spherical_grid grid_exact(nTheta, nPhi_exact, gl.nodes());
        sylph::Spherical_grid grid_pow2(nTheta, nPhi_pow2, gl.nodes());

        for (int i = 0; i < nTheta; ++i) {
            for (int j = 0; j < nPhi_exact; ++j) grid_exact(i, j) = dist(gen);
            for (int j = 0; j < nPhi_pow2;  ++j) grid_pow2(i, j)  = dist(gen);
        }

        std::vector<double> s_direct, s_bluestein, s_radix2, s_direct_pow2;
        sylph::AzimuthalFFT::project_rows(grid_exact, order, s_bluestein);
        sylph::AzimuthalFFT::project_rows(grid_pow2, order, s_radix2);

        const int iters = (order <= 16) ? 20 : 5;

        auto td0 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) direct_azimuthal_project(grid_exact, order, s_direct);
        auto td1 = std::chrono::high_resolution_clock::now();
        double ms_direct = std::chrono::duration<double, std::milli>(td1 - td0).count() / iters;

        auto tb0 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) sylph::AzimuthalFFT::project_rows(grid_exact, order, s_bluestein);
        auto tb1 = std::chrono::high_resolution_clock::now();
        double ms_bluestein = std::chrono::duration<double, std::milli>(tb1 - tb0).count() / iters;

        auto tr0 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) sylph::AzimuthalFFT::project_rows(grid_pow2, order, s_radix2);
        auto tr1 = std::chrono::high_resolution_clock::now();
        double ms_radix2 = std::chrono::duration<double, std::milli>(tr1 - tr0).count() / iters;

        direct_azimuthal_project(grid_pow2, order, s_direct_pow2);

        double max_err = 0.0;
        for (std::size_t idx = 0; idx < s_direct.size(); ++idx) {
            max_err = std::max(max_err, std::abs(s_direct[idx] - s_bluestein[idx]));
            max_err = std::max(max_err, std::abs(s_direct_pow2[idx] - s_radix2[idx]));
        }

        std::cout << std::left
                  << std::setw(7)  << order
                  << std::setw(12) << nPhi_exact
                  << std::setw(12) << nPhi_pow2
                  << std::fixed << std::setprecision(4)
                  << std::setw(10) << ms_direct    << " ms "
                  << std::setw(10) << ms_bluestein << " ms "
                  << std::setw(10) << ms_radix2    << " ms "
                  << std::setprecision(2)
                  << std::setw(10) << (ms_direct / ms_radix2) << "x  "
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << max_err
                  << "\n" << std::flush;
    }

    // ========================================================================
    // Table 3: High-Order Scaling & Numerical Stability (Order 16..64)
    // ========================================================================
    std::cout << "\n================================================================================================================\n";
    std::cout << "      Optimization #2: High-Order Scaling & Numerical Stability (Normalized Recurrence + Parity)\n";
    std::cout << "================================================================================================================\n";
    std::cout << std::left
              << std::setw(8)  << "Order"
              << std::setw(15) << "Product Time"
              << std::setw(18) << "Round-Trip Err"
              << std::setw(18) << "Commutativity"
              << std::setw(18) << "Unit-1 Product"
              << "\n";
    std::cout << "----------------------------------------------------------------------------------------------------------------\n" << std::flush;

    for (int order : {16, 24, 32, 48, 64}) {
        sylph::SHCoefficients a(order), b(order), unit(order);
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                a(l, m) = dist(gen);
                b(l, m) = dist(gen);
            }
        }
        unit(0, 0) = std::sqrt(4.0 * sylph::pi);

        (void)grid_transform.product(a, b, order);

        const int iters = (order <= 32) ? 10 : 3;
        auto t0 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients ab(order);
        for (int it = 0; it < iters; ++it) ab = grid_transform.product(a, b, order);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms_prod = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;

        auto ba = grid_transform.product(b, a, order);
        auto a1 = grid_transform.product(a, unit, order);

        const int tt = 2 * order - 1;
        const int nTheta = (tt + 1) / 2;
        const int nPhi = sylph::AzimuthalFFT::next_pow2(tt);
        sylph::gauss_legendre gl(nTheta);
        sylph::Spherical_grid g(nTheta, nPhi, gl.nodes());
        auto g_fwd = grid_transform.forward_separable(a, g);
        auto a_rt  = grid_transform.inverse(g_fwd, gl.weights(), order, order);

        double err_rt = 0.0, err_comm = 0.0, err_unit = 0.0;
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                err_rt   = std::max(err_rt,   std::abs(a(l, m) - a_rt(l, m)));
                err_comm = std::max(err_comm, std::abs(ab(l, m) - ba(l, m)));
                err_unit = std::max(err_unit, std::abs(a(l, m) - a1(l, m)));
            }
        }

        std::cout << std::left
                  << std::setw(8)  << order
                  << std::fixed << std::setprecision(4)
                  << std::setw(11) << ms_prod << " ms "
                  << std::scientific << std::setprecision(3)
                  << std::setw(18) << err_rt
                  << std::setw(18) << err_comm
                  << std::setw(18) << err_unit
                  << "\n" << std::flush;
    }

    std::cout << "================================================================================================================\n";

    // ========================================================================
    // Table 4: Optimization #9 - Zonal Harmonic (ZH) Fast Path
    // ========================================================================
    std::cout << "\n================================================================================================================\n";
    std::cout << "      Optimization #9: Zonal Harmonic (ZH) Fast Path - General (#4+#2) vs Fused ZH x SH & O(N^2) ZH x ZH\n";
    std::cout << "================================================================================================================\n";
    std::cout << std::left
              << std::setw(7)  << "Order"
              << std::setw(14) << "Gen(#4+#2)"
              << std::setw(14) << "ZH x SH"
              << std::setw(14) << "ZH x ZH"
              << std::setw(12) << "ZH*SH Gain"
              << std::setw(12) << "ZH*ZH Gain"
              << std::setw(15) << "Err(ZH*SH)"
              << std::setw(15) << "Err(ZH*ZH)"
              << "\n";
    std::cout << "----------------------------------------------------------------------------------------------------------------\n" << std::flush;

    for (int order : {4, 8, 12, 16, 24, 32, 48, 64}) {
        sylph::SHCoefficients za(order), zb(order), f(order);
        for (int l = 0; l < order; ++l) {
            za(l, 0) = dist(gen);
            zb(l, 0) = dist(gen);
            for (int m = -l; m <= l; ++m) f(l, m) = dist(gen);
        }

        (void)grid_transform.product_general(za, f, order);
        (void)grid_transform.product(za, f, order);
        (void)grid_transform.product(za, zb, order);

        const int iters = (order <= 16) ? 20 : ((order <= 32) ? 8 : 3);

        auto t0 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_gen_sh(order);
        for (int it = 0; it < iters; ++it) res_gen_sh = grid_transform.product_general(za, f, order);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms_gen = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;

        auto t2 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_zh_sh(order);
        for (int it = 0; it < iters; ++it) res_zh_sh = grid_transform.product(za, f, order);
        auto t3 = std::chrono::high_resolution_clock::now();
        double ms_zh_sh = std::chrono::duration<double, std::milli>(t3 - t2).count() / iters;

        auto t4 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_zh_zh(order);
        for (int it = 0; it < iters; ++it) res_zh_zh = grid_transform.product(za, zb, order);
        auto t5 = std::chrono::high_resolution_clock::now();
        double ms_zh_zh = std::chrono::duration<double, std::milli>(t5 - t4).count() / iters;

        auto res_gen_zh = grid_transform.product_general(za, zb, order);

        double err_zh_sh = 0.0, err_zh_zh = 0.0;
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                err_zh_sh = std::max(err_zh_sh, std::abs(res_gen_sh(l, m) - res_zh_sh(l, m)));
                err_zh_zh = std::max(err_zh_zh, std::abs(res_gen_zh(l, m) - res_zh_zh(l, m)));
            }
        }

        std::cout << std::left
                  << std::setw(7)  << order
                  << std::fixed << std::setprecision(4)
                  << std::setw(10) << ms_gen   << " ms "
                  << std::setw(10) << ms_zh_sh << " ms "
                  << std::setw(10) << ms_zh_zh << " ms "
                  << std::setprecision(2)
                  << std::setw(9)  << (ms_gen / ms_zh_sh) << "x  "
                  << std::setw(9)  << (ms_gen / ms_zh_zh) << "x  "
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << err_zh_sh
                  << std::setw(15) << err_zh_zh
                  << "\n" << std::flush;
    }

    std::cout << "================================================================================================================\n";

    // ========================================================================
    // Table 5: Optimization #8 - Cached GL + Cached Polar Basis + Zero-Alloc
    // ========================================================================
    std::cout << "\n================================================================================================================\n";
    std::cout << "      Optimization #8: Uncached (#4+#2) vs Cached GL + Polar Basis + Zero-Alloc Scratch\n";
    std::cout << "================================================================================================================\n";
    std::cout << std::left
              << std::setw(7)  << "Order"
              << std::setw(14) << "Before #8"
              << std::setw(14) << "After #8"
              << std::setw(11) << "Gain(#8)"
              << std::setw(14) << "ZH*SH(#8)"
              << std::setw(14) << "ZH*ZH(#8)"
              << std::setw(15) << "Err(SH*SH)"
              << "\n";
    std::cout << "----------------------------------------------------------------------------------------------------------------\n" << std::flush;

    for (int order : {4, 8, 12, 16, 24, 32, 48, 64}) {
        sylph::SHCoefficients a(order), b(order), za(order), zb(order);
        for (int l = 0; l < order; ++l) {
            za(l, 0) = dist(gen);
            zb(l, 0) = dist(gen);
            for (int m = -l; m <= l; ++m) { a(l, m) = dist(gen); b(l, m) = dist(gen); }
        }

        (void)grid_transform.product(a, b, order);
        (void)grid_transform.product(za, b, order);
        (void)grid_transform.product(za, zb, order);

        const int iters = (order <= 16) ? 30 : ((order <= 32) ? 10 : 4);

        auto t0 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_uncached(order);
        for (int it = 0; it < iters; ++it) res_uncached = grid_transform.product_uncached_before8(a, b, order);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms_before8 = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;

        auto t2 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_cached(order);
        for (int it = 0; it < iters; ++it) res_cached = grid_transform.product(a, b, order);
        auto t3 = std::chrono::high_resolution_clock::now();
        double ms_after8 = std::chrono::duration<double, std::milli>(t3 - t2).count() / iters;

        auto t4 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) (void)grid_transform.product(za, b, order);
        auto t5 = std::chrono::high_resolution_clock::now();
        double ms_zh_sh = std::chrono::duration<double, std::milli>(t5 - t4).count() / iters;

        auto t6 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) (void)grid_transform.product(za, zb, order);
        auto t7 = std::chrono::high_resolution_clock::now();
        double ms_zh_zh = std::chrono::duration<double, std::milli>(t7 - t6).count() / iters;

        double max_err = 0.0;
        for (int l = 0; l < order; ++l)
            for (int m = -l; m <= l; ++m)
                max_err = std::max(max_err, std::abs(res_uncached(l, m) - res_cached(l, m)));

        std::cout << std::left
                  << std::setw(7)  << order
                  << std::fixed << std::setprecision(4)
                  << std::setw(10) << ms_before8 << " ms "
                  << std::setw(10) << ms_after8  << " ms "
                  << std::setprecision(2)
                  << std::setw(8)  << (ms_before8 / ms_after8) << "x  "
                  << std::setprecision(4)
                  << std::setw(10) << ms_zh_sh   << " ms "
                  << std::setw(10) << ms_zh_zh   << " ms "
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << max_err
                  << "\n" << std::flush;
    }

    std::cout << "================================================================================================================\n\n" << std::flush;
}

} // namespace examples
