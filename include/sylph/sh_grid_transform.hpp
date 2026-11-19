/**
 * @file sh_grid_transform.hpp
 * @brief The True Paper Algorithm (SH to Spatial Grid Transform).
 * 
 * Transforms SH coefficients back into the spatial domain precisely at the 
 * Gauss-Legendre quadrature nodes. This allows complex SH products (O(N^5)) 
 * to be computed as simple point-wise multiplications (O(N^3)).
 */
#pragma once

#include "sylph/sh_coefficients.hpp"
#include "sylph/zonal_harmonics.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/spherical_grid.hpp"
#include "sylph/legendre.hpp"
#include "sylph/gauss_legendre.hpp"
#include "sylph/polar_basis_cache.hpp"
#include "sylph/azimuthal_fft.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include <types/constants.hpp>
namespace sylph{

class SHGridTransform{
public:


// equation 19
Spherical_grid  forward(
   const  SHCoefficients& coefficients,
    Spherical_grid grid
) const
{
    for (int i = 0; i < grid.theta_size(); ++i) {
        for (int j = 0; j < grid.phi_size(); ++j) {

            grid(i, j) =
                sh_.reconstruct(
                    coefficients,
                    grid.theta(i),
                    grid.phi(j)
                );
        }
    }

    return grid;
}


//eq 20 -25 (Normalized 3-term Legendre Recurrence + Equator Parity + FFT Azimuthal Synthesis)
Spherical_grid forward_separable(
    const SHCoefficients& coefficients,
    Spherical_grid grid
) const {
    return forward_separable_impl(coefficients, std::move(grid), true, true);
}

//eq 20 -25 (Legacy tgamma Polar + Direct O(N_theta * N_phi * M) Azimuthal reference)
Spherical_grid forward_separable_direct(
    const SHCoefficients& coefficients,
    Spherical_grid grid
) const {
    return forward_separable_impl(coefficients, std::move(grid), false, false);
}

    SHCoefficients inverse(
        const Spherical_grid& grid,
        const std::vector<double>& weights,
        int input_order,
        int output_order
    ) const{
return inverse_impl(
        grid,
        weights,
        input_order,
        output_order,
        true,
        true,
        true
    );
}

    SHCoefficients inverse_direct(
        const Spherical_grid& grid,
        const std::vector<double>& weights,
        int input_order,
        int output_order
    ) const{
return inverse_impl(
        grid,
        weights,
        input_order,
        output_order,
        true,
        false,
        false
    );
}

Spherical_grid product_grid(
    std::span<const SHCoefficients* const> inputs,
    const Spherical_grid& grid
) const
{
    if (inputs.size() < 2)
    {
        throw std::invalid_argument(
            "product_grid requires at least two inputs"
        );
    }

    for (const SHCoefficients* input : inputs)
    {
        if (input == nullptr)
        {
            throw std::invalid_argument(
                "product_grid received null input"
            );
        }

        if (input->order() <= 0)
        {
            throw std::invalid_argument(
                "all input orders must be positive"
            );
        }
    }

    Spherical_grid result =
        forward_separable(*inputs[0], grid);


    for (std::size_t i = 1; i < inputs.size(); ++i)
    {
        const Spherical_grid next =
            forward_separable(*inputs[i], grid);

        result.multiply_inplace(next);
    }

    return result;
}

SHCoefficients product(
    std::span<const SHCoefficients* const> inputs,
    int output_order
) const{
    
        if (inputs.size() < 2)
        {
            throw std::invalid_argument(
                "product requires at least "
                "two inputs"
            );
        }



        int full_order = 1;

        for (const SHCoefficients* input :
             inputs)
        {
            if (input == nullptr)
            {
                throw std::invalid_argument(
                    "product received null input"
                );
            }

            if (input->order() <= 0)
            {
                throw std::invalid_argument(
                    "all input orders must be positive"
                );
            }

            full_order +=
                input->order() - 1;
        }

    
        if (output_order <= 0 ||
            output_order > full_order)
        {
            throw std::invalid_argument(
                "output_order must satisfy "
                "1 <= output_order <= full product order"
            );
        }
        const int tt=full_order + output_order-1;
        const int nTheta = (tt + 1) / 2;

        // ------------------------------------------------------------
        // Optimization #9 + #8: Zero-alloc ZH dispatch + Cached Tables
        // ------------------------------------------------------------
        thread_local std::vector<const SHCoefficients*> zonal_inputs;
        thread_local std::vector<const SHCoefficients*> general_inputs;
        zonal_inputs.clear();
        general_inputs.clear();

        for (const SHCoefficients* input : inputs) {
            if (input->is_zonal()) {
                zonal_inputs.push_back(input);
            } else {
                general_inputs.push_back(input);
            }
        }

        // Case 1: All inputs are Zonal Harmonics -> pure O(N^2) cached 1D m=0 product
        if (general_inputs.empty()) {
            return product_zonal_all_impl(zonal_inputs, nTheta, output_order);
        }

        // Case 2: (K-1) Zonal inputs + 1 general SH input -> cached fused 1D-per-m polar product
        if (!zonal_inputs.empty() && general_inputs.size() == 1) {
            return product_zonal_sh_impl(zonal_inputs, *general_inputs[0], nTheta, output_order);
        }

        // Case 3 / General: >= 2 general SH inputs -> cached GL + cached polar + zero-alloc FFT
        const int nPhi = AzimuthalFFT::next_pow2(tt);
        return product_2d_cached(
            std::span<const SHCoefficients* const>(general_inputs),
            std::span<const SHCoefficients* const>(zonal_inputs),
            nTheta,
            nPhi,
            output_order
        );
}

/**
 * @brief Explicit ZH x SH fast product (cached GL + cached polar, zero FFT, zero 2D grid).
 */
SHCoefficients product(
    const ZHCoefficients& z,
    const SHCoefficients& f,
    int output_order
) const {
    const int full_order = z.order() + f.order() - 1;
    if (output_order <= 0 || output_order > full_order) {
        throw std::invalid_argument(
            "output_order must satisfy 1 <= output_order <= full product order"
        );
    }
    const int tt = full_order + output_order - 1;
    const int nTheta = (tt + 1) / 2;

    const auto& rule = gauss_legendre::get(nTheta);
    thread_local std::vector<double> z_prof;
    eval_zonal_profile_zh(z, nTheta, z_prof);

    if (f.is_zonal()) {
        thread_local std::vector<double> f_prof;
        eval_zonal_profile_sh(f, nTheta, f_prof);
        for (int i = 0; i < nTheta; ++i) {
            z_prof[static_cast<std::size_t>(i)] *= (2.0 * pi) * f_prof[static_cast<std::size_t>(i)];
        }
        return project_zonal_profile(z_prof, nTheta, rule.weights, output_order);
    }

    for (int i = 0; i < nTheta; ++i) {
        z_prof[static_cast<std::size_t>(i)] *= 2.0 * pi;
    }
    return product_profile_sh_fused(z_prof, f, nTheta, rule.weights, output_order);
}

SHCoefficients product(
    const SHCoefficients& f,
    const ZHCoefficients& z,
    int output_order
) const {
    return product(z, f, output_order);
}

/**
 * @brief Explicit ZH x ZH fast product (pure O(N^2) cached, returns ZHCoefficients).
 */
ZHCoefficients product(
    const ZHCoefficients& za,
    const ZHCoefficients& zb,
    int output_order
) const {
    const int full_order = za.order() + zb.order() - 1;
    if (output_order <= 0 || output_order > full_order) {
        throw std::invalid_argument(
            "output_order must satisfy 1 <= output_order <= full product order"
        );
    }
    const int tt = full_order + output_order - 1;
    const int nTheta = (tt + 1) / 2;

    const auto& rule = gauss_legendre::get(nTheta);
    thread_local std::vector<double> pa;
    thread_local std::vector<double> pb;
    eval_zonal_profile_zh(za, nTheta, pa);
    eval_zonal_profile_zh(zb, nTheta, pb);

    for (int i = 0; i < nTheta; ++i) {
        pa[static_cast<std::size_t>(i)] *= (2.0 * pi) * pb[static_cast<std::size_t>(i)];
    }

    const SHCoefficients sh_res = project_zonal_profile(pa, nTheta, rule.weights, output_order);
    return ZHCoefficients::from_sh(sh_res);
}

/**
 * @brief Benchmark helper: forces general 2D cached FFT + Polar path (bypasses #9 ZH fast path).
 */
SHCoefficients product_general(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order
) const {
    const int full_order = a.order() + b.order() - 1;
    if (output_order <= 0 || output_order > full_order) {
        throw std::invalid_argument(
            "output_order must satisfy 1 <= output_order <= full product order"
        );
    }
    const int tt = full_order + output_order - 1;
    const int nTheta = (tt + 1) / 2;
    const int nPhi   = AzimuthalFFT::next_pow2(tt);

    const std::array<const SHCoefficients*, 2> general_inputs = { &a, &b };
    return product_2d_cached(
        std::span<const SHCoefficients* const>(general_inputs),
        std::span<const SHCoefficients* const>{},
        nTheta,
        nPhi,
        output_order
    );
}

/**
 * @brief Benchmark helper: uncached GL Newton + on-the-fly recurrence (before #8).
 */
SHCoefficients product_uncached_before8(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order
) const {
    const int full_order = a.order() + b.order() - 1;
    if (output_order <= 0 || output_order > full_order) {
        throw std::invalid_argument(
            "output_order must satisfy 1 <= output_order <= full product order"
        );
    }
    const int tt = full_order + output_order - 1;
    const int nTheta = (tt + 1) / 2;
    const int nPhi   = AzimuthalFFT::next_pow2(tt);

    const auto rule = gauss_legendre::compute_uncached(nTheta);
    Spherical_grid grid(nTheta, nPhi, rule.nodes);

    const std::array<const SHCoefficients*, 2> inputs = { &a, &b };
    const Spherical_grid productGrid =
        product_grid(std::span<const SHCoefficients* const>(inputs), grid);

    return inverse(productGrid, rule.weights, full_order, output_order);
}

SHCoefficients product_direct(
    std::span<const SHCoefficients* const> inputs,
    int output_order
) const {
    if (inputs.size() < 2) {
        throw std::invalid_argument("product_direct requires at least two inputs");
    }

    int full_order = 1;
    for (const SHCoefficients* input : inputs) {
        if (input == nullptr) {
            throw std::invalid_argument("product_direct received null input");
        }
        if (input->order() <= 0) {
            throw std::invalid_argument("all input orders must be positive");
        }
        full_order += input->order() - 1;
    }

    if (output_order <= 0 || output_order > full_order) {
        throw std::invalid_argument(
            "output_order must satisfy 1 <= output_order <= full product order"
        );
    }

    const int tt = full_order + output_order - 1;
    const int nTheta = (tt + 1) / 2;
    const int nPhi   = tt;

    gauss_legendre gl(nTheta);
    Spherical_grid grid(nTheta, nPhi, gl.nodes());

    Spherical_grid productGrid = forward_separable_direct(*inputs[0], grid);
    for (std::size_t i = 1; i < inputs.size(); ++i) {
        const Spherical_grid next = forward_separable_direct(*inputs[i], grid);
        productGrid.multiply_inplace(next);
    }

    return inverse_direct(productGrid, gl.weights(), full_order, output_order);
}

SHCoefficients product_direct(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order
) const {
    const std::array<const SHCoefficients*, 2> inputs = { &a, &b };
    return product_direct(std::span<const SHCoefficients* const>(inputs), output_order);
}

/**
 * @brief Benchmark helper: FFT Azimuthal (#4) + legacy tgamma Polar (before #2).
 */
SHCoefficients product_fft_legacy_polar(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order
) const {
    const int full_order = a.order() + b.order() - 1;
    const int tt = full_order + output_order - 1;
    const int nTheta = (tt + 1) / 2;
    const int nPhi   = AzimuthalFFT::next_pow2(tt);

    gauss_legendre gl(nTheta);
    Spherical_grid grid(nTheta, nPhi, gl.nodes());

    Spherical_grid productGrid = forward_separable_impl(a, grid, true, false);
    const Spherical_grid next  = forward_separable_impl(b, grid, true, false);
    productGrid.multiply_inplace(next);

    return inverse_impl(productGrid, gl.weights(), full_order, output_order, true, true, false);
}

SHCoefficients product(
    std::span<const SHCoefficients* const> inputs,
    int output_order,
    int t
) const
{
    if (inputs.size() < 2)
    {
        throw std::invalid_argument(
            "product requires at least two inputs"
        );
    }

    if (t <= 0)
    {
        throw std::invalid_argument(
            "t must be positive"
        );
    }



    int full_order = 1;

    for (const SHCoefficients* input : inputs)
    {
        if (input == nullptr)
        {
            throw std::invalid_argument(
                "product received null input"
            );
        }

        if (input->order() <= 0)
        {
            throw std::invalid_argument(
                "all input orders must be positive"
            );
        }

        full_order += input->order() - 1;
    }

    if (output_order <= 0 ||
        output_order > full_order)
    {
        throw std::invalid_argument(
            "output_order must satisfy "
            "1 <= output_order <= full product order"
        );
    }



    const int nTheta = (t + 1) / 2;
    const int nPhi   = t;

    return product_2d_cached(
        inputs,
        std::span<const SHCoefficients* const>{},
        nTheta,
        nPhi,
        output_order
    );
}

double product_integral(
    std::span<const SHCoefficients* const> inputs
) const
{
    if (inputs.size() < 2)
    {
        throw std::invalid_argument(
            "product_integral requires at least "
            "two inputs"
        );
    }

    const SHCoefficients result =
        product(
            inputs,
            1
        );

    return
        2.0 *
        std::sqrt(pi) *
        result(0, 0);
}

double product_integral(
    const SHCoefficients& a,
    const SHCoefficients& b
) const
{
    const std::array<
        const SHCoefficients*,
        2
    > inputs = {
        &a,
        &b
    };

    return product_integral(
        std::span<
            const SHCoefficients* const
        >(inputs)
    );
}

 SHCoefficients product(
        const SHCoefficients& a,
        const SHCoefficients& b,
        int output_order
    ) const
    {
        const std::array<
            const SHCoefficients*,
            2
        > inputs = {
            &a,
            &b
        };

        return product(
            std::span<
                const SHCoefficients* const
            >(inputs),
            output_order
        );
    }
SHCoefficients product(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order,
    int t
) const
{
    const std::array<const SHCoefficients*, 2> inputs = {
        &a, &b
    };

    return product(
        std::span<const SHCoefficients* const>(inputs),
        output_order,
        t
    );
}

private:
SphericalHarmonics sh_;
Legendre  L_;

 double C_legacy(
        int l,
        int m,
        double theta
    )
  const  {
        const int am = std::abs(m);
        return L_.normalized_tgamma_legacy(l, am, std::cos(theta));
    }

    double R(
        int m,
        double phi
    )
   const  {
        if (m == 0) {
            return 1.0;
        }

        if (m > 0) {
            return std::sqrt(2.0) *
                   std::cos(m * phi);
        }

        return std::sqrt(2.0) *
               std::sin(std::abs(m) * phi);
    }

// -----------------------------------------------------------------------------
// Optimization #8 + #9: Cached Gauss-Legendre, Cached Polar Basis & ZH Helpers
// -----------------------------------------------------------------------------

/**
 * @brief Cached Forward Polar Synthesis (uses PolarBasisCache + Equator Parity).
 *        Zero recurrence, zero trig, zero heap allocation (reuses s).
 */
void forward_polar_cached(
    const SHCoefficients& coefficients,
    int nTheta,
    std::vector<double>& s
) const {
    const int order  = coefficients.order();
    const int mCount = 2 * order - 1;
    const int mZero  = order - 1;

    s.resize(static_cast<std::size_t>(nTheta) * static_cast<std::size_t>(mCount));

    const auto& ptbl = PolarBasisCache::get(nTheta, order);
    const double* c  = coefficients.data();

    auto eval_pair = [&](int i, int ip, bool paired) {
        double* row_n = s.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount);
        double* row_s = paired
            ? (s.data() + static_cast<std::size_t>(ip) * static_cast<std::size_t>(mCount))
            : nullptr;

        // m = 0
        {
            const double* p_col = ptbl.col(i, 0);
            double even_p = 0.0, odd_p = 0.0;
            int l = 0;
            for (; l + 1 < order; l += 2) {
                even_p += p_col[l]     * c[l * l + l];
                odd_p  += p_col[l + 1] * c[(l + 1) * (l + 1) + (l + 1)];
            }
            if (l < order) {
                even_p += p_col[l] * c[l * l + l];
            }
            row_n[mZero] = even_p + odd_p;
            if (paired) {
                row_s[mZero] = even_p - odd_p;
            }
        }

        // m = 1 .. order - 1
        for (int m = 1; m < order; ++m) {
            const double* p_col = ptbl.col(i, m);
            const int len = order - m;
            double even_pos = 0.0, odd_pos = 0.0;
            double even_neg = 0.0, odd_neg = 0.0;

            int k = 0;
            for (; k + 1 < len; k += 2) {
                const int le = m + k;
                const int be = le * le + le;
                const double pe = p_col[k];
                even_pos += pe * c[be + m];
                even_neg += pe * c[be - m];

                const int lo = le + 1;
                const int bo = lo * lo + lo;
                const double po = p_col[k + 1];
                odd_pos += po * c[bo + m];
                odd_neg += po * c[bo - m];
            }
            if (k < len) {
                const int le = m + k;
                const int be = le * le + le;
                const double pe = p_col[k];
                even_pos += pe * c[be + m];
                even_neg += pe * c[be - m];
            }

            row_n[mZero + m] = even_pos + odd_pos;
            row_n[mZero - m] = even_neg + odd_neg;
            if (paired) {
                row_s[mZero + m] = even_pos - odd_pos;
                row_s[mZero - m] = even_neg - odd_neg;
            }
        }
    };

    const int half = nTheta / 2;
    for (int i = 0; i < half; ++i) {
        eval_pair(i, nTheta - 1 - i, true);
    }
    if ((nTheta & 1) != 0) {
        eval_pair(half, half, false);
    }
}

/**
 * @brief Cached Inverse Polar Projection (Eq. 30, PolarBasisCache + Equator Parity).
 */
void inverse_polar_cached(
    int nTheta,
    const std::vector<double>& weights,
    const std::vector<double>& s,
    int mCount,
    int mZero,
    int output_order,
    SHCoefficients& result
) const {
    const auto& ptbl = PolarBasisCache::get(nTheta, output_order);
    double* res = result.data();

    auto accum_pair = [&](int i, int ip, bool paired) {
        const double wi = weights[static_cast<std::size_t>(i)];
        const double* row_n = s.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount);
        const double* row_s = paired
            ? (s.data() + static_cast<std::size_t>(ip) * static_cast<std::size_t>(mCount))
            : nullptr;

        // m = 0
        {
            const double* p_col = ptbl.col(i, 0);
            const double sp = paired ? wi * (row_n[mZero] + row_s[mZero]) : wi * row_n[mZero];
            const double sm = paired ? wi * (row_n[mZero] - row_s[mZero]) : wi * row_n[mZero];

            int l = 0;
            for (; l + 1 < output_order; l += 2) {
                res[l * l + l]                   += p_col[l]     * sp;
                res[(l + 1) * (l + 1) + (l + 1)] += p_col[l + 1] * sm;
            }
            if (l < output_order) {
                res[l * l + l] += p_col[l] * sp;
            }
        }

        // m = 1 .. output_order - 1
        for (int m = 1; m < output_order; ++m) {
            const double* p_col = ptbl.col(i, m);
            const int len = output_order - m;

            const double sp_pos = paired ? wi * (row_n[mZero + m] + row_s[mZero + m]) : wi * row_n[mZero + m];
            const double sm_pos = paired ? wi * (row_n[mZero + m] - row_s[mZero + m]) : wi * row_n[mZero + m];
            const double sp_neg = paired ? wi * (row_n[mZero - m] + row_s[mZero - m]) : wi * row_n[mZero - m];
            const double sm_neg = paired ? wi * (row_n[mZero - m] - row_s[mZero - m]) : wi * row_n[mZero - m];

            int k = 0;
            for (; k + 1 < len; k += 2) {
                const int le = m + k;
                const int be = le * le + le;
                const double pe = p_col[k];
                res[be + m] += pe * sp_pos;
                res[be - m] += pe * sp_neg;

                const int lo = le + 1;
                const int bo = lo * lo + lo;
                const double po = p_col[k + 1];
                res[bo + m] += po * sm_pos;
                res[bo - m] += po * sm_neg;
            }
            if (k < len) {
                const int le = m + k;
                const int be = le * le + le;
                const double pe = p_col[k];
                res[be + m] += pe * sp_pos;
                res[be - m] += pe * sp_neg;
            }
        }
    };

    const int half = nTheta / 2;
    for (int i = 0; i < half; ++i) {
        accum_pair(i, nTheta - 1 - i, true);
    }
    if ((nTheta & 1) != 0) {
        accum_pair(half, half, false);
    }
}

/**
 * @brief Zero-allocation 2D cached SH product pipeline (Optimization #8).
 */
SHCoefficients product_2d_cached(
    std::span<const SHCoefficients* const> general_inputs,
    std::span<const SHCoefficients* const> zonal_inputs,
    int nTheta,
    int nPhi,
    int output_order
) const {
    const auto& rule = gauss_legendre::get(nTheta);
    const std::size_t grid_size =
        static_cast<std::size_t>(nTheta) * static_cast<std::size_t>(nPhi);

    thread_local std::vector<double> s_buf;
    thread_local std::vector<double> grid_prod;
    thread_local std::vector<double> grid_tmp;

    grid_prod.resize(grid_size);

    // First general input -> grid_prod
    forward_polar_cached(*general_inputs[0], nTheta, s_buf);
    AzimuthalFFT::synthesize_rows_raw(s_buf.data(), nTheta, nPhi, general_inputs[0]->order(), grid_prod.data());

    // Remaining general inputs -> grid_tmp, multiply in-place
    if (general_inputs.size() > 1) {
        grid_tmp.resize(grid_size);
        double* gp = grid_prod.data();
        const double* gt = grid_tmp.data();

        for (std::size_t idx = 1; idx < general_inputs.size(); ++idx) {
            forward_polar_cached(*general_inputs[idx], nTheta, s_buf);
            AzimuthalFFT::synthesize_rows_raw(s_buf.data(), nTheta, nPhi, general_inputs[idx]->order(), grid_tmp.data());
            for (std::size_t p = 0; p < grid_size; ++p) {
                gp[p] *= gt[p];
            }
        }
    }

    // Optional zonal inputs: evaluate 1D cached profile and scale rows
    if (!zonal_inputs.empty()) {
        thread_local std::vector<double> z_prod;
        thread_local std::vector<double> z_tmp;
        z_prod.assign(static_cast<std::size_t>(nTheta), 1.0);
        for (const SHCoefficients* z : zonal_inputs) {
            eval_zonal_profile_sh(*z, nTheta, z_tmp);
            for (int i = 0; i < nTheta; ++i) {
                z_prod[static_cast<std::size_t>(i)] *= z_tmp[static_cast<std::size_t>(i)];
            }
        }
        double* gp = grid_prod.data();
        for (int i = 0; i < nTheta; ++i) {
            const double zi = z_prod[static_cast<std::size_t>(i)];
            double* row = gp + static_cast<std::size_t>(i) * static_cast<std::size_t>(nPhi);
            for (int j = 0; j < nPhi; ++j) {
                row[j] *= zi;
            }
        }
    }

    // Inverse FFT azimuthal projection + cached polar projection
    AzimuthalFFT::project_rows_raw(grid_prod.data(), nTheta, nPhi, output_order, s_buf);

    SHCoefficients result(output_order);
    inverse_polar_cached(
        nTheta,
        rule.weights,
        s_buf,
        2 * output_order - 1,
        output_order - 1,
        output_order,
        result
    );
    return result;
}

template <typename Accessor>
void eval_zonal_profile_core(
    int order,
    Accessor&& get_z,
    int nTheta,
    std::vector<double>& out
) const {
    out.resize(static_cast<std::size_t>(nTheta));
    const auto& ptbl = PolarBasisCache::get(nTheta, order);

    auto eval_pair = [&](int i, int ip, bool paired) {
        const double* p_col = ptbl.col(i, 0);
        double even_s = 0.0, odd_s = 0.0;
        int l = 0;
        for (; l + 1 < order; l += 2) {
            even_s += p_col[l]     * get_z(l);
            odd_s  += p_col[l + 1] * get_z(l + 1);
        }
        if (l < order) {
            even_s += p_col[l] * get_z(l);
        }
        out[static_cast<std::size_t>(i)] = even_s + odd_s;
        if (paired) {
            out[static_cast<std::size_t>(ip)] = even_s - odd_s;
        }
    };

    const int half = nTheta / 2;
    for (int i = 0; i < half; ++i) {
        eval_pair(i, nTheta - 1 - i, true);
    }
    if ((nTheta & 1) != 0) {
        eval_pair(half, half, false);
    }
}

void eval_zonal_profile_sh(
    const SHCoefficients& z,
    int nTheta,
    std::vector<double>& out
) const {
    const double* c = z.data();
    eval_zonal_profile_core(
        z.order(),
        [c](int l) { return c[l * l + l]; },
        nTheta,
        out
    );
}

void eval_zonal_profile_zh(
    const ZHCoefficients& z,
    int nTheta,
    std::vector<double>& out
) const {
    const double* d = z.data();
    eval_zonal_profile_core(
        z.order(),
        [d](int l) { return d[l]; },
        nTheta,
        out
    );
}

/**
 * @brief Project a 1D azimuthal-invariant profile (already scaled by 2*pi)
 *        onto m = 0 SH bands in O(N_theta * output_order) = O(N^2) using PolarBasisCache.
 */
SHCoefficients project_zonal_profile(
    const std::vector<double>& profile,
    int nTheta,
    const std::vector<double>& weights,
    int output_order
) const {
    SHCoefficients result(output_order);
    double* res = result.data();

    const auto& ptbl = PolarBasisCache::get(nTheta, output_order);

    auto accum_pair = [&](int i, int ip, bool paired) {
        const double* p_col = ptbl.col(i, 0);
        const double wi   = weights[static_cast<std::size_t>(i)];
        const double pi_n = profile[static_cast<std::size_t>(i)];
        const double pi_s = paired ? profile[static_cast<std::size_t>(ip)] : 0.0;

        const double sp = paired ? wi * (pi_n + pi_s) : wi * pi_n;
        const double sm = paired ? wi * (pi_n - pi_s) : wi * pi_n;

        int l = 0;
        for (; l + 1 < output_order; l += 2) {
            res[l * l + l]                   += p_col[l]     * sp;
            res[(l + 1) * (l + 1) + (l + 1)] += p_col[l + 1] * sm;
        }
        if (l < output_order) {
            res[l * l + l] += p_col[l] * sp;
        }
    };

    const int half = nTheta / 2;
    for (int i = 0; i < half; ++i) {
        accum_pair(i, nTheta - 1 - i, true);
    }
    if ((nTheta & 1) != 0) {
        accum_pair(half, half, false);
    }

    return result;
}

/**
 * @brief Pure O(N^2) cached product when ALL inputs are Zonal Harmonics.
 */
SHCoefficients product_zonal_all_impl(
    const std::vector<const SHCoefficients*>& zonal_inputs,
    int nTheta,
    int output_order
) const {
    const auto& rule = gauss_legendre::get(nTheta);

    thread_local std::vector<double> prod;
    thread_local std::vector<double> tmp;
    prod.assign(static_cast<std::size_t>(nTheta), 2.0 * pi);

    for (const SHCoefficients* z : zonal_inputs) {
        eval_zonal_profile_sh(*z, nTheta, tmp);
        for (int i = 0; i < nTheta; ++i) {
            prod[static_cast<std::size_t>(i)] *= tmp[static_cast<std::size_t>(i)];
        }
    }

    return project_zonal_profile(prod, nTheta, rule.weights, output_order);
}

/**
 * @brief Cached fused single-pass polar product of a 1D zonal profile (scaled by 2*pi)
 *        with a general SH function f. Zero Newton, zero recurrence, zero FFT, zero allocs.
 */
SHCoefficients product_profile_sh_fused(
    const std::vector<double>& z_prof_2pi,
    const SHCoefficients& f,
    int nTheta,
    const std::vector<double>& weights,
    int output_order
) const {
    SHCoefficients result(output_order);
    double* res = result.data();
    const double* cf = f.data();

    const int f_order = f.order();
    const int max_l   = std::max(f_order, output_order);
    const int m_limit = std::min(f_order, output_order);

    const auto& ptbl = PolarBasisCache::get(nTheta, max_l);

    auto process_pair = [&](int i, int ip, bool paired) {
        const double wi = weights[static_cast<std::size_t>(i)];
        const double zn = z_prof_2pi[static_cast<std::size_t>(i)];
        const double zs = paired ? z_prof_2pi[static_cast<std::size_t>(ip)] : 0.0;

        for (int m = 0; m < m_limit; ++m) {
            const double* p_col = ptbl.col(i, m);
            const int f_len = f_order - m;
            const int o_len = output_order - m;

            if (m == 0) {
                double even_p = 0.0, odd_p = 0.0;
                int k = 0;
                for (; k + 1 < f_len; k += 2) {
                    even_p += p_col[k]     * cf[k * k + k];
                    odd_p  += p_col[k + 1] * cf[(k + 1) * (k + 1) + (k + 1)];
                }
                if (k < f_len) {
                    even_p += p_col[k] * cf[k * k + k];
                }

                double sp, sm;
                if (paired) {
                    const double tn = zn * (even_p + odd_p);
                    const double ts = zs * (even_p - odd_p);
                    sp = wi * (tn + ts);
                    sm = wi * (tn - ts);
                } else {
                    const double t = zn * (even_p + odd_p);
                    sp = sm = wi * t;
                }

                k = 0;
                for (; k + 1 < o_len; k += 2) {
                    res[k * k + k]                   += p_col[k]     * sp;
                    res[(k + 1) * (k + 1) + (k + 1)] += p_col[k + 1] * sm;
                }
                if (k < o_len) {
                    res[k * k + k] += p_col[k] * sp;
                }
            } else {
                double even_pos = 0.0, odd_pos = 0.0;
                double even_neg = 0.0, odd_neg = 0.0;
                int k = 0;
                for (; k + 1 < f_len; k += 2) {
                    const int le = m + k;
                    const int be = le * le + le;
                    const double pe = p_col[k];
                    even_pos += pe * cf[be + m];
                    even_neg += pe * cf[be - m];

                    const int lo = le + 1;
                    const int bo = lo * lo + lo;
                    const double po = p_col[k + 1];
                    odd_pos += po * cf[bo + m];
                    odd_neg += po * cf[bo - m];
                }
                if (k < f_len) {
                    const int le = m + k;
                    const int be = le * le + le;
                    const double pe = p_col[k];
                    even_pos += pe * cf[be + m];
                    even_neg += pe * cf[be - m];
                }

                double sp_pos, sm_pos, sp_neg, sm_neg;
                if (paired) {
                    const double tn_pos = zn * (even_pos + odd_pos);
                    const double ts_pos = zs * (even_pos - odd_pos);
                    sp_pos = wi * (tn_pos + ts_pos);
                    sm_pos = wi * (tn_pos - ts_pos);

                    const double tn_neg = zn * (even_neg + odd_neg);
                    const double ts_neg = zs * (even_neg - odd_neg);
                    sp_neg = wi * (tn_neg + ts_neg);
                    sm_neg = wi * (tn_neg - ts_neg);
                } else {
                    const double t_pos = zn * (even_pos + odd_pos);
                    const double t_neg = zn * (even_neg + odd_neg);
                    sp_pos = sm_pos = wi * t_pos;
                    sp_neg = sm_neg = wi * t_neg;
                }

                k = 0;
                for (; k + 1 < o_len; k += 2) {
                    const int le = m + k;
                    const int be = le * le + le;
                    const double pe = p_col[k];
                    res[be + m] += pe * sp_pos;
                    res[be - m] += pe * sp_neg;

                    const int lo = le + 1;
                    const int bo = lo * lo + lo;
                    const double po = p_col[k + 1];
                    res[bo + m] += po * sm_pos;
                    res[bo - m] += po * sm_neg;
                }
                if (k < o_len) {
                    const int le = m + k;
                    const int be = le * le + le;
                    const double pe = p_col[k];
                    res[be + m] += pe * sp_pos;
                    res[be - m] += pe * sp_neg;
                }
            }
        }
    };

    const int half = nTheta / 2;
    for (int i = 0; i < half; ++i) {
        process_pair(i, nTheta - 1 - i, true);
    }
    if ((nTheta & 1) != 0) {
        process_pair(half, half, false);
    }

    return result;
}

SHCoefficients product_zonal_sh_impl(
    const std::vector<const SHCoefficients*>& zonal_inputs,
    const SHCoefficients& general_input,
    int nTheta,
    int output_order
) const {
    const auto& rule = gauss_legendre::get(nTheta);

    thread_local std::vector<double> z_prod;
    thread_local std::vector<double> tmp;
    z_prod.assign(static_cast<std::size_t>(nTheta), 2.0 * pi);

    for (const SHCoefficients* z : zonal_inputs) {
        eval_zonal_profile_sh(*z, nTheta, tmp);
        for (int i = 0; i < nTheta; ++i) {
            z_prod[static_cast<std::size_t>(i)] *= tmp[static_cast<std::size_t>(i)];
        }
    }

    return product_profile_sh_fused(z_prod, general_input, nTheta, rule.weights, output_order);
}

/**
 * @brief Forward Polar Synthesis using directly normalized 3-term Legendre
 *        recurrence + simultaneous (+/-m) reuse + North-South Equator Parity Folding.
 *
 * For symmetric Gauss-Legendre pairs (i, nTheta - 1 - i) with x_{i'} = -x_i:
 *   P~_l^m(-x) = (-1)^(l - m) P~_l^m(x)
 * Splitting l = m..order-1 into even (E) and odd (O) parity sums yields:
 *   s[i,  m] = E + O
 *   s[i', m] = E - O
 * halving both Legendre evaluations and polar multiply-adds.
 */
void forward_polar_recurrence(
    const SHCoefficients& coefficients,
    const Spherical_grid& grid,
    std::vector<double>& s
) const {
    const int order  = coefficients.order();
    const int nTheta = grid.theta_size();
    const int mCount = 2 * order - 1;
    const int mZero  = order - 1;

    s.assign(
        static_cast<std::size_t>(nTheta) * static_cast<std::size_t>(mCount),
        0.0
    );

    const auto& rec = Legendre::get_recurrence_table(order);
    const double* diag = rec.diag.data();
    const double* a_tbl = rec.a.data();
    const double* b_tbl = rec.b.data();
    const std::size_t* col_off = rec.col_offset.data();
    const double* c = coefficients.data();
    const double p00 = Legendre::inv_sqrt_4pi();

    std::vector<double> cos_t(static_cast<std::size_t>(nTheta));
    std::vector<double> sin_t(static_cast<std::size_t>(nTheta));
    for (int i = 0; i < nTheta; ++i) {
        const double th = grid.theta(i);
        cos_t[static_cast<std::size_t>(i)] = std::cos(th);
        sin_t[static_cast<std::size_t>(i)] = std::max(0.0, std::sin(th));
    }

    const int half = nTheta / 2;
    bool symmetric = (nTheta > 1);
    for (int i = 0; i < half; ++i) {
        const int ip = nTheta - 1 - i;
        if (std::abs(cos_t[static_cast<std::size_t>(i)] + cos_t[static_cast<std::size_t>(ip)]) > 1e-12) {
            symmetric = false;
            break;
        }
    }

    auto eval_row = [&](int i, int ip, bool paired) {
        const double x  = cos_t[static_cast<std::size_t>(i)];
        const double st = sin_t[static_cast<std::size_t>(i)];

        double* row_n = s.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount);
        double* row_s = paired
            ? (s.data() + static_cast<std::size_t>(ip) * static_cast<std::size_t>(mCount))
            : nullptr;

        double pmm = p00;

        // m = 0
        {
            const double* a_col = a_tbl + col_off[0];
            const double* b_col = b_tbl + col_off[0];

            double even_p = pmm * c[0];
            double odd_p  = 0.0;

            if (order > 1) {
                double p0 = pmm;
                double p1 = a_col[1] * x * p0;
                odd_p += p1 * c[2]; // index(1, 0) = 2

                int l = 2;
                for (; l + 1 < order; l += 2) {
                    const double pe = a_col[l] * x * p1 - b_col[l] * p0;
                    even_p += pe * c[l * l + l];

                    const double po = a_col[l + 1] * x * pe - b_col[l + 1] * p1;
                    odd_p += po * c[(l + 1) * (l + 1) + (l + 1)];

                    p0 = pe;
                    p1 = po;
                }
                if (l < order) {
                    const double pe = a_col[l] * x * p1 - b_col[l] * p0;
                    even_p += pe * c[l * l + l];
                }
            }

            row_n[mZero] = even_p + odd_p;
            if (paired) {
                row_s[mZero] = even_p - odd_p;
            }
        }

        // m = 1 .. order - 1 (reusing the exact same P~_l^m for both +m and -m)
        for (int m = 1; m < order; ++m) {
            pmm = -diag[static_cast<std::size_t>(m)] * st * pmm;

            const double* a_col = a_tbl + col_off[static_cast<std::size_t>(m)];
            const double* b_col = b_tbl + col_off[static_cast<std::size_t>(m)];

            const int base0 = m * m + m;
            double even_pos = pmm * c[base0 + m];
            double even_neg = pmm * c[base0 - m];
            double odd_pos  = 0.0;
            double odd_neg  = 0.0;

            if (m + 1 < order) {
                double p0 = pmm;
                double p1 = a_col[1] * x * p0;
                const int base1 = (m + 1) * (m + 1) + (m + 1);
                odd_pos += p1 * c[base1 + m];
                odd_neg += p1 * c[base1 - m];

                int l = m + 2;
                for (; l + 1 < order; l += 2) {
                    const int k = l - m;
                    const double pe = a_col[k] * x * p1 - b_col[k] * p0;
                    const int be = l * l + l;
                    even_pos += pe * c[be + m];
                    even_neg += pe * c[be - m];

                    const double po = a_col[k + 1] * x * pe - b_col[k + 1] * p1;
                    const int bo = (l + 1) * (l + 1) + (l + 1);
                    odd_pos += po * c[bo + m];
                    odd_neg += po * c[bo - m];

                    p0 = pe;
                    p1 = po;
                }
                if (l < order) {
                    const int k = l - m;
                    const double pe = a_col[k] * x * p1 - b_col[k] * p0;
                    const int be = l * l + l;
                    even_pos += pe * c[be + m];
                    even_neg += pe * c[be - m];
                }
            }

            row_n[mZero + m] = even_pos + odd_pos;
            row_n[mZero - m] = even_neg + odd_neg;
            if (paired) {
                row_s[mZero + m] = even_pos - odd_pos;
                row_s[mZero - m] = even_neg - odd_neg;
            }
        }
    };

    if (symmetric) {
        for (int i = 0; i < half; ++i) {
            eval_row(i, nTheta - 1 - i, true);
        }
        if ((nTheta & 1) != 0) {
            eval_row(half, half, false);
        }
    } else {
        for (int i = 0; i < nTheta; ++i) {
            eval_row(i, i, false);
        }
    }
}

/**
 * @brief Inverse Polar Projection (Eq. 30) using directly normalized 3-term
 *        Legendre recurrence + (+/-m) reuse + North-South Equator Parity Folding.
 *
 * For symmetric Gauss-Legendre pairs (i, i') with x_{i'} = -x_i and w_{i'} = w_i:
 *   s^+ = w_i * (s[i, m] + s[i', m])   (accumulated into even l - |m|)
 *   s^- = w_i * (s[i, m] - s[i', m])   (accumulated into odd  l - |m|)
 */
void inverse_polar_recurrence(
    const Spherical_grid& grid,
    const std::vector<double>& weights,
    const std::vector<double>& s,
    int mCount,
    int mZero,
    int output_order,
    SHCoefficients& result
) const {
    const int nTheta = grid.theta_size();

    const auto& rec = Legendre::get_recurrence_table(output_order);
    const double* diag = rec.diag.data();
    const double* a_tbl = rec.a.data();
    const double* b_tbl = rec.b.data();
    const std::size_t* col_off = rec.col_offset.data();
    double* res = result.data();
    const double p00 = Legendre::inv_sqrt_4pi();

    std::vector<double> cos_t(static_cast<std::size_t>(nTheta));
    std::vector<double> sin_t(static_cast<std::size_t>(nTheta));
    for (int i = 0; i < nTheta; ++i) {
        const double th = grid.theta(i);
        cos_t[static_cast<std::size_t>(i)] = std::cos(th);
        sin_t[static_cast<std::size_t>(i)] = std::max(0.0, std::sin(th));
    }

    const int half = nTheta / 2;
    bool symmetric = (nTheta > 1);
    for (int i = 0; i < half; ++i) {
        const int ip = nTheta - 1 - i;
        if (std::abs(cos_t[static_cast<std::size_t>(i)] + cos_t[static_cast<std::size_t>(ip)]) > 1e-12 ||
            std::abs(weights[static_cast<std::size_t>(i)] - weights[static_cast<std::size_t>(ip)]) > 1e-12) {
            symmetric = false;
            break;
        }
    }

    auto accum_row = [&](int i, int ip, bool paired) {
        const double x  = cos_t[static_cast<std::size_t>(i)];
        const double st = sin_t[static_cast<std::size_t>(i)];
        const double wi = weights[static_cast<std::size_t>(i)];

        const double* row_n = s.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount);
        const double* row_s = paired
            ? (s.data() + static_cast<std::size_t>(ip) * static_cast<std::size_t>(mCount))
            : nullptr;

        double pmm = p00;

        // m = 0
        {
            const double* a_col = a_tbl + col_off[0];
            const double* b_col = b_tbl + col_off[0];

            const double sp = paired ? wi * (row_n[mZero] + row_s[mZero]) : wi * row_n[mZero];
            const double sm = paired ? wi * (row_n[mZero] - row_s[mZero]) : wi * row_n[mZero];

            res[0] += pmm * sp;

            if (output_order > 1) {
                double p0 = pmm;
                double p1 = a_col[1] * x * p0;
                res[2] += p1 * sm; // index(1, 0) = 2

                int l = 2;
                for (; l + 1 < output_order; l += 2) {
                    const double pe = a_col[l] * x * p1 - b_col[l] * p0;
                    res[l * l + l] += pe * sp;

                    const double po = a_col[l + 1] * x * pe - b_col[l + 1] * p1;
                    res[(l + 1) * (l + 1) + (l + 1)] += po * sm;

                    p0 = pe;
                    p1 = po;
                }
                if (l < output_order) {
                    const double pe = a_col[l] * x * p1 - b_col[l] * p0;
                    res[l * l + l] += pe * sp;
                }
            }
        }

        // m = 1 .. output_order - 1 (reusing P~_l^m for both +m and -m)
        for (int m = 1; m < output_order; ++m) {
            pmm = -diag[static_cast<std::size_t>(m)] * st * pmm;

            const double* a_col = a_tbl + col_off[static_cast<std::size_t>(m)];
            const double* b_col = b_tbl + col_off[static_cast<std::size_t>(m)];

            const double sp_pos = paired ? wi * (row_n[mZero + m] + row_s[mZero + m]) : wi * row_n[mZero + m];
            const double sm_pos = paired ? wi * (row_n[mZero + m] - row_s[mZero + m]) : wi * row_n[mZero + m];
            const double sp_neg = paired ? wi * (row_n[mZero - m] + row_s[mZero - m]) : wi * row_n[mZero - m];
            const double sm_neg = paired ? wi * (row_n[mZero - m] - row_s[mZero - m]) : wi * row_n[mZero - m];

            const int base0 = m * m + m;
            res[base0 + m] += pmm * sp_pos;
            res[base0 - m] += pmm * sp_neg;

            if (m + 1 < output_order) {
                double p0 = pmm;
                double p1 = a_col[1] * x * p0;
                const int base1 = (m + 1) * (m + 1) + (m + 1);
                res[base1 + m] += p1 * sm_pos;
                res[base1 - m] += p1 * sm_neg;

                int l = m + 2;
                for (; l + 1 < output_order; l += 2) {
                    const int k = l - m;
                    const double pe = a_col[k] * x * p1 - b_col[k] * p0;
                    const int be = l * l + l;
                    res[be + m] += pe * sp_pos;
                    res[be - m] += pe * sp_neg;

                    const double po = a_col[k + 1] * x * pe - b_col[k + 1] * p1;
                    const int bo = (l + 1) * (l + 1) + (l + 1);
                    res[bo + m] += po * sm_pos;
                    res[bo - m] += po * sm_neg;

                    p0 = pe;
                    p1 = po;
                }
                if (l < output_order) {
                    const int k = l - m;
                    const double pe = a_col[k] * x * p1 - b_col[k] * p0;
                    const int be = l * l + l;
                    res[be + m] += pe * sp_pos;
                    res[be - m] += pe * sp_neg;
                }
            }
        }
    };

    if (symmetric) {
        for (int i = 0; i < half; ++i) {
            accum_row(i, nTheta - 1 - i, true);
        }
        if ((nTheta & 1) != 0) {
            accum_row(half, half, false);
        }
    } else {
        for (int i = 0; i < nTheta; ++i) {
            accum_row(i, i, false);
        }
    }
}

Spherical_grid forward_separable_impl(
    const SHCoefficients& coefficients,
    Spherical_grid grid,
    bool use_fft,
    bool use_fast_polar = true
) const {
    const int order = coefficients.order();
    const int mCount = 2 * order - 1;
    std::vector<double> s;

    // Polar synthesis: sum over l for each (theta_i, m)
    if (use_fast_polar) {
        forward_polar_recurrence(coefficients, grid, s);
    } else {
        s.assign(
            static_cast<std::size_t>(grid.theta_size()) *
            static_cast<std::size_t>(mCount),
            0.0
        );
        auto mIndex = [order](int m) {
            return m + order - 1;
        };

        for (int i = 0; i < grid.theta_size(); ++i) {
            const double theta = grid.theta(i);

            for (int m = -order + 1; m <= order - 1; ++m) {
                double sum = 0.0;
                for (int l = std::abs(m); l < order; ++l) {
                    sum += coefficients(l, m) * C_legacy(l, m, theta);
                }
                s[static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount)
                  + static_cast<std::size_t>(mIndex(m))] = sum;
            }
        }
    }

    // Azimuthal synthesis: sum over m for each (theta_i, phi_j)
    if (use_fft) {
        AzimuthalFFT::synthesize_rows(s, order, grid);
    } else {
        auto mIndex = [order](int m) {
            return m + order - 1;
        };
        for (int i = 0; i < grid.theta_size(); ++i) {
            for (int j = 0; j < grid.phi_size(); ++j) {
                const double phi = grid.phi(j);
                double value = 0.0;
                for (int m = -order + 1; m <= order - 1; ++m) {
                    value += R(m, phi) *
                        s[static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount)
                          + static_cast<std::size_t>(mIndex(m))];
                }
                grid(i, j) = value;
            }
        }
    }

    return grid;
}

SHCoefficients inverse_impl(
    const Spherical_grid& grid,
    const std::vector<double>& weights,
    int input_order,
    int output_order,
    bool validate_resolution,
    bool use_fft = true,
    bool use_fast_polar = true
) const
{
    if (input_order <= 0)
    {
        throw std::invalid_argument(
            "input_order must be positive"
        );
    }

    if (output_order <= 0 ||
        output_order > input_order)
    {
        throw std::invalid_argument(
            "output_order must satisfy "
            "1 <= output_order <= input_order"
        );
    }

    if (weights.size() !=
        static_cast<std::size_t>(
            grid.theta_size()
        ))
    {
        throw std::invalid_argument(
            "weights size must equal "
            "grid theta size"
        );
    }

    // ------------------------------------------------------------
    // Exact-resolution validation is optional here.
    //
    // Public inverse() passes true.
    // Approximate product(..., t) passes false.
    // ------------------------------------------------------------

    if (validate_resolution)
    {
        const int min_theta =
            (input_order + output_order) / 2;

        const int min_phi =
            input_order + output_order - 1;

        if (grid.theta_size() < min_theta)
        {
            throw std::invalid_argument(
                "theta resolution is insufficient "
                "for inverse transform"
            );
        }

        if (grid.phi_size() < min_phi)
        {
            throw std::invalid_argument(
                "phi resolution is insufficient "
                "for inverse transform"
            );
        }
    }

    // ------------------------------------------------------------
    // Eq. 37: azimuthal projection
    // ------------------------------------------------------------

    std::vector<double> s;
    int mCount = 0;
    int mZero  = 0;

    if (use_fft)
    {
        mCount = 2 * output_order - 1;
        mZero  = output_order - 1;
        AzimuthalFFT::project_rows(grid, output_order, s);
    }
    else
    {
        mCount = 2 * input_order - 1;
        mZero  = input_order - 1;

        s.assign(
            static_cast<std::size_t>(grid.theta_size()) *
            static_cast<std::size_t>(mCount),
            0.0
        );

        const double phiOrder =
            2.0 * pi / static_cast<double>(grid.phi_size());

        for (int i = 0; i < grid.theta_size(); ++i)
        {
            for (int m = -input_order + 1; m <= input_order - 1; ++m)
            {
                double sum = 0.0;
                for (int j = 0; j < grid.phi_size(); ++j)
                {
                    sum += grid(i, j) * R(m, grid.phi(j));
                }

                s[static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount) +
                  static_cast<std::size_t>(m + mZero)] = phiOrder * sum;
            }
        }
    }

    // ------------------------------------------------------------
    // Eq. 30: polar projection
    // ------------------------------------------------------------

    SHCoefficients result(
        output_order
    );

    if (use_fast_polar)
    {
        inverse_polar_recurrence(
            grid,
            weights,
            s,
            mCount,
            mZero,
            output_order,
            result
        );
    }
    else
    {
        for (int l = 0;
             l < output_order;
             ++l)
        {
            for (int m = -l;
                 m <= l;
                 ++m)
            {
                double coefficient = 0.0;

                for (int i = 0;
                     i < grid.theta_size();
                     ++i)
                {
                    coefficient +=
                        weights[i] *
                        C_legacy(
                            l,
                            m,
                            grid.theta(i)
                        ) *
                        s[
                            static_cast<std::size_t>(i) *
                            static_cast<std::size_t>(mCount) +
                            static_cast<std::size_t>(m + mZero)
                        ];
                }

                result(l, m) =
                    coefficient;
            }
        }
    }

    return result;
}
};

}
