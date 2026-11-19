/**
 * @file zonal_harmonics.hpp
 * @brief Zonal Harmonic (ZH) container and O(N^2) operations.
 *
 * A Zonal Harmonic is rotationally symmetric around the z-axis (only m = 0
 * coefficients are non-zero), storing N coefficients instead of N^2.
 *
 * Provides:
 *  - Compact O(N) storage of z_l = c(l, 0)
 *  - Exact O(N^2) rotation to an arbitrary axis (theta_d, phi_d) via the
 *    Spherical Harmonic Addition Theorem:
 *      c(l, m) = sqrt(4*pi / (2l + 1)) * z_l * Y_l^m(theta_d, phi_d)
 *  - Exact O(N^2) Funk-Hecke spherical convolution with any SH function:
 *      (Z * F)(l, m) = sqrt(4*pi / (2l + 1)) * z_l * F(l, m)
 */
#pragma once

#include "sylph/legendre.hpp"
#include "sylph/sh_coefficients.hpp"
#include "types/constants.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace sylph {

class ZHCoefficients {
public:
    explicit ZHCoefficients(int order)
        : order_(order),
          data_(static_cast<std::size_t>(order > 0 ? order : 0), 0.0)
    {
        if (order <= 0) {
            throw std::invalid_argument("ZHCoefficients order must be positive");
        }
    }

    static ZHCoefficients from_sh(const SHCoefficients& sh) {
        ZHCoefficients zh(sh.order());
        for (int l = 0; l < sh.order(); ++l) {
            zh(l) = sh(l, 0);
        }
        return zh;
    }

    double& operator()(int l) {
        validate(l);
        return data_[static_cast<std::size_t>(l)];
    }

    double operator()(int l) const {
        validate(l);
        return data_[static_cast<std::size_t>(l)];
    }

    int order() const noexcept { return order_; }
    const double* data() const noexcept { return data_.data(); }
    double* data() noexcept { return data_.data(); }

    /**
     * @brief Expand into an axis-aligned (+z) SHCoefficients object.
     */
    SHCoefficients to_sh() const {
        SHCoefficients out(order_);
        for (int l = 0; l < order_; ++l) {
            out(l, 0) = data_[static_cast<std::size_t>(l)];
        }
        return out;
    }

    /**
     * @brief Exact O(N^2) rotation of this ZH lobe to an arbitrary axis (theta_d, phi_d)
     *        using the Spherical Harmonic Addition Theorem:
     *          c_{l, m} = sqrt(4*pi / (2l + 1)) * z_l * Y_l^m(theta_d, phi_d)
     */
    SHCoefficients rotate_to_sh(double theta_d, double phi_d) const {
        SHCoefficients out(order_);
        double* dst = out.data();

        const double x  = std::cos(theta_d);
        const double st = std::max(0.0, std::sin(theta_d));

        const auto& rec = Legendre::get_recurrence_table(order_);
        const double* diag = rec.diag.data();
        const double* a_tbl = rec.a.data();
        const double* b_tbl = rec.b.data();
        const std::size_t* col_off = rec.col_offset.data();
        const double p00 = Legendre::inv_sqrt_4pi();

        // Precompute band factors k_l = sqrt(4*pi / (2l + 1)) * z_l
        std::vector<double> band_scale(static_cast<std::size_t>(order_));
        for (int l = 0; l < order_; ++l) {
            band_scale[static_cast<std::size_t>(l)] =
                std::sqrt(4.0 * pi / static_cast<double>(2 * l + 1)) *
                data_[static_cast<std::size_t>(l)];
        }

        constexpr double sqrt2 = std::numbers::sqrt2_v<double>;
        double pmm = p00;

        // m = 0
        {
            const double* a_col = a_tbl + col_off[0];
            const double* b_col = b_tbl + col_off[0];

            dst[0] = band_scale[0] * pmm;
            if (order_ > 1) {
                double p0 = pmm;
                double p1 = a_col[1] * x * p0;
                dst[2] = band_scale[1] * p1;
                for (int l = 2; l < order_; ++l) {
                    const double p2 = a_col[l] * x * p1 - b_col[l] * p0;
                    dst[l * l + l] = band_scale[static_cast<std::size_t>(l)] * p2;
                    p0 = p1;
                    p1 = p2;
                }
            }
        }

        // m = 1 .. order_ - 1
        for (int m = 1; m < order_; ++m) {
            pmm = -diag[static_cast<std::size_t>(m)] * st * pmm;

            const double cosm = sqrt2 * std::cos(static_cast<double>(m) * phi_d);
            const double sinm = sqrt2 * std::sin(static_cast<double>(m) * phi_d);

            const double* a_col = a_tbl + col_off[static_cast<std::size_t>(m)];
            const double* b_col = b_tbl + col_off[static_cast<std::size_t>(m)];

            {
                const double s0 = band_scale[static_cast<std::size_t>(m)] * pmm;
                const int base0 = m * m + m;
                dst[base0 + m] = s0 * cosm;
                dst[base0 - m] = s0 * sinm;
            }

            if (m + 1 < order_) {
                double p0 = pmm;
                double p1 = a_col[1] * x * p0;
                {
                    const double s1 = band_scale[static_cast<std::size_t>(m + 1)] * p1;
                    const int base1 = (m + 1) * (m + 1) + (m + 1);
                    dst[base1 + m] = s1 * cosm;
                    dst[base1 - m] = s1 * sinm;
                }

                for (int l = m + 2; l < order_; ++l) {
                    const int k = l - m;
                    const double p2 = a_col[k] * x * p1 - b_col[k] * p0;
                    const double sl = band_scale[static_cast<std::size_t>(l)] * p2;
                    const int base = l * l + l;
                    dst[base + m] = sl * cosm;
                    dst[base - m] = sl * sinm;
                    p0 = p1;
                    p1 = p2;
                }
            }
        }

        return out;
    }

    /**
     * @brief Exact O(N^2) Funk-Hecke spherical convolution of this ZH kernel with SH function f.
     */
    SHCoefficients convolve(const SHCoefficients& f) const {
        const int out_order = std::min(order_, f.order());
        SHCoefficients out(out_order);
        const double* src = f.data();
        double* dst = out.data();

        for (int l = 0; l < out_order; ++l) {
            const double factor =
                std::sqrt(4.0 * pi / static_cast<double>(2 * l + 1)) *
                data_[static_cast<std::size_t>(l)];
            const int base = l * l + l;
            for (int m = -l; m <= l; ++m) {
                dst[base + m] = factor * src[base + m];
            }
        }
        return out;
    }

private:
    int order_;
    std::vector<double> data_;

    void validate(int l) const {
        if (l < 0 || l >= order_) {
            throw std::out_of_range("Invalid ZH band index l");
        }
    }
};

} // namespace sylph
