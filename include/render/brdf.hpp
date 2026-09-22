/**
 * @file brdf.hpp
 * @brief Manages the Bidirectional Reflectance Distribution Function (BRDF) data.
 * 
 * In this implementation, the BRDF is assumed to be circularly symmetric (e.g., Phong).
 * This allows it to be represented efficiently using Zonal Harmonics rather than full 
 * Spherical Harmonics, drastically reducing memory and evaluation costs.
 */
#pragma once

#include "sylph/sh_coefficients.hpp"

#include <stdexcept>

namespace render {

class BRDF {
public:
    explicit BRDF(int order)
        : coefficients_(order)
    {
        if (order <= 0) {
            throw std::invalid_argument(
                "BRDF order must be positive"
            );
        }
    }

    [[nodiscard]]
    int order() const noexcept
    {
        return coefficients_.order();
    }

    [[nodiscard]]
    const sylph::SHCoefficients& coefficients() const noexcept
    {
        return coefficients_;
    }

    [[nodiscard]]
    sylph::SHCoefficients& coefficients() noexcept
    {
        return coefficients_;
    }

    void set_coefficients(
        const sylph::SHCoefficients& coefficients
    )
    {
        coefficients_ = coefficients;
    }

    void clear() noexcept
    {
        for (int l = 0; l < order(); ++l) {
            for (int m = -l; m <= l; ++m) {
                coefficients_(l, m) = 0.0;
            }
        }
    }

private:
    sylph::SHCoefficients coefficients_;
};

} // namespace render