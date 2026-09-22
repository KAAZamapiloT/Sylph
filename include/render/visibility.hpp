#pragma once

#include "sylph/sh_coefficients.hpp"

#include <stdexcept>

namespace render {

class Visibility {
public:
    explicit Visibility(int order)
        : coefficients_(order)
    {
        if (order <= 0) {
            throw std::invalid_argument(
                "Visibility order must be positive"
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