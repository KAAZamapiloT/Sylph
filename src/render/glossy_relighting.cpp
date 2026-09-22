#include "render/glossy_relighting.hpp"

#include "sylph/sh_grid_transform.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <span>

namespace render {

Eigen::Vector3f GlossyRelighting::evaluate(
    const Lighting& lighting,
    const BRDF& brdf,
    const Visibility& visibility
) const
{
    sylph::SHGridTransform transform;

    Eigen::Vector3f result =
        Eigen::Vector3f::Zero();

    for (int channel = 0; channel < 3; ++channel)
    {
        const auto& lighting_coefficients =
            lighting.coefficients(
                static_cast<LightChannel>(channel)
            );

        const std::array<
            const sylph::SHCoefficients*,
            3
        > inputs = {
            &lighting_coefficients,
            &brdf.coefficients(),
            &visibility.coefficients()
        };

        result[channel] =
            static_cast<float>(
                transform.product_integral(
                    std::span<
                        const sylph::SHCoefficients* const
                    >(inputs)
                )
            );
    }

    return result;
}


Eigen::Vector3f GlossyRelighting::evaluate(
    const Lighting& lighting,
    const BRDF& brdf,
    const Visibility& visibility,
    int t
) const
{
    if (t <= 0)
    {
        throw std::invalid_argument(
            "GlossyRelighting parameter t "
            "must be positive"
        );
    }

    sylph::SHGridTransform transform;

    Eigen::Vector3f result =
        Eigen::Vector3f::Zero();

    for (int channel = 0; channel < 3; ++channel)
    {
        const auto& lighting_coefficients =
            lighting.coefficients(
                static_cast<LightChannel>(channel)
            );

        const std::array<
            const sylph::SHCoefficients*,
            3
        > inputs = {
            &lighting_coefficients,
            &brdf.coefficients(),
            &visibility.coefficients()
        };

        const auto product =
            transform.product(
                std::span<
                    const sylph::SHCoefficients* const
                >(inputs),
                1,
                t
            );

        result[channel] =
            static_cast<float>(
                2.0 *
                std::sqrt(sylph::pi) *
                product(0, 0)
            );
    }

    return result;
}

} // namespace render