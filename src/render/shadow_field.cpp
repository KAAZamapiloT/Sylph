#include "render/shadow_field.hpp"

#include "sylph/sh_grid_transform.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace render {

Eigen::Vector3f ShadowField::evaluate(
    const Lighting& lighting,
    const Visibility& self_visibility,
    const Eigen::Vector3f& position,
    std::span<const OOF* const> oofs
) const
{
    if (oofs.empty()) {
        throw std::invalid_argument(
            "ShadowField requires at least one OOF"
        );
    }

    sylph::SHGridTransform transform;

    Eigen::Vector3f result =
        Eigen::Vector3f::Zero();

    for (int channel = 0; channel < 3; ++channel) {

        const auto& lighting_coefficients =
            lighting.coefficients(
                static_cast<LightChannel>(channel)
            );

        std::vector<
            const sylph::SHCoefficients*
        > inputs;

        inputs.reserve(
            2 + oofs.size()
        );

        // Lighting
        inputs.push_back(
            &lighting_coefficients
        );

        // Self visibility
        inputs.push_back(
            &self_visibility.coefficients()
        );

        // Visibility from dynamic occluders
        for (const OOF* oof : oofs) {

            if (oof == nullptr) {
                throw std::invalid_argument(
                    "ShadowField received a null OOF"
                );
            }

            const Visibility& visibility =
                oof->lookup(position);

            inputs.push_back(
                &visibility.coefficients()
            );
        }

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


Eigen::Vector3f ShadowField::evaluate(
    const Lighting& lighting,
    const Visibility& self_visibility,
    const Eigen::Vector3f& position,
    std::span<const OOF* const> oofs,
    int t
) const
{
    if (oofs.empty()) {
        throw std::invalid_argument(
            "ShadowField requires at least one OOF"
        );
    }

    if (t <= 0) {
        throw std::invalid_argument(
            "ShadowField resolution parameter t "
            "must be positive"
        );
    }

    sylph::SHGridTransform transform;

    Eigen::Vector3f result =
        Eigen::Vector3f::Zero();

    for (int channel = 0; channel < 3; ++channel) {

        const auto& lighting_coefficients =
            lighting.coefficients(
                static_cast<LightChannel>(channel)
            );

        std::vector<
            const sylph::SHCoefficients*
        > inputs;

        inputs.reserve(
            2 + oofs.size()
        );

        // Lighting
        inputs.push_back(
            &lighting_coefficients
        );

        // Self visibility
        inputs.push_back(
            &self_visibility.coefficients()
        );

        // Visibility from dynamic occluders
        for (const OOF* oof : oofs) {

            if (oof == nullptr) {
                throw std::invalid_argument(
                    "ShadowField received a null OOF"
                );
            }

            const Visibility& visibility =
                oof->lookup(position);

            inputs.push_back(
                &visibility.coefficients()
            );
        }

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


