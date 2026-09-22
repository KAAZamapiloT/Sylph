#include "render/sh_rendering.hpp"

#include "sylph/spherical_harmonics.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace render
{

namespace
{

constexpr double PI =
    std::numbers::pi_v<double>;

} // namespace

// ------------------------------------------------------------
// Scalar SH evaluation
// ------------------------------------------------------------

double SHRendering::evaluate(
    const sylph::SHCoefficients& coefficients,
    const Eigen::Vector3f& direction)
{
    if (direction.squaredNorm() == 0.0f) {
        throw std::invalid_argument(
            "SH evaluation direction cannot be zero."
        );
    }

    if (coefficients.order() <= 0) {
        throw std::invalid_argument(
            "SH coefficient order must be positive."
        );
    }

    const Eigen::Vector2d spherical =
        direction_to_spherical(direction);

    const double theta = spherical.x();
    const double phi = spherical.y();

    sylph::SphericalHarmonics sh;

    return sh.reconstruct(
        coefficients,
        theta,
        phi
    );
}

// ------------------------------------------------------------
// RGB SH evaluation
// ------------------------------------------------------------

Eigen::Vector3f SHRendering::evaluate(
    const Lighting& lighting,
    const Eigen::Vector3f& direction)
{
    const double red =
        evaluate(
            lighting.coefficients(
                LightChannel::Red
            ),
            direction
        );

    const double green =
        evaluate(
            lighting.coefficients(
                LightChannel::Green
            ),
            direction
        );

    const double blue =
        evaluate(
            lighting.coefficients(
                LightChannel::Blue
            ),
            direction
        );

    return Eigen::Vector3f(
        static_cast<float>(red),
        static_cast<float>(green),
        static_cast<float>(blue)
    );
}

// ------------------------------------------------------------
// Direction → spherical coordinates
//
// Convention used throughout the renderer:
//
//     x = sin(theta) cos(phi)
//     y = cos(theta)
//     z = sin(theta) sin(phi)
//
// Therefore:
//
//     theta = acos(y)
//     phi   = atan2(z, x)
//
// Returned range:
//
//     theta ∈ [0, pi]
//     phi   ∈ [-pi, pi]
// ------------------------------------------------------------

Eigen::Vector2d
SHRendering::direction_to_spherical(
    const Eigen::Vector3f& direction)
{
    const Eigen::Vector3f normalized =
        direction.normalized();

    const double y =
        std::clamp(
            static_cast<double>(normalized.y()),
            -1.0,
            1.0
        );

    const double theta =
        std::acos(y);

    double phi =
        std::atan2(
            static_cast<double>(normalized.z()),
            static_cast<double>(normalized.x())
        );

    // Convert to [0, 2pi).
    if (phi < 0.0) {
        phi += 2.0 * PI;
    }

    return Eigen::Vector2d(
        theta,
        phi
    );
}

} // namespace render