#include "render/enviroment.hpp"

#include "sylph/spherical_harmonics.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace render
{

namespace
{
constexpr float TWO_PI =
    2.0f * std::numbers::pi_v<float>;

} // namespace

// ------------------------------------------------------------
// Construction
// ------------------------------------------------------------

Environment::Environment(
    std::size_t width,
    std::size_t height)
{
    resize(width, height);
}

// ------------------------------------------------------------
// Dimensions
// ------------------------------------------------------------

std::size_t Environment::width() const
{
    return width_;
}

std::size_t Environment::height() const
{
    return height_;
}

// ------------------------------------------------------------
// Resize
// ------------------------------------------------------------

void Environment::resize(
    std::size_t width,
    std::size_t height)
{
    if (width == 0 || height == 0) {
        throw std::invalid_argument(
            "Environment dimensions must be greater than zero."
        );
    }

    width_ = width;
    height_ = height;

    data_.assign(
        width_ * height_,
        Eigen::Vector3f::Zero()
    );
}

// ------------------------------------------------------------
// Clear
// ------------------------------------------------------------

void Environment::clear()
{
    std::fill(
        data_.begin(),
        data_.end(),
        Eigen::Vector3f::Zero()
    );
}

// ------------------------------------------------------------
// Texel access
// ------------------------------------------------------------

void Environment::set_texel(
    std::size_t x,
    std::size_t y,
    const Eigen::Vector3f& radiance)
{
    if (x >= width_ || y >= height_) {
        throw std::out_of_range(
            "Environment texel coordinates are out of range."
        );
    }

    data_[index(x, y)] = radiance;
}

const Eigen::Vector3f& Environment::texel(
    std::size_t x,
    std::size_t y) const
{
    if (x >= width_ || y >= height_) {
        throw std::out_of_range(
            "Environment texel coordinates are out of range."
        );
    }

    return data_[index(x, y)];
}

// ------------------------------------------------------------
// Directional sampling
// ------------------------------------------------------------

Eigen::Vector3f Environment::sample(
    const Eigen::Vector3f& direction) const
{
    if (width_ == 0 || height_ == 0) {
        throw std::logic_error(
            "Cannot sample an empty environment."
        );
    }

    if (direction.squaredNorm() == 0.0f) {
        throw std::invalid_argument(
            "Environment sampling direction cannot be zero."
        );
    }

    const Eigen::Vector3f normalized =
        direction.normalized();

    const Eigen::Vector3f spherical =
        direction_to_spherical(normalized);

    const float theta = spherical.x();
    float phi = spherical.y();

    // [0, 2pi)
    if (phi < 0.0f) {
        phi += TWO_PI;
    }

    // Map spherical coordinates to image coordinates.
    //
    // u = 0..1 around longitude
    // v = 0..1 from north to south
    //
    // We use texel-center coordinates for interpolation.

    const float u =
        phi / TWO_PI;

    const float v =
        theta / std::numbers::pi_v<float>;

    const float fx =
        u * static_cast<float>(width_) - 0.5f;

    const float fy =
        v * static_cast<float>(height_) - 0.5f;

    const auto floor_x =
        static_cast<long long>(std::floor(fx));

    const auto floor_y =
        static_cast<long long>(std::floor(fy));

    const long long x0 = floor_x;
    const long long x1 = floor_x + 1;

    const long long y0 = floor_y;
    const long long y1 = floor_y + 1;

    const float tx = fx - static_cast<float>(x0);
    const float ty = fy - static_cast<float>(y0);

    // Longitude wraps around.
    const auto wrap_x =
        [this](long long x) -> std::size_t
        {
            const long long w =
                static_cast<long long>(width_);

            x %= w;

            if (x < 0) {
                x += w;
            }

            return static_cast<std::size_t>(x);
        };

    // Latitude is clamped at the poles.
    const auto clamp_y =
        [this](long long y) -> std::size_t
        {
            const long long h =
                static_cast<long long>(height_);

            y = std::clamp(
                y,
                0LL,
                h - 1
            );

            return static_cast<std::size_t>(y);
        };

    const Eigen::Vector3f& c00 =
        data_[index(
            wrap_x(x0),
            clamp_y(y0)
        )];

    const Eigen::Vector3f& c10 =
        data_[index(
            wrap_x(x1),
            clamp_y(y0)
        )];

    const Eigen::Vector3f& c01 =
        data_[index(
            wrap_x(x0),
            clamp_y(y1)
        )];

    const Eigen::Vector3f& c11 =
        data_[index(
            wrap_x(x1),
            clamp_y(y1)
        )];

    const Eigen::Vector3f top =
        c00 + tx * (c10 - c00);

    const Eigen::Vector3f bottom =
        c01 + tx * (c11 - c01);

    return top + ty * (bottom - top);
}

// ------------------------------------------------------------
// Project environment to SH
// ------------------------------------------------------------

Lighting Environment::project_to_sh(int order) const
{
    if (width_ == 0 || height_ == 0) {
        throw std::logic_error(
            "Cannot project an empty environment."
        );
    }

    if (order <= 0) {
        throw std::invalid_argument(
            "SH order must be positive."
        );
    }

    Lighting lighting(order);

    sylph::SphericalHarmonics sh;

    const double dtheta =
        std::numbers::pi_v<double> /
        static_cast<double>(height_);

    const double dphi =
        2.0 * std::numbers::pi_v<double> /
        static_cast<double>(width_);

    /*
        Midpoint quadrature on the equirectangular domain.

        f_l^m =
            integral F(theta,phi) Y_l^m(theta,phi)
            sin(theta) dtheta dphi

        Each texel is sampled at its center:

            theta = (y + 0.5) * dtheta
            phi   = (x + 0.5) * dphi
    */

    for (std::size_t y = 0; y < height_; ++y) {

        const double theta =
            (static_cast<double>(y) + 0.5) * dtheta;

        const double sin_theta =
            std::sin(theta);

        for (std::size_t x = 0; x < width_; ++x) {

            const double phi =
                (static_cast<double>(x) + 0.5) * dphi;

            const Eigen::Vector3f& radiance =
                data_[index(x, y)];

            const double weight =
                sin_theta * dtheta * dphi;

            for (int l = 0; l < order; ++l) {

                for (int m = -l; m <= l; ++m) {

                    const double basis =
                        sh.evaluate(
                            l,
                            m,
                            theta,
                            phi
                        );

                    const double contribution =
                        weight * basis;

                    lighting.coefficients(
                        LightChannel::Red
                    )(l, m) +=
                        static_cast<double>(
                            radiance.x()
                        ) *
                        contribution;

                    lighting.coefficients(
                        LightChannel::Green
                    )(l, m) +=
                        static_cast<double>(
                            radiance.y()
                        ) *
                        contribution;

                    lighting.coefficients(
                        LightChannel::Blue
                    )(l, m) +=
                        static_cast<double>(
                            radiance.z()
                        ) *
                        contribution;
                }
            }
        }
    }

    return lighting;
}

// ------------------------------------------------------------
// Linear index
// ------------------------------------------------------------

std::size_t Environment::index(
    std::size_t x,
    std::size_t y) const
{
    return y * width_ + x;
}

// ------------------------------------------------------------
// Direction → spherical coordinates
//
// Returns:
//
//   x = theta
//   y = phi
//   z = unused
//
// Convention:
//
//   x = sin(theta) cos(phi)
//   y = cos(theta)
//   z = sin(theta) sin(phi)
// ------------------------------------------------------------

Eigen::Vector3f
Environment::direction_to_spherical(
    const Eigen::Vector3f& direction)
{
    const float theta =
        std::acos(
            std::clamp(
                direction.y(),
                -1.0f,
                1.0f
            )
        );

    const float phi =
        std::atan2(
            direction.z(),
            direction.x()
        );

    return Eigen::Vector3f(
        theta,
        phi,
        0.0f
    );
}

} // namespace render