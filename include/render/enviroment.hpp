#pragma once

#include "render/lighting.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <vector>

namespace render
{

class Environment
{
public:
    Environment() = default;

    Environment(
        std::size_t width,
        std::size_t height
    );

    // --------------------------------------------------------
    // Environment dimensions
    // --------------------------------------------------------

    [[nodiscard]]
    std::size_t width() const;

    [[nodiscard]]
    std::size_t height() const;

    // --------------------------------------------------------
    // Data
    // --------------------------------------------------------

    void resize(
        std::size_t width,
        std::size_t height
    );

    void clear();

    void set_texel(
        std::size_t x,
        std::size_t y,
        const Eigen::Vector3f& radiance
    );

    [[nodiscard]]
    const Eigen::Vector3f& texel(
        std::size_t x,
        std::size_t y
    ) const;

    // --------------------------------------------------------
    // Directional sampling
    //
    // Direction must be normalized.
    //
    // Spherical convention:
    //
    //   x = sin(theta) cos(phi)
    //   y = cos(theta)
    //   z = sin(theta) sin(phi)
    // --------------------------------------------------------

    [[nodiscard]]
    Eigen::Vector3f sample(
        const Eigen::Vector3f& direction
    ) const;

    // --------------------------------------------------------
    // Convert environment lighting into SH.
    // --------------------------------------------------------

    [[nodiscard]]
    Lighting project_to_sh(int order) const;

private:
    [[nodiscard]]
    std::size_t index(
        std::size_t x,
        std::size_t y
    ) const;

    static Eigen::Vector3f
    direction_to_spherical(
        const Eigen::Vector3f& direction
    );

private:
    std::size_t width_{0};
    std::size_t height_{0};

    // Row-major equirectangular image.
    //
    // y = 0     -> north pole
    // y = h - 1 -> south pole
    //
    // x wraps around the sphere.
    std::vector<Eigen::Vector3f> data_;
};

} // namespace render