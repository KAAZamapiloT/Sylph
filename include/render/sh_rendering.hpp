/**
 * @file sh_rendering.hpp
 * @brief High-level Rendering logic for SH-based objects.
 * 
 * Handles the specific shader setups and rendering passes required to draw PRT meshes 
 * and SH Lobes, abstracting the complex uniform bindings away from the main loop.
 */
#pragma once

#include "render/lighting.hpp"
#include "sylph/sh_coefficients.hpp"

#include <Eigen/Dense>

namespace render
{

class SHRendering
{
public:
    // --------------------------------------------------------
    // Evaluate a scalar SH function at a 3D direction.
    //
    // direction does not need to be normalized.
    // A zero-length direction is invalid.
    // --------------------------------------------------------

    static double evaluate(
        const sylph::SHCoefficients& coefficients,
        const Eigen::Vector3f& direction
    );

    // --------------------------------------------------------
    // Evaluate RGB SH lighting at a 3D direction.
    // --------------------------------------------------------

    static Eigen::Vector3f evaluate(
        const Lighting& lighting,
        const Eigen::Vector3f& direction
    );

private:
    static Eigen::Vector2d direction_to_spherical(
        const Eigen::Vector3f& direction
    );
};

} // namespace render