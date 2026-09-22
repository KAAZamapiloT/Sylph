/**
 * @file glossy_relighting.hpp
 * @brief Implementation of Glossy PRT (Section 7.3.1 of the paper).
 * 
 * Calculates the Triple Product Integral: Incident Light * Material BRDF * Visibility.
 * By combining these three Spherical Harmonic functions on the fly, it achieves 
 * real-time glossy reflections with self-shadowing.
 */
#pragma once

#include "render/brdf.hpp"
#include "render/lighting.hpp"
#include "render/visibility.hpp"

#include <Eigen/Core>

namespace render {

class GlossyRelighting {
public:
    [[nodiscard]]
    Eigen::Vector3f evaluate(
        const Lighting& lighting,
        const BRDF& brdf,
        const Visibility& visibility
    ) const;

    [[nodiscard]]
    Eigen::Vector3f evaluate(
        const Lighting& lighting,
        const BRDF& brdf,
        const Visibility& visibility,
        int t
    ) const;
};

} // namespace render