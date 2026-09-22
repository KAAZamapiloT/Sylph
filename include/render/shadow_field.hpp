/**
 * @file shadow_field.hpp
 * @brief Implementation of Shadow Fields (Section 7.3.2 of the paper).
 * 
 * Evaluates the N-way Multiple Product Integral of SH functions. Used to combine 
 * Incident Lighting, Self-Visibility, and multiple Object Occlusion Fields (OOFs) 
 * to cast mathematically accurate soft shadows between dynamic objects.
 */
#pragma once

#include "render/lighting.hpp"
#include "render/oof.hpp"
#include "render/visibility.hpp"

#include <Eigen/Core>

#include <span>

namespace render {

class ShadowField {
public:
    [[nodiscard]]
    Eigen::Vector3f evaluate(
        const Lighting& lighting,
        const Visibility& self_visibility,
        const Eigen::Vector3f& position,
        std::span<const OOF* const> oofs
    ) const;

    [[nodiscard]]
    Eigen::Vector3f evaluate(
        const Lighting& lighting,
        const Visibility& self_visibility,
        const Eigen::Vector3f& position,
        std::span<const OOF* const> oofs,
        int t
    ) const;
};

} // namespace render