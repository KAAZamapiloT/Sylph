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