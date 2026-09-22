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