#include "render/material.hpp"

#include <algorithm>

namespace render
{

Material::Material()
    : base_color_(1.0f, 1.0f, 1.0f),
      metallic_(0.0f),
      roughness_(0.5f)
{
}

// ------------------------------------------------------------
// Base color
// ------------------------------------------------------------

void Material::set_base_color(
    const Eigen::Vector3f& color)
{
    base_color_ = color;
}

const Eigen::Vector3f& Material::base_color() const
{
    return base_color_;
}

// ------------------------------------------------------------
// Metallic
// ------------------------------------------------------------

void Material::set_metallic(float metallic)
{
    metallic_ = std::clamp(
        metallic,
        0.0f,
        1.0f
    );
}

float Material::metallic() const
{
    return metallic_;
}

// ------------------------------------------------------------
// Roughness
// ------------------------------------------------------------

void Material::set_roughness(float roughness)
{
    roughness_ = std::clamp(
        roughness,
        0.0f,
        1.0f
    );
}

float Material::roughness() const
{
    return roughness_;
}

} // namespace render