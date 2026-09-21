#pragma once

#include <Eigen/Dense>

namespace render
{

class Material
{
public:
    Material();

    // Base surface color.
    void set_base_color(const Eigen::Vector3f& color);
    const Eigen::Vector3f& base_color() const;

    // Metallic factor.
    void set_metallic(float metallic);
    float metallic() const;

    // Roughness factor.
    void set_roughness(float roughness);
    float roughness() const;

private:
    Eigen::Vector3f base_color_;
    float metallic_;
    float roughness_;
};

} // namespace render