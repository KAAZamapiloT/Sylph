/**
 * @file sh_lobe_mesh.hpp
 * @brief Procedural Spherical Harmonic Lobe Generator.
 * 
 * Dynamically generates a 3D mesh that visualizes an SH function in real space. 
 * The mesh's radius at any given angle is determined by the evaluated SH magnitude, 
 * with colors indicating positive (green) or negative (red) values.
 */
#pragma once

#include "render/mesh.hpp"
#include "sylph/sh_coefficients.hpp"

#include <memory>

namespace render
{

class SHLobeMesh
{
public:
    explicit SHLobeMesh(
        const sylph::SHCoefficients& coefficients
    );

    void rebuild(
        const sylph::SHCoefficients& coefficients
    );

    [[nodiscard]]
    const Mesh& mesh() const noexcept;

private:
    static std::unique_ptr<Mesh>
    build_mesh(
        const sylph::SHCoefficients& coefficients
    );

private:
    std::unique_ptr<Mesh> mesh_;
};

} // namespace render