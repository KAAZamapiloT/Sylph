#include "render/sh_lobe_mesh.hpp"

#include "render/sh_rendering.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace render
{

namespace
{

struct SHSample
{
    Eigen::Vector3f direction;
    float value;
};

// -----------------------------------------------------------------------------
// Sample the SH function on a UV sphere.
// -----------------------------------------------------------------------------

std::vector<SHSample> sample_sh(
    const sylph::SHCoefficients& coefficients,
    int segments,
    int rings)
{
    if (segments < 3)
    {
        throw std::invalid_argument(
            "SH lobe requires at least 3 segments"
        );
    }

    if (rings < 2)
    {
        throw std::invalid_argument(
            "SH lobe requires at least 2 rings"
        );
    }

    std::vector<SHSample> samples;

    samples.reserve(
        static_cast<std::size_t>(rings + 1) *
        static_cast<std::size_t>(segments)
    );

    constexpr double pi =
        std::numbers::pi_v<double>;

    // -------------------------------------------------------------------------
    // theta:
    //
    // 0   -> north pole
    // pi  -> south pole
    //
    // phi:
    //
    // 0 .. 2pi
    // -------------------------------------------------------------------------

    for (int y = 0; y <= rings; ++y)
    {
        const double theta =
            pi *
            static_cast<double>(y) /
            static_cast<double>(rings);

        const double sin_theta =
            std::sin(theta);

        const double cos_theta =
            std::cos(theta);

        for (int x = 0; x < segments; ++x)
        {
            const double phi =
                2.0 *
                pi *
                static_cast<double>(x) /
                static_cast<double>(segments);

            const double sin_phi =
                std::sin(phi);

            const double cos_phi =
                std::cos(phi);

            // Same spherical convention used by Environment:
            //
            // x = sin(theta) cos(phi)
            // y = cos(theta)
            // z = sin(theta) sin(phi)

            Eigen::Vector3f direction(
                static_cast<float>(
                    sin_theta * cos_phi
                ),
                static_cast<float>(
                    cos_theta
                ),
                static_cast<float>(
                    sin_theta * sin_phi
                )
            );

            // SHRendering evaluates the spherical function represented
            // by the coefficient vector.
            const double value =
                SHRendering::evaluate(
                    coefficients,
                    direction
                );

            samples.push_back(
                SHSample{
                    direction,
                    static_cast<float>(value)
                }
            );
        }
    }

    return samples;
}

// -----------------------------------------------------------------------------
// Build triangle indices for the UV sphere.
// -----------------------------------------------------------------------------

std::vector<std::uint32_t> build_indices(
    int segments,
    int rings)
{
    std::vector<std::uint32_t> indices;

    indices.reserve(
        static_cast<std::size_t>(rings) *
        static_cast<std::size_t>(segments) *
        6
    );

    for (int y = 0; y < rings; ++y)
    {
        for (int x = 0; x < segments; ++x)
        {
            const int x_next =
                (x + 1) % segments;

            const int current =
                y * segments + x;

            const int next =
                (y + 1) * segments + x;

            const int current_next =
                y * segments + x_next;

            const int next_next =
                (y + 1) * segments + x_next;

            // First triangle
            indices.push_back(
                static_cast<std::uint32_t>(current)
            );

            indices.push_back(
                static_cast<std::uint32_t>(next)
            );

            indices.push_back(
                static_cast<std::uint32_t>(current_next)
            );

            // Second triangle
            indices.push_back(
                static_cast<std::uint32_t>(current_next)
            );

            indices.push_back(
                static_cast<std::uint32_t>(next)
            );

            indices.push_back(
                static_cast<std::uint32_t>(next_next)
            );
        }
    }

    return indices;
}

} // namespace

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------

SHLobeMesh::SHLobeMesh(
    const sylph::SHCoefficients& coefficients)
{
    rebuild(coefficients);
}

// -----------------------------------------------------------------------------
// Rebuild the geometry from a new SH coefficient set.
// -----------------------------------------------------------------------------

void SHLobeMesh::rebuild(
    const sylph::SHCoefficients& coefficients)
{
    mesh_ = build_mesh(coefficients);
}

// -----------------------------------------------------------------------------
// Return the generated OpenGL mesh.
// -----------------------------------------------------------------------------

const Mesh&
SHLobeMesh::mesh() const noexcept
{
    return *mesh_;
}

// -----------------------------------------------------------------------------
// Build the actual deformed sphere.
// -----------------------------------------------------------------------------

std::unique_ptr<Mesh>
SHLobeMesh::build_mesh(
    const sylph::SHCoefficients& coefficients)
{
    constexpr int segments = 40;
    constexpr int rings = 64;

    // -------------------------------------------------------------------------
    // Evaluate SH at every spherical direction.
    // -------------------------------------------------------------------------

    const std::vector<SHSample> samples =
        sample_sh(
            coefficients,
            segments,
            rings
        );

    // -------------------------------------------------------------------------
    // Find maximum absolute SH amplitude.
    //
    // This lets us normalize arbitrary coefficient magnitudes into
    // roughly [-1, 1] before deforming the sphere.
    // -------------------------------------------------------------------------

    float max_abs = 0.0f;

    for (const SHSample& sample : samples)
    {
        max_abs =
            std::max(
                max_abs,
                std::abs(sample.value)
            );
    }

    // Avoid division by zero for the all-zero SH function.
    if (max_abs < 1e-8f)
    {
        max_abs = 1.0f;
    }

    // -------------------------------------------------------------------------
    // Sphere deformation.
    //
    // normalized = -1 ... +1
    //
    // radius:
    //
    // value = -1  -> radius = 0.25
    // value =  0  -> radius = 1.00
    // value = +1  -> radius = 1.75
    // -------------------------------------------------------------------------

    constexpr float base_radius = 1.0f;
    constexpr float deformation = 0.75f;

    std::vector<Vertex> vertices;

    vertices.reserve(samples.size());

    for (const SHSample& sample : samples)
    {
        const float normalized =
            sample.value / max_abs;

        const float radius =
            base_radius +
            deformation * normalized;

        Vertex vertex{};

        vertex.position =
            sample.direction * radius;

        // For the first visualization, use the original sphere
        // direction as the normal.
        //
        // This is not the exact geometric normal of the deformed
        // surface, but is perfectly adequate for the first lobe viewer.
        vertex.normal =
            sample.direction;

        // UVs are not needed for this demo.
        vertex.uv =
            Eigen::Vector2f::Zero();

        vertices.push_back(vertex);
    }

    // -------------------------------------------------------------------------
    // Build topology.
    // -------------------------------------------------------------------------

    std::vector<std::uint32_t> indices =
        build_indices(
            segments,
            rings
        );

    // -------------------------------------------------------------------------
    // Create the OpenGL mesh.
    //
    // Mesh accepts spans, and vectors provide contiguous storage.
    // -------------------------------------------------------------------------

    return std::make_unique<Mesh>(
        vertices,
        indices
    );
}

} // namespace render
