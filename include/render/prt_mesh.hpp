/**
 * @file prt_mesh.hpp
 * @brief Precomputed Radiance Transfer Mesh.
 * 
 * A specialized mesh that stores Spherical Harmonic coefficients per-vertex. 
 * This allows each vertex to carry its own precalculated self-shadowing or 
 * inter-reflection data for advanced lighting techniques.
 */
#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <span>
#include <glad/gl.h>

namespace render
{

// Order 4 SH requires 16 floats. We pack them into 4 vec4s.
struct PRTVertex
{
    Eigen::Vector3f position;
    Eigen::Vector3f normal;
    float sh[16];
};

class PRTMesh
{
public:
    PRTMesh(
        std::span<const PRTVertex> vertices,
        std::span<const std::uint32_t> indices
    );

    ~PRTMesh();

    PRTMesh(const PRTMesh&) = delete;
    PRTMesh& operator=(const PRTMesh&) = delete;

    PRTMesh(PRTMesh&& other) noexcept;
    PRTMesh& operator=(PRTMesh&& other) noexcept;

    void draw() const;

private:
    GLuint vao_{0};
    GLuint vbo_{0};
    GLuint ebo_{0};

    std::size_t vertex_count_{0};
    std::size_t index_count_{0};
};

} // namespace render
