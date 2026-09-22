/**
 * @file mesh.hpp
 * @brief Standard 3D Mesh wrapper.
 * 
 * Manages the OpenGL Vertex Array Objects (VAOs) and Buffers (VBOs) for standard 
 * geometry (Position, Normal, UV). Used as the base for more complex mesh types.
 */
#pragma once

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <span>

#include <glad/gl.h>

namespace render
{

struct Vertex
{
    Eigen::Vector3f position;
    Eigen::Vector3f normal;
    Eigen::Vector2f uv;
};

class Mesh
{
public:
    Mesh(
        std::span<const Vertex> vertices,
        std::span<const std::uint32_t> indices
    );

    ~Mesh();

    // OpenGL resources are owned by Mesh.
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    // Movable so Mesh can be stored in containers.
    Mesh(Mesh&& other) noexcept;
    Mesh& operator=(Mesh&& other) noexcept;

    // Draw the mesh.
    void draw() const;

    // Information about the geometry.
    std::size_t vertex_count() const;
    std::size_t index_count() const;

private:
    GLuint vao_{0};
    GLuint vbo_{0};
    GLuint ebo_{0};

    std::size_t vertex_count_{0};
    std::size_t index_count_{0};
};

} 