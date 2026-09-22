/**
 * @file color_mesh.hpp
 * @brief Simple 3D Mesh wrapper that supports per-vertex colors.
 * 
 * Used primarily for basic rendering, prototyping, or visualizations where baked 
 * PRT data is unnecessary. The mesh handles its own OpenGL VAO/VBO lifecycle.
 */
#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <span>
#include <glad/gl.h>

namespace render
{

struct ColorVertex
{
    Eigen::Vector3f position;
    Eigen::Vector3f color;
};

class ColorMesh
{
public:
    ColorMesh(
        std::span<const ColorVertex> vertices,
        std::span<const std::uint32_t> indices
    );

    ~ColorMesh();

    ColorMesh(const ColorMesh&) = delete;
    ColorMesh& operator=(const ColorMesh&) = delete;

    ColorMesh(ColorMesh&& other) noexcept;
    ColorMesh& operator=(ColorMesh&& other) noexcept;

    void update_colors(std::span<const ColorVertex> vertices);

    void draw() const;

private:
    GLuint vao_{0};
    GLuint vbo_{0};
    GLuint ebo_{0};

    std::size_t vertex_count_{0};
    std::size_t index_count_{0};
};

} // namespace render
