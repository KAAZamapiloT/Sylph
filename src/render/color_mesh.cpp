#include "render/color_mesh.hpp"
#include <utility>
#include <stdexcept>

namespace render
{

ColorMesh::ColorMesh(
    std::span<const ColorVertex> vertices,
    std::span<const std::uint32_t> indices)
    : vertex_count_(vertices.size()),
      index_count_(indices.size())
{
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size_bytes(), vertices.data(), GL_DYNAMIC_DRAW);

    glGenBuffers(1, &ebo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size_bytes(), indices.data(), GL_STATIC_DRAW);

    // Position - location 0
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ColorVertex), reinterpret_cast<const void*>(offsetof(ColorVertex, position)));

    // Color - location 1
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ColorVertex), reinterpret_cast<const void*>(offsetof(ColorVertex, color)));

    glBindVertexArray(0);
}

ColorMesh::~ColorMesh()
{
    if (ebo_) glDeleteBuffers(1, &ebo_);
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
}

ColorMesh::ColorMesh(ColorMesh&& other) noexcept
    : vao_(other.vao_), vbo_(other.vbo_), ebo_(other.ebo_),
      vertex_count_(other.vertex_count_), index_count_(other.index_count_)
{
    other.vao_ = 0; other.vbo_ = 0; other.ebo_ = 0;
}

ColorMesh& ColorMesh::operator=(ColorMesh&& other) noexcept
{
    if (this != &other) {
        this->~ColorMesh();
        vao_ = other.vao_; vbo_ = other.vbo_; ebo_ = other.ebo_;
        vertex_count_ = other.vertex_count_; index_count_ = other.index_count_;
        other.vao_ = 0; other.vbo_ = 0; other.ebo_ = 0;
    }
    return *this;
}

void ColorMesh::update_colors(std::span<const ColorVertex> vertices)
{
    if (vertices.size() != vertex_count_) return;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, vertices.size_bytes(), vertices.data());
}

void ColorMesh::draw() const
{
    if (vao_ == 0 || index_count_ == 0) return;
    glBindVertexArray(vao_);
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(index_count_), GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
}

}
