#include "render/prt_mesh.hpp"
#include <utility>
#include <cstddef>

namespace render
{

PRTMesh::PRTMesh(
    std::span<const PRTVertex> vertices,
    std::span<const std::uint32_t> indices)
    : vertex_count_(vertices.size()),
      index_count_(indices.size())
{
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size_bytes(), vertices.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &ebo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size_bytes(), indices.data(), GL_STATIC_DRAW);

    // Location 0: position
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(PRTVertex), reinterpret_cast<const void*>(offsetof(PRTVertex, position)));

    // Location 1: normal
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(PRTVertex), reinterpret_cast<const void*>(offsetof(PRTVertex, normal)));

    // Locations 2,3,4,5: SH coefficients (4 vec4s = 16 floats)
    std::size_t sh_offset = offsetof(PRTVertex, sh);
    for (int i = 0; i < 4; ++i) {
        glEnableVertexAttribArray(2 + i);
        glVertexAttribPointer(2 + i, 4, GL_FLOAT, GL_FALSE, sizeof(PRTVertex), reinterpret_cast<const void*>(sh_offset + i * 4 * sizeof(float)));
    }

    glBindVertexArray(0);
}

PRTMesh::~PRTMesh()
{
    if (ebo_) glDeleteBuffers(1, &ebo_);
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
}

PRTMesh::PRTMesh(PRTMesh&& other) noexcept
    : vao_(other.vao_), vbo_(other.vbo_), ebo_(other.ebo_),
      vertex_count_(other.vertex_count_), index_count_(other.index_count_)
{
    other.vao_ = 0; other.vbo_ = 0; other.ebo_ = 0;
}

PRTMesh& PRTMesh::operator=(PRTMesh&& other) noexcept
{
    if (this != &other) {
        this->~PRTMesh();
        vao_ = other.vao_; vbo_ = other.vbo_; ebo_ = other.ebo_;
        vertex_count_ = other.vertex_count_; index_count_ = other.index_count_;
        other.vao_ = 0; other.vbo_ = 0; other.ebo_ = 0;
    }
    return *this;
}

void PRTMesh::draw() const
{
    if (vao_ == 0 || index_count_ == 0) return;
    glBindVertexArray(vao_);
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(index_count_), GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
}

} // namespace render
