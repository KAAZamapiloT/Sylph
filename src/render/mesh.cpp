#include "render/mesh.hpp"

#include <utility>
#include <stdexcept>

namespace render
{

Mesh::Mesh(
    std::span<const Vertex> vertices,
    std::span<const std::uint32_t> indices)
    : vertex_count_(vertices.size()),
      index_count_(indices.size())
{
    if (vertices.empty()) {
        throw std::invalid_argument(
            "Mesh requires at least one vertex."
        );
    }

    if (indices.empty()) {
        throw std::invalid_argument(
            "Mesh requires at least one index."
        );
    }

    // --------------------------------------------------------
    // Vertex Array Object
    // --------------------------------------------------------

    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    // --------------------------------------------------------
    // Vertex Buffer Object
    // --------------------------------------------------------

    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);

    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(
            vertices.size_bytes()
        ),
        vertices.data(),
        GL_STATIC_DRAW
    );

    // --------------------------------------------------------
    // Element Buffer Object
    // --------------------------------------------------------

    glGenBuffers(1, &ebo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);

    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(
            indices.size_bytes()
        ),
        indices.data(),
        GL_STATIC_DRAW
    );

    // --------------------------------------------------------
    // Vertex layout
    //
    // Vertex:
    //   Vector3f position
    //   Vector3f normal
    //   Vector2f uv
    // --------------------------------------------------------

    // Position - location 0
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(Vertex),
        reinterpret_cast<const void*>(
            offsetof(Vertex, position)
        )
    );

    // Normal - location 1
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(
        1,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(Vertex),
        reinterpret_cast<const void*>(
            offsetof(Vertex, normal)
        )
    );

    // UV - location 2
    glEnableVertexAttribArray(2);

    glVertexAttribPointer(
        2,
        2,
        GL_FLOAT,
        GL_FALSE,
        sizeof(Vertex),
        reinterpret_cast<const void*>(
            offsetof(Vertex, uv)
        )
    );

    // --------------------------------------------------------
    // Unbind VAO
    // --------------------------------------------------------

    glBindVertexArray(0);

    // Do NOT unbind GL_ELEMENT_ARRAY_BUFFER here.
    //
    // The EBO binding is stored inside the VAO.
}

Mesh::~Mesh()
{
    if (ebo_ != 0) {
        glDeleteBuffers(1, &ebo_);
    }

    if (vbo_ != 0) {
        glDeleteBuffers(1, &vbo_);
    }

    if (vao_ != 0) {
        glDeleteVertexArrays(1, &vao_);
    }
}

// ------------------------------------------------------------
// Move constructor
// ------------------------------------------------------------

Mesh::Mesh(Mesh&& other) noexcept
    : vao_(other.vao_),
      vbo_(other.vbo_),
      ebo_(other.ebo_),
      vertex_count_(other.vertex_count_),
      index_count_(other.index_count_)
{
    other.vao_ = 0;
    other.vbo_ = 0;
    other.ebo_ = 0;

    other.vertex_count_ = 0;
    other.index_count_ = 0;
}

// ------------------------------------------------------------
// Move assignment
// ------------------------------------------------------------

Mesh& Mesh::operator=(Mesh&& other) noexcept
{
    if (this == &other) {
        return *this;
    }

    // Release current resources.
    if (ebo_ != 0) {
        glDeleteBuffers(1, &ebo_);
    }

    if (vbo_ != 0) {
        glDeleteBuffers(1, &vbo_);
    }

    if (vao_ != 0) {
        glDeleteVertexArrays(1, &vao_);
    }

    // Take ownership.
    vao_ = other.vao_;
    vbo_ = other.vbo_;
    ebo_ = other.ebo_;

    vertex_count_ = other.vertex_count_;
    index_count_ = other.index_count_;

    // Leave source empty.
    other.vao_ = 0;
    other.vbo_ = 0;
    other.ebo_ = 0;

    other.vertex_count_ = 0;
    other.index_count_ = 0;

    return *this;
}

// ------------------------------------------------------------
// Draw
// ------------------------------------------------------------

void Mesh::draw() const
{
    if (vao_ == 0 || index_count_ == 0) {
        return;
    }

    glBindVertexArray(vao_);

    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(index_count_),
        GL_UNSIGNED_INT,
        nullptr
    );

    glBindVertexArray(0);
}

// ------------------------------------------------------------
// Information
// ------------------------------------------------------------

std::size_t Mesh::vertex_count() const
{
    return vertex_count_;
}

std::size_t Mesh::index_count() const
{
    return index_count_;
}

} 