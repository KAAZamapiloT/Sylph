#pragma once

#include "render/mesh.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace render
{

// CPU-side representation of geometry.
// It intentionally reuses render::Vertex so it can be passed
// directly to the existing Mesh constructor.
struct MeshData
{
    std::vector<Vertex> vertices;
    std::vector<std::uint32_t> indices;

    void clear() noexcept
    {
        vertices.clear();
        indices.clear();
    }

    [[nodiscard]]
    bool empty() const noexcept
    {
        return vertices.empty() || indices.empty();
    }

    [[nodiscard]]
    std::size_t vertex_count() const noexcept
    {
        return vertices.size();
    }

    [[nodiscard]]
    std::size_t index_count() const noexcept
    {
        return indices.size();
    }

    [[nodiscard]]
    std::size_t triangle_count() const noexcept
    {
        return indices.size() / 3;
    }
};

} // namespace render
