#include "render/obj_loader.hpp"

#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>

#include <Eigen/Core>

#include <cstdint>
#include <filesystem>
#include <format>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace render
{

namespace
{

struct VertexKey
{
    int vertex_index = -1;
    int normal_index = -1;
    int texcoord_index = -1;

    bool operator==(const VertexKey& other) const noexcept
    {
        return vertex_index == other.vertex_index &&
               normal_index == other.normal_index &&
               texcoord_index == other.texcoord_index;
    }
};

struct VertexKeyHash
{
    std::size_t operator()(const VertexKey& key) const noexcept
    {
        std::size_t h =
            std::hash<int>{}(key.vertex_index);

        h ^=
            std::hash<int>{}(key.normal_index) +
            0x9e3779b9u +
            (h << 6u) +
            (h >> 2u);

        h ^=
            std::hash<int>{}(key.texcoord_index) +
            0x9e3779b9u +
            (h << 6u) +
            (h >> 2u);

        return h;
    }
};

} // namespace

void ObjLoader::validate_path(
    const std::filesystem::path& path)
{
    if (path.empty())
    {
        throw std::invalid_argument(
            "ObjLoader: empty path."
        );
    }

    if (!std::filesystem::exists(path))
    {
        throw std::runtime_error(
            std::format(
                "ObjLoader: file does not exist: {}",
                path.string()
            )
        );
    }

    if (!std::filesystem::is_regular_file(path))
    {
        throw std::runtime_error(
            std::format(
                "ObjLoader: path is not a regular file: {}",
                path.string()
            )
        );
    }
}

MeshData ObjLoader::load(
    const std::filesystem::path& path)
{
    validate_path(path);

    // ---------------------------------------------------------
    // tinyobjloader legacy / v1.x API
    // ---------------------------------------------------------

    tinyobj::attrib_t attrib;

    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string err;

    const std::string basedir =
        path.parent_path().string();

    const bool success =
        tinyobj::LoadObj(
            &attrib,
            &shapes,
            &materials,
            &err,
            path.string().c_str(),
            basedir.c_str(),
            true
        );

    if (!success)
    {
        throw std::runtime_error(
            std::format(
                "ObjLoader: failed to parse '{}': {}",
                path.string(),
                err
            )
        );
    }

    MeshData result;

    std::unordered_map<
        VertexKey,
        std::uint32_t,
        VertexKeyHash
    > vertex_cache;

    vertex_cache.reserve(4096);

    auto create_vertex =
        [&](const tinyobj::index_t& index)
        -> std::uint32_t
    {
        const VertexKey key{
            index.vertex_index,
            index.normal_index,
            index.texcoord_index
        };

        const auto found =
            vertex_cache.find(key);

        if (found != vertex_cache.end())
        {
            return found->second;
        }

        // IMPORTANT:
        // MeshData uses render::Vertex directly.
        // There is no MeshData::Vertex nested type.
        Vertex vertex{};

        // -----------------------------------------------------
        // Position
        // -----------------------------------------------------

        if (index.vertex_index < 0)
        {
            throw std::runtime_error(
                "ObjLoader: invalid negative vertex index."
            );
        }

        const std::size_t position_base =
            static_cast<std::size_t>(
                index.vertex_index
            ) * 3u;

        if (position_base + 2u >=
            attrib.vertices.size())
        {
            throw std::runtime_error(
                "ObjLoader: vertex index out of range."
            );
        }

        vertex.position.x() =
            attrib.vertices[position_base + 0u];

        vertex.position.y() =
            attrib.vertices[position_base + 1u];

        vertex.position.z() =
            attrib.vertices[position_base + 2u];

        // -----------------------------------------------------
        // Normal
        // -----------------------------------------------------

        if (index.normal_index >= 0)
        {
            const std::size_t normal_base =
                static_cast<std::size_t>(
                    index.normal_index
                ) * 3u;

            if (normal_base + 2u >=
                attrib.normals.size())
            {
                throw std::runtime_error(
                    "ObjLoader: normal index out of range."
                );
            }

            vertex.normal.x() =
                attrib.normals[normal_base + 0u];

            vertex.normal.y() =
                attrib.normals[normal_base + 1u];

            vertex.normal.z() =
                attrib.normals[normal_base + 2u];
        }

        // -----------------------------------------------------
        // UV
        // -----------------------------------------------------

        if (index.texcoord_index >= 0)
        {
            const std::size_t uv_base =
                static_cast<std::size_t>(
                    index.texcoord_index
                ) * 2u;

            if (uv_base + 1u >=
                attrib.texcoords.size())
            {
                throw std::runtime_error(
                    "ObjLoader: texcoord index out of range."
                );
            }

            vertex.uv.x() =
                attrib.texcoords[uv_base + 0u];

            vertex.uv.y() =
                attrib.texcoords[uv_base + 1u];
        }

        const std::uint32_t new_index =
            static_cast<std::uint32_t>(
                result.vertices.size()
            );

        result.vertices.push_back(vertex);
        vertex_cache.emplace(key, new_index);

        return new_index;
    };

    // ---------------------------------------------------------
    // OBJ shapes -> indexed triangles
    //
    // LoadObj(..., true) triangulates faces, so every sequence
    // below is emitted as triangles.
    // ---------------------------------------------------------

    for (const tinyobj::shape_t& shape : shapes)
    {
        for (const tinyobj::index_t& index :
             shape.mesh.indices)
        {
            result.indices.push_back(
                create_vertex(index)
            );
        }
    }

    // ---------------------------------------------------------
    // Generate normals if the OBJ did not provide them.
    // ---------------------------------------------------------

    bool missing_normals = false;

    for (const Vertex& vertex :
         result.vertices)
    {
        if (vertex.normal.squaredNorm() < 1e-12f)
        {
            missing_normals = true;
            break;
        }
    }

    if (missing_normals)
    {
        for (Vertex& vertex :
             result.vertices)
        {
            vertex.normal.setZero();
        }

        for (std::size_t i = 0;
             i + 2u < result.indices.size();
             i += 3u)
        {
            const std::uint32_t ia =
                result.indices[i + 0u];

            const std::uint32_t ib =
                result.indices[i + 1u];

            const std::uint32_t ic =
                result.indices[i + 2u];

            const Eigen::Vector3f& a =
                result.vertices[ia].position;

            const Eigen::Vector3f& b =
                result.vertices[ib].position;

            const Eigen::Vector3f& c =
                result.vertices[ic].position;

            // Area-weighted face normal.
            const Eigen::Vector3f face_normal =
                (b - a).cross(c - a);

            result.vertices[ia].normal += face_normal;
            result.vertices[ib].normal += face_normal;
            result.vertices[ic].normal += face_normal;
        }

        for (Vertex& vertex :
             result.vertices)
        {
            const float length_sq =
                vertex.normal.squaredNorm();

            if (length_sq > 1e-12f)
            {
                vertex.normal.normalize();
            }
            else
            {
                // Fallback for isolated / degenerate vertices.
                vertex.normal =
                    Eigen::Vector3f::UnitY();
            }
        }
    }

    validate_mesh(
        result,
        path
    );

    return result;
}

void ObjLoader::validate_mesh(
    const MeshData& mesh,
    const std::filesystem::path& path)
{
    if (mesh.vertices.empty())
    {
        throw std::runtime_error(
            std::format(
                "ObjLoader: '{}' contains no vertices.",
                path.string()
            )
        );
    }

    if (mesh.indices.empty())
    {
        throw std::runtime_error(
            std::format(
                "ObjLoader: '{}' contains no indices.",
                path.string()
            )
        );
    }

    if (mesh.indices.size() % 3u != 0u)
    {
        throw std::runtime_error(
            std::format(
                "ObjLoader: '{}' index count is not divisible by 3.",
                path.string()
            )
        );
    }

    for (const std::uint32_t index :
         mesh.indices)
    {
        if (index >= mesh.vertices.size())
        {
            throw std::runtime_error(
                std::format(
                    "ObjLoader: '{}' contains out-of-range "
                    "index {}.",
                    path.string(),
                    index
                )
            );
        }
    }
}

} // namespace render
