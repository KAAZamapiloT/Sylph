/**
 * @file obj_loader.hpp
 * @brief Wavefront .obj file parsing utility.
 * 
 * Reads 3D model files from disk, extracts vertex positions, normals, and texture 
 * coordinates, and packs them into the MeshData structure for the engine to use.
 */
#pragma once

#include "render/mesh_data.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>

namespace render {

class ObjLoader
{
public:
    ObjLoader() = delete;

    [[nodiscard]]
    static MeshData load(
        const std::filesystem::path& path
    );

private:
    static void validate_path(
        const std::filesystem::path& path
    );

    static void validate_mesh(
        const MeshData& mesh,
        const std::filesystem::path& path
    );
};

} // namespace render