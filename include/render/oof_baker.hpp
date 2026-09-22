#pragma once

#include "render/oof.hpp"
#include "render/mesh_data.hpp"
#include <memory>

namespace render {

class OOFBaker {
public:
    /**
     * @brief Bakes an Object Occlusion Field (OOF) using CPU raytracing.
     * 
     * @param mesh The target mesh to act as an occluder.
     * @param resolution The resolution of the 3D grid (e.g. 10 for a 10x10x10 grid).
     * @param visibility_order The SH order for the baked visibility (e.g. 4 or 8).
     * @param num_rays Number of spherical rays to cast per voxel (higher = better quality).
     * @return A unique_ptr to the baked OOF.
     */
    static std::unique_ptr<OOF> bake(
        const MeshData& mesh,
        int resolution,
        int visibility_order,
        int num_rays = 500
    );
};

} // namespace render
