#pragma once
#include "render/mesh_data.hpp"
#include "render/prt_mesh.hpp"
#include <memory>

namespace render {

class PRTBaker {
public:
    // Bakes diffuse unshadowed OR shadowed transfer function into a PRTMesh
    // order: SH order (e.g. 4 for 16 coefficients)
    // num_rays: Rays per vertex for Monte Carlo integration
    static std::shared_ptr<PRTMesh> bake(
        const MeshData& mesh, 
        int order = 4, 
        int num_rays = 1000,
        bool shadows = true
    );
};

}
