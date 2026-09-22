#include "render/oof_baker.hpp"
#include "sylph/spherical_harmonics.hpp"
#include <iostream>
#include <cmath>
#include <numbers>
#include <limits>
#include <vector>

namespace render {

// Helper: Möller-Trumbore ray-triangle intersection
static bool ray_triangle_intersect(
    const Eigen::Vector3f& orig, const Eigen::Vector3f& dir,
    const Eigen::Vector3f& v0, const Eigen::Vector3f& v1, const Eigen::Vector3f& v2)
{
    const float EPSILON = 1e-6f;
    Eigen::Vector3f edge1 = v1 - v0;
    Eigen::Vector3f edge2 = v2 - v0;
    Eigen::Vector3f h = dir.cross(edge2);
    float a = edge1.dot(h);

    if (a > -EPSILON && a < EPSILON)
        return false; // Ray is parallel to triangle

    float f = 1.0f / a;
    Eigen::Vector3f s = orig - v0;
    float u = f * s.dot(h);

    if (u < 0.0f || u > 1.0f)
        return false;

    Eigen::Vector3f q = s.cross(edge1);
    float v = f * dir.dot(q);

    if (v < 0.0f || u + v > 1.0f)
        return false;

    float t = f * edge2.dot(q);
    if (t > EPSILON) 
        return true; // Hit
    
    return false;
}

std::unique_ptr<OOF> OOFBaker::bake(
    const MeshData& mesh, 
    int resolution, 
    int visibility_order, 
    int num_rays)
{
    std::cout << "Starting OOF Baking...\n";
    std::cout << "Mesh Triangles: " << (mesh.indices.size() / 3) << "\n";
    std::cout << "Grid: " << resolution << "^3, Rays per voxel: " << num_rays << "\n";

    // 1. Calculate AABB of the mesh
    Eigen::Vector3f min_bound = mesh.vertices.empty() ? Eigen::Vector3f::Zero() : mesh.vertices[0].position;
    Eigen::Vector3f max_bound = min_bound;
    
    for(const auto& v : mesh.vertices) {
        min_bound = min_bound.cwiseMin(v.position);
        max_bound = max_bound.cwiseMax(v.position);
    }
    
    // Pad the bounds slightly so the mesh is completely inside the volume
    Eigen::Vector3f padding = (max_bound - min_bound) * 0.1f;
    min_bound -= padding;
    max_bound += padding;
    
    std::cout << "OOF Bounds: Min(" << min_bound.transpose() << ") Max(" << max_bound.transpose() << ")\n";

    auto oof = std::make_unique<OOF>(resolution, resolution, resolution, min_bound, max_bound, visibility_order);
    
    // 2. Generate Fibonacci sphere directions
    std::vector<Eigen::Vector3f> ray_dirs;
    std::vector<std::pair<double, double>> ray_angles; // theta, phi for SH eval
    ray_dirs.reserve(num_rays);
    ray_angles.reserve(num_rays);
    
    const double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
    for (int i = 0; i < num_rays; ++i) {
        double theta = 2.0 * std::numbers::pi * i / golden_ratio; // Azimuthal
        double phi = std::acos(1.0 - 2.0 * (i + 0.5) / num_rays); // Polar
        
        float x = static_cast<float>(std::sin(phi) * std::cos(theta));
        float y = static_cast<float>(std::cos(phi));
        float z = static_cast<float>(std::sin(phi) * std::sin(theta));
        
        ray_dirs.push_back(Eigen::Vector3f(x, y, z));
        ray_angles.push_back({phi, theta}); // Note: Sylph SH uses theta=polar, phi=azimuthal
    }
    
    sylph::SphericalHarmonics sh;
    Eigen::Vector3f step = (max_bound - min_bound) / static_cast<float>(resolution);
    
    // 3. Raytrace every voxel
    int total_voxels = resolution * resolution * resolution;
    int baked_count = 0;
    
    for (int z = 0; z < resolution; ++z) {
        for (int y = 0; y < resolution; ++y) {
            for (int x = 0; x < resolution; ++x) {
                // Voxel center position
                Eigen::Vector3f pos = min_bound + Eigen::Vector3f(
                    (x + 0.5f) * step.x(),
                    (y + 0.5f) * step.y(),
                    (z + 0.5f) * step.z()
                );
                
                std::vector<double> sh_coeffs(visibility_order * visibility_order, 0.0);
                
                // Cast rays
                for(int r = 0; r < num_rays; ++r) {
                    const Eigen::Vector3f& dir = ray_dirs[r];
                    bool hit = false;
                    
                    // Check intersection against all triangles
                    for(size_t i = 0; i < mesh.indices.size(); i += 3) {
                        const Eigen::Vector3f& v0 = mesh.vertices[mesh.indices[i]].position;
                        const Eigen::Vector3f& v1 = mesh.vertices[mesh.indices[i+1]].position;
                        const Eigen::Vector3f& v2 = mesh.vertices[mesh.indices[i+2]].position;
                        
                        if (ray_triangle_intersect(pos, dir, v0, v1, v2)) {
                            hit = true;
                            break;
                        }
                    }
                    
                    double visibility = hit ? 0.0 : 1.0;
                    
                    // Accumulate SH Projection
                    // Integral(V(w) * Y_lm(w)) dw  ~~ (4PI / N) * Sum(V_i * Y_lm(w_i))
                    double weight = (4.0 * std::numbers::pi) / num_rays;
                    
                    double polar = ray_angles[r].first;
                    double azimuthal = ray_angles[r].second;
                    
                    int sh_idx = 0;
                    for(int l = 0; l < visibility_order; ++l) {
                        for(int m = -l; m <= l; ++m) {
                            sh_coeffs[sh_idx++] += visibility * sh.evaluate(l, m, polar, azimuthal) * weight;
                        }
                    }
                }
                
                // Write to OOF
                auto& vis = oof->at(x, y, z);
                int sh_idx = 0;
                for(int l = 0; l < visibility_order; ++l) {
                    for(int m = -l; m <= l; ++m) {
                        vis.coefficients()(l, m) = sh_coeffs[sh_idx++];
                    }
                }
                
                baked_count++;
                if (baked_count % (total_voxels / 10) == 0) {
                    std::cout << "Baking... " << (100 * baked_count / total_voxels) << "%\n";
                }
            }
        }
    }
    
    std::cout << "OOF Baking Complete!\n";
    return oof;
}

} // namespace render
