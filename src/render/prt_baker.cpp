#include "render/prt_baker.hpp"
#include "sylph/spherical_harmonics.hpp"
#include <iostream>
#include <cmath>
#include <numbers>
#include <limits>
#include <vector>

namespace render {

static bool ray_triangle_intersect(
    const Eigen::Vector3f& orig, const Eigen::Vector3f& dir,
    const Eigen::Vector3f& v0, const Eigen::Vector3f& v1, const Eigen::Vector3f& v2)
{
    const float EPSILON = 1e-5f;
    Eigen::Vector3f edge1 = v1 - v0;
    Eigen::Vector3f edge2 = v2 - v0;
    Eigen::Vector3f h = dir.cross(edge2);
    float a = edge1.dot(h);
    if (a > -EPSILON && a < EPSILON) return false;
    float f = 1.0f / a;
    Eigen::Vector3f s = orig - v0;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) return false;
    Eigen::Vector3f q = s.cross(edge1);
    float v = f * dir.dot(q);
    if (v < 0.0f || u + v > 1.0f) return false;
    float t = f * edge2.dot(q);
    return t > EPSILON;
}

std::shared_ptr<PRTMesh> PRTBaker::bake(const MeshData& mesh, int order, int num_rays, bool shadows)
{
    std::cout << "Baking PRT Mesh... Vertices: " << mesh.vertices.size() << ", Rays/Vert: " << num_rays << "\n";
    
    std::vector<PRTVertex> prt_verts;
    prt_verts.reserve(mesh.vertices.size());
    
    sylph::SphericalHarmonics sh;
    int num_coeffs = order * order;
    
    // Generate uniform ray directions on a sphere using Fibonacci lattice
    std::vector<Eigen::Vector3f> ray_dirs;
    std::vector<std::vector<double>> ray_sh_evals;
    ray_dirs.reserve(num_rays);
    ray_sh_evals.reserve(num_rays);
    
    const double phi_golden = std::numbers::pi_v<double> * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < num_rays; ++i) {
        double y = 1.0 - (i / float(num_rays - 1)) * 2.0;
        double radius = std::sqrt(1.0 - y * y);
        double theta_spiral = phi_golden * i;
        double x = std::cos(theta_spiral) * radius;
        double z = std::sin(theta_spiral) * radius;
        
        Eigen::Vector3f dir(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
        ray_dirs.push_back(dir.normalized());
        
        double sh_theta = std::acos(y);
        double sh_phi = std::atan2(z, x);
        if (sh_phi < 0) sh_phi += 2.0 * std::numbers::pi_v<double>;
        
        std::vector<double> evals(num_coeffs, 0.0);
        int idx = 0;
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                evals[idx++] = sh.evaluate(l, m, sh_theta, sh_phi);
            }
        }
        ray_sh_evals.push_back(evals);
    }
    
    float mc_weight = (4.0f * std::numbers::pi_v<float>) / num_rays;
    
    int vert_idx = 0;
    for (const auto& v : mesh.vertices) {
        PRTVertex prt_v;
        prt_v.position = v.position;
        prt_v.normal = v.normal.normalized();
        for (int c = 0; c < 16; ++c) prt_v.sh[c] = 0.0f;
        
        Eigen::Vector3f orig = prt_v.position + prt_v.normal * 1e-4f;
        
        std::vector<double> coeffs(num_coeffs, 0.0);
        
        for (int r = 0; r < num_rays; ++r) {
            const auto& dir = ray_dirs[r];
            float ndotl = prt_v.normal.dot(dir);
            
            if (ndotl <= 0.0f) continue;
            
            float visibility = 1.0f;
            if (shadows) {
                for (size_t i = 0; i < mesh.indices.size(); i += 3) {
                    const auto& v0 = mesh.vertices[mesh.indices[i+0]].position;
                    const auto& v1 = mesh.vertices[mesh.indices[i+1]].position;
                    const auto& v2 = mesh.vertices[mesh.indices[i+2]].position;
                    
                    if (ray_triangle_intersect(orig, dir, v0, v1, v2)) {
                        visibility = 0.0f;
                        break;
                    }
                }
            }
            
            if (visibility > 0.0f) {
                for (int c = 0; c < num_coeffs; ++c) {
                    coeffs[c] += ray_sh_evals[r][c] * ndotl * mc_weight;
                }
            }
        }
        
        for (int c = 0; c < std::min(num_coeffs, 16); ++c) {
            prt_v.sh[c] = static_cast<float>(coeffs[c]);
        }
        
        prt_verts.push_back(prt_v);
        
        vert_idx++;
        if (vert_idx % 250 == 0) {
            std::cout << "Baked " << vert_idx << " / " << mesh.vertices.size() << " vertices...\r";
        }
    }
    std::cout << "\nPRT Bake Complete!\n";
    
    return std::make_shared<PRTMesh>(prt_verts, mesh.indices);
}

}

