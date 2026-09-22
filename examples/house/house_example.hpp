#pragma once

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/transform.hpp"
#include "render/color_mesh.hpp"
#include "render/lighting.hpp"
#include "render/brdf.hpp"
#include "render/visibility.hpp"
#include "render/glossy_relighting.hpp"
#include "sylph/sh_coefficients.hpp"

#include <Eigen/Dense>
#include <SDL3/SDL.h>
#include <vector>
#include <memory>

namespace examples
{

/**
 * @brief GPU-Accelerated Glossy PRT Demo using the True Spherical Grid Algorithm.
 * 
 * This example correctly implements Section 7.3.1 ("Unshadowed / Shadowed Glossy PRT") 
 * of the paper "Fast and Accurate Spherical Harmonics Products".
 * 
 * THE PROBLEM:
 * To render glossy materials with shadows under dynamic lighting, we must evaluate 
 * the Triple Product Integral for every pixel:  Integral( L * B * V ) d(omega)
 * where:
 *   L = Incident Lighting Environment (Dynamic, View-independent)
 *   B = BRDF (Dynamic, View-dependent, centered around reflection vector R)
 *   V = Visibility / Self-Shadowing (Static, View-independent, baked per-vertex)
 * 
 * Evaluating this integral using traditional Clebsch-Gordan SH coefficients is 
 * extraordinarily slow.
 * 
 * THE PAPER'S SOLUTION (Spherical Grids):
 * Instead of computing the product in SH frequency space, we transform the functions 
 * onto a carefully chosen "Spherical Grid" (Gauss-Legendre nodes for theta, uniform for phi), 
 * do a simple point-wise multiplication, and sum the results using quadrature weights 
 * to instantly extract the DC (Direct Current) component, which perfectly equals the integral!
 * 
 * OUR HYBRID CPU-GPU ARCHITECTURE:
 * To achieve thousands of FPS, we split the Grid algorithm across the CPU and GPU:
 * 
 * 1. CPU (Grid Generation & Caching):
 *    - Because 'L' and 'V' do not depend on the camera, the CPU evaluates them on the 
 *      Spherical Grid exactly once per frame (or at load time for V).
 *    - The grid size is dynamically determined by the 't' parameter (t_param_). 
 *      Grid Size = (t+1)/2 * t. For t=8, the grid is only 32 to 128 points!
 *    - We send L_grid, V_grid, Grid_Directions, and Grid_Weights to the GPU as uniform arrays.
 * 
 * 2. GPU (Fragment Shader - Per Pixel):
 *    - The GPU calculates the Reflection vector 'R' based on the camera view.
 *    - It loops through the small uniform array of Grid points.
 *    - For each point, it dynamically evaluates the BRDF (B_grid) using an ultra-fast 
 *      Legendre polynomial recurrence relation.
 *    - Finally, it computes the point-wise product: L_grid * B_grid * V_grid * Weight
 *    - The sum of these products across the grid is the exact Triple Product Integral!
 */
class HouseExample
{
public:
    HouseExample(int width, int height);

    void process_event(const SDL_Event& event, SDL_Window* window);
    void update(float dt);
    void resize(int width, int height);
    void render(render::Renderer& renderer);

private:
    void generate_lighting(float time);
    void update_roughness();
    void update_grid_data();

    render::Camera camera_;
    render::CameraController camera_controller_;
    render::Shader shader_;

    struct SceneObject {
        std::shared_ptr<render::Mesh> mesh;
        render::Transform transform;
        float roughness = 1.0f;
        Eigen::Vector3f base_color = Eigen::Vector3f(1.0f, 1.0f, 1.0f);
        
        double base_brdf(int l, int m) const {
            if (m != 0) return 0.0;
            double factor = std::exp(-roughness * l * (l + 1));
            return factor * std::sqrt((2 * l + 1) / (4.0 * std::numbers::pi_v<double>));
        }
    };
    std::vector<SceneObject> objects_;

    render::Lighting lighting_;
    
    // Grid Data for GPU
    std::vector<Eigen::Vector3f> grid_light_;
    std::vector<float> grid_weights_;
    std::vector<Eigen::Vector3f> grid_dirs_;
    std::vector<float> grid_vis_;
    
    float time_{0.0f};
    float roughness_{20.0f}; // Phong exponent
    int t_param_{8};         // Truncation parameter t (1-8)
};

} // namespace examples
