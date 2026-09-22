#pragma once

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/sh_lobe_mesh.hpp"
#include "render/transform.hpp"

#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"

#include <Eigen/Dense>

#include <SDL3/SDL.h>

namespace examples
{

/**
 * @brief Interactive Spherical Harmonics Product Visualizer.
 * 
 * This example visually demonstrates the core mathematical operation that powers 
 * the entire "Fast and Accurate Spherical Harmonics Products" paper: the SH Product.
 * 
 * WHAT IT DOES:
 * - It creates two distinct Spherical Harmonic shapes (functions over a sphere):
 *   1. Shape A: A dynamic, rotating directional lobe.
 *   2. Shape B: A static quadrupole (pinched) lobe.
 * - Every frame, it uses the paper's 'SHGridTransform::product' algorithm to 
 *   multiply Shape A and Shape B together.
 * - It dynamically generates and renders a 3D mesh (LobeMesh) representing the 
 *   resulting SH coefficients.
 * 
 * CONTROLS:
 * Press 1: View Shape A (Moving Directional Lobe)
 * Press 2: View Shape B (Static Quadrupole Lobe)
 * Press 3: View the Product (A x B) computed using the Spherical Grid algorithm.
 * 
 * WHY THIS MATTERS:
 * In Glossy Relighting (Section 7.3), Shape A represents the incoming Light, 
 * Shape B represents the BRDF, and the Product represents the final masked 
 * lighting lobe before integration. This visualizer proves that the paper's 
 * Grid Transform successfully multiplies two complex spherical functions 
 * interactively in real-time.
 */
class SHLobeExample
{
public:
    SHLobeExample(
        int width,
        int height
    );

    void process_event(
        const SDL_Event& event,
        SDL_Window* window
    );

    void update(float dt);

    void resize(
        int width,
        int height
    );

    void render(
        render::Renderer& renderer
    );

private:
    sylph::SHCoefficients coefficients_;

    render::SHLobeMesh lobe_;

    render::Camera camera_;
    render::CameraController camera_controller_;

    render::Shader shader_;

    render::Transform transform_;

 private:
    sylph::SHCoefficients coeff_a_;
    sylph::SHCoefficients coeff_b_;
    sylph::SHGridTransform sh_transform_;
    
    int current_view_{3}; // 1 = A, 2 = B, 3 = Product
    float time_{0.0f};   
};

} // namespace examples