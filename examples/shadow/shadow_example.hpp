#pragma once

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/transform.hpp"
#include "render/prt_mesh.hpp"
#include "render/lighting.hpp"
#include <Eigen/Dense>
#include <SDL3/SDL.h>
#include <vector>
#include <memory>

namespace examples {

/**
 * @brief Static Shadowed Diffuse PRT Demo.
 * Demonstrates Precomputed Radiance Transfer by baking a transfer function 
 * (which includes visibility/shadows) into the mesh vertices offline.
 */
class ShadowExample {
public:
    ShadowExample(int width, int height);

    void process_event(const SDL_Event& event, SDL_Window* window);
    void update(float dt);
    void resize(int width, int height);
    void render(render::Renderer& renderer);

private:
    render::Camera camera_;
    render::CameraController camera_controller_;
    render::Shader shader_;
    render::Lighting lighting_;
    
    std::shared_ptr<render::PRTMesh> shadowed_mesh_;
    std::shared_ptr<render::PRTMesh> unshadowed_mesh_;
    render::Transform transform_;
    Eigen::Vector3f base_color_;
    
    bool show_shadows_{true};
    float time_{0.0f};
};

} // namespace examples

