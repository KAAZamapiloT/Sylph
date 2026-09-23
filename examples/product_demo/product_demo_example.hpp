#pragma once

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/sh_lobe_mesh.hpp"
#include "render/transform.hpp"
#include "sylph/sh_grid_transform.hpp"
#include <SDL3/SDL.h>
#include <memory>

namespace examples {

class ProductDemoExample {
public:
    ProductDemoExample(int width, int height);

    void process_event(const SDL_Event& event, SDL_Window* window);
    void update(float dt);
    void resize(int width, int height);
    void render(render::Renderer& renderer);

private:
    void update_meshes();

    render::Camera camera_;
    render::CameraController camera_controller_;
    std::unique_ptr<render::Shader> shader_;
    
    std::unique_ptr<render::SHLobeMesh> lobe_a_;
    std::unique_ptr<render::SHLobeMesh> lobe_b_;
    std::unique_ptr<render::SHLobeMesh> lobe_product_;
    
    render::Transform transform_a_;
    render::Transform transform_b_;
    render::Transform transform_c_;
    
    sylph::SHCoefficients coeff_a_{8};
    sylph::SHCoefficients coeff_b_{8};
    
    int t_param_{8};
    float time_{0.0f};
};

} // namespace examples
