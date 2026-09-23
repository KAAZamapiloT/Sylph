#pragma once

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/sh_lobe_mesh.hpp"
#include "render/transform.hpp"
#include <SDL3/SDL.h>
#include <memory>

namespace examples {

class SHBasisExample {
public:
    SHBasisExample(int width, int height);

    void process_event(const SDL_Event& event, SDL_Window* window);
    void update(float dt);
    void resize(int width, int height);
    void render(render::Renderer& renderer);

private:
    void update_mesh();

    render::Camera camera_;
    render::CameraController camera_controller_;
    std::unique_ptr<render::Shader> shader_;
    
    std::unique_ptr<render::SHLobeMesh> lobe_mesh_;
    render::Transform transform_;
    
    int l_{2};
    int m_{0};
    float time_{0.0f};
};

} // namespace examples

