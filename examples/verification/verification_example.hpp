#pragma once

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/transform.hpp"

#include <SDL3/SDL.h>

namespace examples
{

class VerificationExample
{
public:
    VerificationExample(int width, int height);

    void process_event(const SDL_Event& event, SDL_Window* window);
    void update(float dt);
    void resize(int width, int height);
    void render(render::Renderer& renderer);

private:
    void run_verification();

    render::Camera camera_;
    render::CameraController camera_controller_;
    bool run_once_{false};
};

} // namespace examples
