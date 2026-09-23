#ifndef SYLPH_EXAMPLES_OOF_EXAMPLE_HPP
#define SYLPH_EXAMPLES_OOF_EXAMPLE_HPP

#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/lighting.hpp"
#include "render/color_mesh.hpp"
#include "render/mesh.hpp"
#include "render/transform.hpp"
#include "render/oof_baker.hpp"
#include "render/shadow_field.hpp"
#include "sylph/sh_coefficients.hpp"

#include <vector>
#include <memory>
#include <SDL3/SDL.h>

namespace examples {

class OofExample {
public:
    OofExample(int width, int height);

    void process_event(const SDL_Event& event, SDL_Window* window);
    void update(float delta_time);
    void render(render::Renderer& renderer);
    void resize(int width, int height);

private:
    render::Camera camera_;
    render::CameraController camera_controller_;
    
    render::Shader shader_;
    render::Lighting lighting_;
    
    // Dynamic Occluder
    std::shared_ptr<render::Mesh> occluder_mesh_;
    render::Transform occluder_transform_;
    std::unique_ptr<render::OOF> occluder_oof_;
    
    // Receiver Room
    std::shared_ptr<render::ColorMesh> room_mesh_;
    std::vector<render::ColorVertex> room_vertices_; 
    std::vector<Eigen::Vector3f> room_normals_;
    std::vector<render::Visibility> room_self_visibility_;
    
    render::ShadowField shadow_field_;
    
    float time_ = 0.0f;
};

} // namespace examples

#endif

