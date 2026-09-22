#ifndef SYLPH_EXAMPLES_OOF_EXAMPLE_HPP
#define SYLPH_EXAMPLES_OOF_EXAMPLE_HPP

#include "render/camera_controller.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/lighting.hpp"
#include "render/color_mesh.hpp"
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
    void build_dense_floor();

    render::Camera camera_;
    render::CameraController camera_controller_;
    
    render::Shader shader_;
    render::Lighting lighting_;
    // Environment
    struct SceneObject {
        std::shared_ptr<render::Mesh> mesh;
        render::Transform transform;
        Eigen::Vector3f base_color;
    };
    std::vector<SceneObject> house_objects_;
    
    // Dynamic Occluder
    std::shared_ptr<render::Mesh> occluder_mesh_;
    render::Transform occluder_transform_;
    std::unique_ptr<render::OOF> occluder_oof_;
    
    // Receiver (Dense Ground Plane)
    std::shared_ptr<render::ColorMesh> floor_mesh_;
    std::vector<render::ColorVertex> floor_vertices_; // CPU side copy for dynamic updates
    render::Transform floor_transform_;
    render::Visibility floor_self_visibility_;
    
    render::ShadowField shadow_field_;
    
    float time_ = 0.0f;
};

} // namespace examples

#endif
