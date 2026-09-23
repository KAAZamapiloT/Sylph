#include "shadow_example.hpp"
#include "render/obj_loader.hpp"
#include "render/prt_baker.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace examples {

static std::string load_shader_source(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) return "";
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

ShadowExample::ShadowExample(int width, int height)
    : shader_(load_shader_source("shaders/prt_shadow.vert"), load_shader_source("shaders/prt_shadow.frag")),
      lighting_(4)
{
    std::cout << "Loading Static Shadow PRT Demo...\n";
    
    render::MeshData mesh_data = render::ObjLoader::load("assets/meshes/teapot.obj");
    Eigen::Vector3f min_b = mesh_data.vertices[0].position;
    Eigen::Vector3f max_b = mesh_data.vertices[0].position;
    for (const auto& v : mesh_data.vertices) {
        min_b = min_b.cwiseMin(v.position); max_b = max_b.cwiseMax(v.position);
    }
    Eigen::Vector3f center = (min_b + max_b) * 0.5f;
    float extent = (max_b - min_b).maxCoeff();
    for (auto& v : mesh_data.vertices) {
        v.position = (v.position - center) * (2.0f / extent);
    }

    std::cout << "Baking Unshadowed PRT (Instant)...\n";
    unshadowed_mesh_ = render::PRTBaker::bake(mesh_data, 4, 50, false);
    
    std::cout << "Baking Shadowed PRT (50 rays, may take 5-10s)...\n";
    shadowed_mesh_ = render::PRTBaker::bake(mesh_data, 4, 50, true);

    transform_.set_position(Eigen::Vector3f(0.0f, -0.5f, 0.0f));
    transform_.set_scale(Eigen::Vector3f(2.0f, 2.0f, 2.0f));
    base_color_ = Eigen::Vector3f(1.0f, 0.8f, 0.2f); // Gold

    camera_.set_position(Eigen::Vector3f(0.0f, 0.5f, 4.0f));
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
    
    camera_controller_.set_move_speed(8.0f);

    auto& light_r = lighting_.coefficients(render::LightChannel::Red);
    auto& light_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_b = lighting_.coefficients(render::LightChannel::Blue);
    
    // Key light from top-right
    light_r(0, 0) = 0.5; light_g(0, 0) = 0.5; light_b(0, 0) = 0.6;
    light_r(1, 0) = 1.0; light_g(1, 0) = 1.0; light_b(1, 0) = 0.9;
    light_r(1, 1) = 0.8; light_g(1, 1) = 0.8; light_b(1, 1) = 0.7;
    
    std::cout << "\n[ShadowExample] Ready! Press 'S' to toggle shadows.\n";
}

void ShadowExample::process_event(const SDL_Event& event, SDL_Window* window) {
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
    
    if (event.type == SDL_EVENT_KEY_DOWN && event.key.key == SDLK_S) {
        show_shadows_ = !show_shadows_;
        std::cout << "Shadows " << (show_shadows_ ? "ON" : "OFF") << "\n";
    }
}

void ShadowExample::update(float dt) {
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, dt);
    time_ += dt;
    transform_.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(time_ * 0.5f, Eigen::Vector3f::UnitY())));
}

void ShadowExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void ShadowExample::render(render::Renderer& renderer) {
    (void)renderer;
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    
    shader_.bind();
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix());
    
    float light_coeffs[48]; 
    auto& light_r = lighting_.coefficients(render::LightChannel::Red);
    auto& light_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_b = lighting_.coefficients(render::LightChannel::Blue);
    int idx = 0;
    for (int l = 0; l <= 3; ++l) {
        for (int m = -l; m <= l; ++m) {
            light_coeffs[idx*3 + 0] = static_cast<float>(light_r(l, m));
            light_coeffs[idx*3 + 1] = static_cast<float>(light_g(l, m));
            light_coeffs[idx*3 + 2] = static_cast<float>(light_b(l, m));
            idx++;
        }
    }
    GLint loc = glGetUniformLocation(shader_.id(), "uLightSH");
    glUniform3fv(loc, 16, light_coeffs);

    shader_.set_mat4("uModel", transform_.model_matrix());
    shader_.set_vec3("uBaseColor", base_color_);
    
    if (show_shadows_) {
        shadowed_mesh_->draw();
    } else {
        unshadowed_mesh_->draw();
    }
}

}

