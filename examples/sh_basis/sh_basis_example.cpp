#include "sh_basis_example.hpp"
#include "sylph/sh_coefficients.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace examples {

constexpr const char* SH_BASIS_VERT = R"(
#version 330 core
layout(location = 0) in vec3 a_position;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec2 a_uv;
uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
out vec3 v_normal;
out vec2 v_uv;
void main() {
    v_normal = mat3(uModel) * a_normal;
    v_uv = a_uv;
    gl_Position = uProjection * uView * uModel * vec4(a_position, 1.0);
}
)";

constexpr const char* SH_BASIS_FRAG = R"(
#version 330 core
in vec3 v_normal;
in vec2 v_uv;
out vec4 FragColor;
void main() {
    vec3 n = normalize(v_normal);
    vec3 light_dir = normalize(vec3(1.0, 1.0, 1.0));
    float diff = max(dot(n, light_dir), 0.3);
    vec3 color = (v_uv.x > 0.5) ? vec3(0.2, 0.8, 0.2) : vec3(0.8, 0.2, 0.2);
    FragColor = vec4(color * diff, 1.0);
}
)";

SHBasisExample::SHBasisExample(int width, int height) {
    std::cout << "Loading SH Basis Visualizer...\n";
    std::cout << "Controls: [L/K] change Band (l), [M/N] change Degree (m)\n";
    
    shader_ = std::make_unique<render::Shader>(SH_BASIS_VERT, SH_BASIS_FRAG);
    
    camera_.set_position(Eigen::Vector3f(0.0f, 0.0f, 3.0f));
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
    
    update_mesh();
}

void SHBasisExample::update_mesh() {
    sylph::SHCoefficients coeffs(8); 
    coeffs(l_, m_) = 1.0;
    lobe_mesh_ = std::make_unique<render::SHLobeMesh>(coeffs);
}

void SHBasisExample::process_event(const SDL_Event& event, SDL_Window* window) {
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
    
    if (event.type == SDL_EVENT_KEY_DOWN) {
        bool changed = false;
        if (event.key.key == SDLK_L) { l_ = std::min(l_ + 1, 7); changed = true; }
        if (event.key.key == SDLK_K) { l_ = std::max(l_ - 1, 0); changed = true; }
        
        if (changed) m_ = std::clamp(m_, -l_, l_);
        
        if (event.key.key == SDLK_M) { m_ = std::min(m_ + 1, l_); changed = true; }
        if (event.key.key == SDLK_N) { m_ = std::max(m_ - 1, -l_); changed = true; }
        
        if (changed) {
            std::cout << "SH Basis: l=" << l_ << ", m=" << m_ << "\n";
            update_mesh();
        }
    }
}

void SHBasisExample::update(float dt) {
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, dt);
    time_ += dt;
    transform_.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(time_ * 0.5f, Eigen::Vector3f::UnitY())));
}

void SHBasisExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void SHBasisExample::render(render::Renderer& renderer) {
    (void)renderer;
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    
    shader_->bind();
    shader_->set_mat4("uView", camera_.view_matrix());
    shader_->set_mat4("uProjection", camera_.projection_matrix());
    shader_->set_mat4("uModel", transform_.model_matrix());
    
    lobe_mesh_->mesh().draw();
}

}


