#include "product_demo_example.hpp"
#include "sylph/spherical_harmonics.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace examples {

constexpr const char* PRODUCT_VERT = R"(
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

constexpr const char* PRODUCT_FRAG = R"(
#version 330 core
in vec3 v_normal;
in vec2 v_uv;
out vec4 FragColor;
void main() {
    vec3 n = normalize(v_normal);
    vec3 light_dir = normalize(vec3(0.0, 1.0, 1.0));
    float diff = max(dot(n, light_dir), 0.3);
    vec3 color = (v_uv.x > 0.5) ? vec3(0.2, 0.8, 0.2) : vec3(0.8, 0.2, 0.2);
    FragColor = vec4(color * diff, 1.0);
}
)";

ProductDemoExample::ProductDemoExample(int width, int height) {
    std::cout << "Loading Product Demo...\n";
    std::cout << "Controls: [T/G] change Grid Truncation (t)\n";
    
    shader_ = std::make_unique<render::Shader>(PRODUCT_VERT, PRODUCT_FRAG);
    
    camera_.set_position(Eigen::Vector3f(0.0f, 0.0f, 10.0f));
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
    
    transform_a_.set_position(Eigen::Vector3f(-4.0f, 0.0f, 0.0f));
    transform_b_.set_position(Eigen::Vector3f(0.0f, 0.0f, 0.0f));
    transform_c_.set_position(Eigen::Vector3f(4.0f, 0.0f, 0.0f));

    // B is a static quadrupole
    coeff_b_(2, 0) = 1.0;
    lobe_b_ = std::make_unique<render::SHLobeMesh>(coeff_b_);
}

void ProductDemoExample::update_meshes() {
    // A rotates
    for(int l=0; l<8; ++l) for(int m=-l; m<=l; ++m) coeff_a_(l, m) = 0.0f;
    float theta = std::numbers::pi_v<float> / 2.0f;
    float phi = time_ * 1.5f;
    sylph::SphericalHarmonics sh;
    for(int l=0; l<8; ++l) {
        for(int m=-l; m<=l; ++m) {
            coeff_a_(l, m) = sh.evaluate(l, m, theta, phi) * 2.0; // Cosine-like lobe approximation
        }
    }
    
    lobe_a_ = std::make_unique<render::SHLobeMesh>(coeff_a_);
    
    // Product
    sylph::SHGridTransform transform;
    auto coeff_c = transform.product(coeff_a_, coeff_b_, 8, t_param_);
    lobe_product_ = std::make_unique<render::SHLobeMesh>(coeff_c);
}

void ProductDemoExample::process_event(const SDL_Event& event, SDL_Window* window) {
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
    
    if (event.type == SDL_EVENT_KEY_DOWN) {
        if (event.key.key == SDLK_T) { t_param_ = std::min(t_param_ + 1, 8); std::cout << "t=" << t_param_ << "\n"; }
        if (event.key.key == SDLK_G) { t_param_ = std::max(t_param_ - 1, 1); std::cout << "t=" << t_param_ << "\n"; }
    }
}

void ProductDemoExample::update(float dt) {
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, dt);
    time_ += dt;
    update_meshes();
}

void ProductDemoExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void ProductDemoExample::render(render::Renderer& renderer) {
    (void)renderer;
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    
    shader_->bind();
    shader_->set_mat4("uView", camera_.view_matrix());
    shader_->set_mat4("uProjection", camera_.projection_matrix());
    
    shader_->set_mat4("uModel", transform_a_.model_matrix());
    if (lobe_a_) lobe_a_->mesh().draw();
    
    shader_->set_mat4("uModel", transform_b_.model_matrix());
    if (lobe_b_) lobe_b_->mesh().draw();
    
    shader_->set_mat4("uModel", transform_c_.model_matrix());
    if (lobe_product_) lobe_product_->mesh().draw();
}

}

