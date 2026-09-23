#include "oof_example.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "render/obj_loader.hpp"
#include "render/oof_baker.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <numbers>

namespace examples {

static std::string load_shader_source(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) return "";
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

OofExample::OofExample(int width, int height)
    : shader_(load_shader_source("shaders/oof.vert"), load_shader_source("shaders/oof.frag")),
      lighting_(4),
      shadow_field_()
{
    std::cout << "Loading OOF Shadow Fields Demo...\n";
    std::cout << "Preparing procedural room and dynamic occluder...\n";
    
    render::MeshData mesh_data = render::ObjLoader::load("assets/meshes/bunny.obj");
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

    occluder_mesh_ = std::make_shared<render::Mesh>(mesh_data.vertices, mesh_data.indices);
    occluder_oof_ = render::OOFBaker::bake(mesh_data, 8, 4, 100);
    
    // Procedural Room (Floor + 2 Walls)
    const int segments = 20; 
    const float size = 10.0f;
    const float step = size / segments;
    
    std::vector<uint32_t> indices;
    
    auto add_plane = [&](Eigen::Vector3f origin, Eigen::Vector3f right, Eigen::Vector3f up, Eigen::Vector3f normal) {
        uint32_t base_idx = room_vertices_.size();
        for (int v = 0; v <= segments; ++v) {
            for (int u = 0; u <= segments; ++u) {
                render::ColorVertex vert;
                vert.position = origin + right * (u * step) + up * (v * step);
                vert.color = Eigen::Vector3f(0.8f, 0.8f, 0.8f);
                room_vertices_.push_back(vert);
                room_normals_.push_back(normal);
            }
        }
        for (int v = 0; v < segments; ++v) {
            for (int u = 0; u < segments; ++u) {
                uint32_t top_left = base_idx + v * (segments + 1) + u;
                uint32_t top_right = top_left + 1;
                uint32_t bottom_left = base_idx + (v + 1) * (segments + 1) + u;
                uint32_t bottom_right = bottom_left + 1;
                indices.push_back(top_left); indices.push_back(bottom_left); indices.push_back(top_right);
                indices.push_back(top_right); indices.push_back(bottom_left); indices.push_back(bottom_right);
            }
        }
    };
    
    // Floor
    add_plane(Eigen::Vector3f(-size/2, 0.0f, -size/2), Eigen::Vector3f(1.0f, 0, 0), Eigen::Vector3f(0, 0, 1.0f), Eigen::Vector3f(0, 1, 0));
    // Back Wall
    add_plane(Eigen::Vector3f(-size/2, 0.0f, size/2), Eigen::Vector3f(1.0f, 0, 0), Eigen::Vector3f(0, 1.0f, 0), Eigen::Vector3f(0, 0, -1));
    // Left Wall
    add_plane(Eigen::Vector3f(-size/2, 0.0f, -size/2), Eigen::Vector3f(0, 0, 1.0f), Eigen::Vector3f(0, 1.0f, 0), Eigen::Vector3f(1, 0, 0));

    room_mesh_ = std::make_shared<render::ColorMesh>(room_vertices_, indices);
    
    // Evaluate self-visibility for room (unshadowed = 1.0 everywhere)
    for (size_t i = 0; i < room_vertices_.size(); ++i) { room_self_visibility_.emplace_back(4); }
    sylph::SphericalHarmonics sh;
    for (auto& v : room_self_visibility_) {
        sylph::SHCoefficients c(4);
        c(0, 0) = std::sqrt(4.0f * std::numbers::pi_v<float>); // Constant 1.0 function
        v.set_coefficients(c);
    }

    camera_.set_position(Eigen::Vector3f(5.0f, 4.0f, -8.0f));
    
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
    
    camera_controller_.set_move_speed(8.0f);

    auto& light_r = lighting_.coefficients(render::LightChannel::Red);
    auto& light_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_b = lighting_.coefficients(render::LightChannel::Blue);
    
    // Directional-ish light from top right
    float theta = std::numbers::pi_v<float> / 4.0f;
    float phi = 0.5f;
    for (int l = 0; l < 4; ++l) {
        for (int m = -l; m <= l; ++m) {
            float val = sh.evaluate(l, m, theta, phi) * 2.0f;
            light_r(l, m) = val; light_g(l, m) = val * 0.9f; light_b(l, m) = val * 0.8f;
        }
    }
    light_r(0, 0) += 0.5f; light_g(0, 0) += 0.5f; light_b(0, 0) += 0.5f; // ambient
}

void OofExample::process_event(const SDL_Event& event, SDL_Window* window) {
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
}

void OofExample::update(float delta_time) {
    time_ += delta_time;
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, delta_time);
    
    // Animate bunny in a circle near the walls
    float ox = std::sin(time_ * 0.8f) * 3.0f;
    float oz = std::cos(time_ * 0.8f) * 3.0f;
    float oy = 1.0f + std::sin(time_ * 2.0f) * 0.5f;
    occluder_transform_.set_position(Eigen::Vector3f(ox, oy, oz));
    
    // Update OOF shadows on the room vertices
    std::vector<const render::OOF*> active_oofs = { occluder_oof_.get() };
    Eigen::Matrix4f occluder_inv = occluder_transform_.model_matrix().inverse();
    
    for (size_t i = 0; i < room_vertices_.size(); ++i) {
        auto& v = room_vertices_[i];
        Eigen::Vector4f world_pos(v.position.x(), v.position.y(), v.position.z(), 1.0f);
        Eigen::Vector4f local_pos = occluder_inv * world_pos;
        
        Eigen::Vector3f shadowed_color;
        try {
            shadowed_color = shadow_field_.evaluate(
                lighting_, room_self_visibility_[i], local_pos.head<3>(), active_oofs);
        } catch (const std::out_of_range&) {
            shadowed_color = shadow_field_.evaluate(
                lighting_, room_self_visibility_[i], local_pos.head<3>(), {});
        }
        
        // Simple dot product diffuse for the room
        Eigen::Vector3f n = room_normals_[i];
        Eigen::Vector3f l_dir = Eigen::Vector3f(1.0f, 1.0f, 1.0f).normalized();
        float diff = std::max(0.0f, n.dot(l_dir));
        
        Eigen::Vector3f base_color = Eigen::Vector3f(0.8f, 0.8f, 0.8f); 
        v.color = base_color.cwiseProduct(shadowed_color) * diff;
    }
    room_mesh_->update_colors(room_vertices_);
}

void OofExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void OofExample::render(render::Renderer& renderer) {
    (void)renderer;
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    
    shader_.bind();
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix());
    
    // Draw Room (CPU shadow evaluated)
    shader_.set_mat4("uModel", Eigen::Matrix4f::Identity());
    shader_.set_int("uUseShading", 0);
    room_mesh_->draw();
    
    // Draw Occluder (GPU diffuse SH)
    shader_.set_mat4("uModel", occluder_transform_.model_matrix());
    shader_.set_int("uUseShading", 1);
    
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
    
    occluder_mesh_->draw();
}

}



