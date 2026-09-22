#include "oof_example.hpp"
#include "render/obj_loader.hpp"
#include "sylph/spherical_harmonics.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
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
      floor_self_visibility_(4)
{
    std::cout << "SHADOW FIELDS DEMO (Section 7.3.2)\n";

    // Setup Lighting
    auto& light_r = lighting_.coefficients(render::LightChannel::Red);
    auto& light_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_b = lighting_.coefficients(render::LightChannel::Blue);
    light_r(0, 0) = 0.5; light_g(0, 0) = 0.5; light_b(0, 0) = 0.6;
    light_r(1, 0) = 1.0; light_g(1, 0) = 1.0; light_b(1, 0) = 0.9;
    
    floor_self_visibility_.coefficients()(0, 0) = 2.0 * std::sqrt(std::numbers::pi);

    // 1. Load House Walls (Using standard Mesh with Normals!)
    auto load_wall = [](const std::string& path) {
        render::MeshData mesh_data = render::ObjLoader::load(path);
        for (auto& v : mesh_data.vertices) {
            v.position *= 4.0f; // Scale house up by 4!
        }
        return std::make_shared<render::Mesh>(mesh_data.vertices, mesh_data.indices);
    };

    Eigen::Vector3f color_wall(0.9f, 0.9f, 0.85f);
    auto mesh_wall = load_wall("assets/house/Wall_Straight_Main.obj");
    auto mesh_window = load_wall("assets/house/Exterior_Window_Double.obj");
    auto mesh_doorway = load_wall("assets/house/Wall_Doorway.obj");

    // Build walls around a 12x12 area
    for (int x = -1; x <= 1; ++x) {
        SceneObject back_wall;
        back_wall.mesh = (x == 0) ? mesh_window : mesh_wall;
        back_wall.transform.set_position(Eigen::Vector3f(x * 4.0f, 0.0f, -6.0f));
        back_wall.base_color = color_wall;
        house_objects_.push_back(std::move(back_wall));
        
        SceneObject front_wall;
        front_wall.mesh = (x == 0) ? mesh_doorway : mesh_wall;
        front_wall.transform.set_position(Eigen::Vector3f(x * 4.0f, 0.0f, 6.0f));
        front_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(std::numbers::pi, Eigen::Vector3f::UnitY())));
        front_wall.base_color = color_wall;
        house_objects_.push_back(std::move(front_wall));
    }
    for (int z = -1; z <= 1; ++z) {
        SceneObject left_wall;
        left_wall.mesh = mesh_wall;
        left_wall.transform.set_position(Eigen::Vector3f(-6.0f, 0.0f, z * 4.0f));
        left_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(-std::numbers::pi/2.0f, Eigen::Vector3f::UnitY())));
        left_wall.base_color = color_wall;
        house_objects_.push_back(std::move(left_wall));
        
        SceneObject right_wall;
        right_wall.mesh = mesh_wall;
        right_wall.transform.set_position(Eigen::Vector3f(6.0f, 0.0f, z * 4.0f));
        right_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(std::numbers::pi/2.0f, Eigen::Vector3f::UnitY())));
        right_wall.base_color = color_wall;
        house_objects_.push_back(std::move(right_wall));
    }

    // 2. Load Bunny Occluder
    render::MeshData mesh_data = render::ObjLoader::load("assets/meshes/bunny.obj");
    Eigen::Vector3f min_b = mesh_data.vertices[0].position, max_b = min_b;
    for (const auto& v : mesh_data.vertices) {
        min_b = min_b.cwiseMin(v.position); max_b = max_b.cwiseMax(v.position);
    }
    Eigen::Vector3f center = (min_b + max_b) * 0.5f;
    float extent = (max_b - min_b).maxCoeff();
    
    for (auto& v : mesh_data.vertices) {
        v.position = (v.position - center) * (1.5f / extent);
    }
    occluder_mesh_ = std::make_shared<render::Mesh>(mesh_data.vertices, mesh_data.indices);
    occluder_transform_.set_position(Eigen::Vector3f(0.0f, 1.5f, 0.0f));

    // Bake OOF (Lower quality grid to save CPU time on laptops)
    occluder_oof_ = render::OOFBaker::bake(mesh_data, 8, 4, 100);

    // 3. Build Dense Floor
    build_dense_floor();

    camera_.set_position(Eigen::Vector3f(0.0f, 2.5f, 5.0f));
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
}

void OofExample::build_dense_floor() {
    const float size = 12.0f; 
    const int segments = 50; 
    const float step = size / segments;
    
    for (int z = 0; z <= segments; ++z) {
        for (int x = 0; x <= segments; ++x) {
            render::ColorVertex v;
            v.position = Eigen::Vector3f(-size/2 + x*step, 0.0f, -size/2 + z*step);
            v.color = Eigen::Vector3f(1.0f, 1.0f, 1.0f);
            floor_vertices_.push_back(v);
        }
    }
    
    std::vector<uint32_t> indices;
    for (int z = 0; z < segments; ++z) {
        for (int x = 0; x < segments; ++x) {
            uint32_t top_left = z * (segments + 1) + x;
            uint32_t top_right = top_left + 1;
            uint32_t bottom_left = (z + 1) * (segments + 1) + x;
            uint32_t bottom_right = bottom_left + 1;
            
            indices.push_back(top_left); indices.push_back(bottom_left); indices.push_back(top_right);
            indices.push_back(top_right); indices.push_back(bottom_left); indices.push_back(bottom_right);
        }
    }
    floor_mesh_ = std::make_shared<render::ColorMesh>(floor_vertices_, indices);
}

void OofExample::process_event(const SDL_Event& event, SDL_Window* window) {
    if (event.type == SDL_EVENT_KEY_DOWN || event.type == SDL_EVENT_KEY_UP) {
        camera_controller_.process_event(event, window);
    }
    if (event.type == SDL_EVENT_MOUSE_MOTION) {
        camera_controller_.process_mouse_motion(camera_, event);
    }
}

void OofExample::update(float delta_time) {
    time_ += delta_time;
    
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, delta_time);
    
    float ox = std::sin(time_) * 1.5f;
    float oz = std::cos(time_) * 1.5f;
    float oy = 1.0f + std::sin(time_ * 2.0f) * 0.5f;
    occluder_transform_.set_position(Eigen::Vector3f(ox, oy, oz));
    
    std::vector<const render::OOF*> active_oofs = { occluder_oof_.get() };
    
    Eigen::Matrix4f floor_model = floor_transform_.model_matrix();
    Eigen::Matrix4f occluder_model = occluder_transform_.model_matrix();
    Eigen::Matrix4f occluder_inv = occluder_model.inverse();
    
    for (auto& v : floor_vertices_) {
        Eigen::Vector4f world_pos = floor_model * Eigen::Vector4f(v.position.x(), v.position.y(), v.position.z(), 1.0f);
        Eigen::Vector4f local_pos = occluder_inv * world_pos;
        
        Eigen::Vector3f shadowed_color;
        try {
            shadowed_color = shadow_field_.evaluate(
                lighting_, floor_self_visibility_, local_pos.head<3>(), active_oofs);
        } catch (const std::out_of_range&) {
            shadowed_color = shadow_field_.evaluate(
                lighting_, floor_self_visibility_, local_pos.head<3>(), {});
        }
        
        Eigen::Vector3f base_color(0.5f, 0.4f, 0.35f); 
        v.color = base_color.cwiseProduct(shadowed_color);
    }
    
    floor_mesh_->update_colors(floor_vertices_);
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
    
    // Pass SH lighting to Shader for diffuse evaluation
    float light_coeffs[27]; // 9 coeffs * 3 channels
    auto& light_r = lighting_.coefficients(render::LightChannel::Red);
    auto& light_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_b = lighting_.coefficients(render::LightChannel::Blue);
    int idx = 0;
    for (int l = 0; l <= 2; ++l) {
        for (int m = -l; m <= l; ++m) {
            light_coeffs[idx*3 + 0] = static_cast<float>(light_r(l, m));
            light_coeffs[idx*3 + 1] = static_cast<float>(light_g(l, m));
            light_coeffs[idx*3 + 2] = static_cast<float>(light_b(l, m));
            idx++;
        }
    }
    GLint loc = glGetUniformLocation(shader_.id(), "uLightSH");
    glUniform3fv(loc, 9, light_coeffs);

    // Render Ground (No SH diffuse on GPU, CPU already evaluated it)
    shader_.set_int("uIsFloor", 1);
    shader_.set_mat4("uModel", floor_transform_.model_matrix());
    floor_mesh_->draw();
    
    shader_.set_int("uIsFloor", 0);
    // Render Walls
    for(const auto& wall : house_objects_) {
        shader_.set_mat4("uModel", wall.transform.model_matrix());
        shader_.set_vec3("uBaseColor", wall.base_color);
        wall.mesh->draw();
    }
    
    // Render Occluder
    shader_.set_mat4("uModel", occluder_transform_.model_matrix());
    shader_.set_vec3("uBaseColor", Eigen::Vector3f(0.8f, 0.4f, 0.2f));
    occluder_mesh_->draw();
}

void OofExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}
}
