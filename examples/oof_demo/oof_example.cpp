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
    if (!file.is_open()) {
        std::cerr << "Failed to open shader: " << path << "\n";
        return "";
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

OofExample::OofExample(int width, int height)
    : shader_(load_shader_source("shaders/oof.vert"), load_shader_source("shaders/oof.frag")),
      lighting_(4),
      floor_self_visibility_(4)
{
    std::cout << "-------------------------------------------\n";
    std::cout << "SHADOW FIELDS DEMO (Section 7.3.2)\n";
    std::cout << "Multiple Product SH Integral for Soft Shadows\n";
    std::cout << "-------------------------------------------\n";

    // 1. Setup Environment Lighting (e.g. overhead directional + ambient)
    auto& light_coeffs = lighting_.coefficients(render::LightChannel::Red);
    auto& light_coeffs_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_coeffs_b = lighting_.coefficients(render::LightChannel::Blue);
    
    // Ambient term
    light_coeffs(0, 0) = 0.5;
    light_coeffs_g(0, 0) = 0.5;
    light_coeffs_b(0, 0) = 0.6;
    
    // Directional term from above (Y-axis)
    light_coeffs(1, 0) = 0.8;
    light_coeffs_g(1, 0) = 0.8;
    light_coeffs_b(1, 0) = 0.7;

    // 2. Setup Floor Self-Visibility (Unoccluded, constant 1.0)
    // Integral of 1.0 * Y_00 is 2*sqrt(pi)
    floor_self_visibility_.coefficients()(0, 0) = 2.0 * std::sqrt(std::numbers::pi);

    // 3. Load Occluder and Bake OOF!
    std::cout << "Loading occluder...\n";
    render::MeshData mesh_data = render::ObjLoader::load("assets/meshes/bunny.obj");
    
    // Center and scale the occluder for better baking
    Eigen::Vector3f min_b = mesh_data.vertices[0].position;
    Eigen::Vector3f max_b = min_b;
    for (const auto& v : mesh_data.vertices) {
        min_b = min_b.cwiseMin(v.position);
        max_b = max_b.cwiseMax(v.position);
    }
    Eigen::Vector3f center = (min_b + max_b) * 0.5f;
    float extent = (max_b - min_b).maxCoeff();
    
    std::vector<render::ColorVertex> occluder_verts;
    for (auto& v : mesh_data.vertices) {
        v.position = (v.position - center) * (2.0f / extent);
        render::ColorVertex cv;
        cv.position = v.position;
        cv.color = Eigen::Vector3f(0.8f, 0.4f, 0.2f); // Terracotta color
        occluder_verts.push_back(cv);
    }
    
    // We pass the modified mesh_data to the baker so bounds match
    mesh_data.vertices.clear();
    for(const auto& cv : occluder_verts) {
        render::Vertex v;
        v.position = cv.position;
        mesh_data.vertices.push_back(v);
    }

    occluder_mesh_ = std::make_shared<render::ColorMesh>(occluder_verts, mesh_data.indices);
    occluder_transform_.set_position(Eigen::Vector3f(0.0f, 1.0f, 0.0f));

    // BAKE THE OOF (Grid 12x12x12, Order 4, 300 rays per voxel)
    occluder_oof_ = render::OOFBaker::bake(mesh_data, 12, 4, 300);

    // 4. Build Dense Ground Plane
    build_dense_floor();

    camera_.set_position(Eigen::Vector3f(0.0f, 3.0f, 8.0f));
    // camera_.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(-0.3f, Eigen::Vector3f::UnitX()))); // handled manually or via lookat if needed
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
}

void OofExample::build_dense_floor() {
    const float size = 10.0f;
    const int segments = 80; // 80x80 grid = 6400 vertices
    const float step = size / segments;
    
    for (int z = 0; z <= segments; ++z) {
        for (int x = 0; x <= segments; ++x) {
            render::ColorVertex v;
            v.position = Eigen::Vector3f(-size/2 + x*step, 0.0f, -size/2 + z*step);
            v.color = Eigen::Vector3f(1.0f, 1.0f, 1.0f); // White, will be modulated by shadow
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
            
            indices.push_back(top_left);
            indices.push_back(bottom_left);
            indices.push_back(top_right);
            
            indices.push_back(top_right);
            indices.push_back(bottom_left);
            indices.push_back(bottom_right);
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
    
    // Animate the occluder (bobbing and moving in a circle)
    float ox = std::sin(time_) * 1.5f;
    float oz = std::cos(time_) * 1.5f;
    float oy = 1.0f + std::sin(time_ * 2.0f) * 0.5f;
    occluder_transform_.set_position(Eigen::Vector3f(ox, oy, oz));
    
    // CPU-side Shadow Field Evaluation!
    // For every vertex on the ground plane, lookup its position in the occluder's OOF
    // and evaluate the multiple product integral.
    
    std::vector<const render::OOF*> active_oofs = { occluder_oof_.get() };
    
    // Transform vertices into world space to check against OOF bounding box
    Eigen::Matrix4f floor_model = floor_transform_.model_matrix();
    Eigen::Matrix4f occluder_model = occluder_transform_.model_matrix();
    Eigen::Matrix4f occluder_inv = occluder_model.inverse();
    
    for (auto& v : floor_vertices_) {
        Eigen::Vector4f world_pos = floor_model * Eigen::Vector4f(v.position.x(), v.position.y(), v.position.z(), 1.0f);
        
        // Convert world pos into occluder's local space to sample OOF
        Eigen::Vector4f local_pos = occluder_inv * world_pos;
        
        // The evaluate function will do the product: Light * SelfVis * OOF_Vis
        Eigen::Vector3f shadowed_color;
        try {
            // OOF throws std::out_of_range if we are outside the bounding box
            shadowed_color = shadow_field_.evaluate(
                lighting_,
                floor_self_visibility_,
                local_pos.head<3>(),
                active_oofs
            );
        } catch (const std::out_of_range&) {
            // If the ground vertex is outside the occluder's OOF, it just receives unshadowed light!
            // Evaluate without the OOF:
            shadowed_color = shadow_field_.evaluate(
                lighting_,
                floor_self_visibility_,
                local_pos.head<3>(), // unused when no oofs
                {} // empty span
            );
        }
        
        // Base ground color is pale blue/gray
        Eigen::Vector3f base_color(0.8f, 0.85f, 0.9f);
        v.color = base_color.cwiseProduct(shadowed_color);
    }
    
    // Upload the modified vertices to the GPU dynamically!
    floor_mesh_->update_colors(floor_vertices_);
}

void OofExample::render() {
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    
    shader_.bind();
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix());
    
    // Render Ground
    shader_.set_mat4("uModel", floor_transform_.model_matrix());
    floor_mesh_->draw();
    
    // Render Occluder
    shader_.set_mat4("uModel", occluder_transform_.model_matrix());
    occluder_mesh_->draw();
}

void OofExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

} // namespace examples
