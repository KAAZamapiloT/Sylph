#include "glossy_example.hpp"
#include "render/obj_loader.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/gauss_legendre.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <glad/gl.h>

namespace examples
{

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

GlossyExample::GlossyExample(int width, int height)
    : shader_(load_shader_source("shaders/glossy.vert"), load_shader_source("shaders/glossy.frag")),
      lighting_(8)
{
    std::cout << "Loading mesh...\n";
    render::MeshData mesh_data = render::ObjLoader::load("assets/meshes/bunny.obj");

    std::vector<std::uint32_t> indices = mesh_data.indices;
    std::vector<render::ColorVertex> vertices;
    vertices.reserve(mesh_data.vertices.size());
    
    Eigen::Vector3f min_b = mesh_data.vertices[0].position;
    Eigen::Vector3f max_b = mesh_data.vertices[0].position;

    for (const auto& v : mesh_data.vertices) {
        min_b = min_b.cwiseMin(v.position);
        max_b = max_b.cwiseMax(v.position);
        
        render::ColorVertex cv;
        cv.position = v.position;
        cv.color = v.normal; 
        vertices.push_back(cv);
    }
    
    Eigen::Vector3f center = (min_b + max_b) * 0.5f;
    float extent = (max_b - min_b).maxCoeff();
    transform_.set_position(-center * (2.0f / extent));
    transform_.set_scale(Eigen::Vector3f::Constant(2.0f / extent));

    mesh_ = std::make_unique<render::ColorMesh>(vertices, indices);

    camera_.set_position(Eigen::Vector3f(0.0f, 0.0f, 3.0f));
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);

    camera_controller_.set_move_speed(4.0f);
    camera_controller_.set_sprint_multiplier(3.0f);

    update_roughness();
    
    std::cout << "-------------------------------------------\n";
    std::cout << "TRUE Paper Algorithm (Spherical Grid Triple Product)\n";
    std::cout << "Controls:\n";
    std::cout << " UP/DOWN: Adjust Roughness\n";
    std::cout << " 1-8: Adjust approximation parameter 't' (truncation)\n";
    std::cout << "-------------------------------------------\n";
}

void GlossyExample::update_grid_data()
{
    int nTheta = (t_param_ + 1) / 2;
    int nPhi = t_param_;
    int grid_size = nTheta * nPhi;
    
    grid_light_.resize(grid_size);
    grid_weights_.resize(grid_size);
    grid_dirs_.resize(grid_size);
    grid_vis_.resize(grid_size);
    
    sylph::gauss_legendre gl(nTheta);
    const auto& nodes = gl.nodes();
    const auto& weights = gl.weights();
    
    sylph::SphericalHarmonics sh;
    
    double pi = std::numbers::pi_v<double>;
    double dPhi = 2.0 * pi / nPhi;
    
    for (int i = 0; i < nTheta; ++i) {
        double theta = std::acos(nodes[i]);
        double w_theta = weights[i];
        
        for (int j = 0; j < nPhi; ++j) {
            double phi = j * dPhi;
            int idx = i * nPhi + j;
            
            // Grid weight = w_theta * (2*pi / N_phi)
            grid_weights_[idx] = static_cast<float>(w_theta * dPhi);
            
            // Direction vector
            float x = static_cast<float>(std::sin(theta) * std::cos(phi));
            float z = static_cast<float>(std::sin(theta) * std::sin(phi));
            float y = static_cast<float>(std::cos(theta));
            grid_dirs_[idx] = Eigen::Vector3f(x, y, z);
            
            // Reconstruct Light at (theta, phi)
            Eigen::Vector3f L(0.0f, 0.0f, 0.0f);
            for(int c=0; c<3; ++c) {
                const auto& coeffs = lighting_.coefficients(static_cast<render::LightChannel>(c));
                double val = 0.0;
                for (int l = 0; l < lighting_.order(); ++l) {
                    for (int m = -l; m <= l; ++m) {
                        val += coeffs(l, m) * sh.evaluate(l, m, theta, phi);
                    }
                }
                L[c] = static_cast<float>(val);
            }
            grid_light_[idx] = L;
            
            // Procedural Visibility (e.g. shadow from below)
            grid_vis_[idx] = (y < -0.2f) ? 0.2f : 1.0f; // Fake ambient occlusion
        }
    }
}

void GlossyExample::update_roughness()
{
    sylph::SphericalHarmonics sh;
    float s = roughness_;
    
    base_brdf_ = sh.project(
        [s](double theta, double phi) {
            (void)phi;
            if (theta >= std::numbers::pi_v<double> / 2.0) return 0.0;
            return std::pow(std::cos(theta), s);
        },
        8, 128, 64
    );
    
    update_grid_data(); // Rebuild grid if t changes (handled in process_event)
}

void GlossyExample::process_event(const SDL_Event& event, SDL_Window* window)
{
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
    
    if (event.type == SDL_EVENT_KEY_DOWN) {
        if (event.key.scancode == SDL_SCANCODE_UP) {
            roughness_ = std::min(100.0f, roughness_ + 5.0f);
            update_roughness();
            std::cout << "Roughness: " << roughness_ << "\n";
        }
        if (event.key.scancode == SDL_SCANCODE_DOWN) {
            roughness_ = std::max(1.0f, roughness_ - 5.0f);
            update_roughness();
            std::cout << "Roughness: " << roughness_ << "\n";
        }
        
        int new_t = t_param_;
        if (event.key.scancode == SDL_SCANCODE_1) new_t = 1;
        if (event.key.scancode == SDL_SCANCODE_2) new_t = 2;
        if (event.key.scancode == SDL_SCANCODE_3) new_t = 3;
        if (event.key.scancode == SDL_SCANCODE_4) new_t = 4;
        if (event.key.scancode == SDL_SCANCODE_5) new_t = 5;
        if (event.key.scancode == SDL_SCANCODE_6) new_t = 6;
        if (event.key.scancode == SDL_SCANCODE_7) new_t = 7;
        if (event.key.scancode == SDL_SCANCODE_8) new_t = 8;
        
        if (new_t != t_param_) {
            t_param_ = new_t;
            std::cout << "t = " << t_param_ << "\n";
            update_grid_data();
        }
    }
}

void GlossyExample::generate_lighting(float time)
{
    lighting_.clear();
    
    Eigen::Vector3f light_dir(
        std::cos(time * 0.5f),
        0.5f,
        std::sin(time * 0.5f)
    );
    light_dir.normalize();
    
    double theta = std::acos(std::clamp(light_dir.y(), -1.0f, 1.0f));
    double phi = std::atan2(light_dir.z(), light_dir.x());
    
    sylph::SphericalHarmonics sh;
    for (int l = 0; l < lighting_.order(); ++l) {
        for (int m = -l; m <= l; ++m) {
            double val = sh.evaluate(l, m, theta, phi);
            
            // We NO LONGER convolve on the CPU! We pass the RAW lighting SH coefficients.
            lighting_.coefficients(render::LightChannel::Red)(l, m) = val * 5.0;
            lighting_.coefficients(render::LightChannel::Green)(l, m) = val * 4.0;
            lighting_.coefficients(render::LightChannel::Blue)(l, m) = val * 3.0;
        }
    }
    
    update_grid_data();
}

void GlossyExample::update(float dt)
{
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, dt);

    time_ += dt;
    generate_lighting(time_);
}

void GlossyExample::resize(int width, int height)
{
    if (height > 0)
        camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void GlossyExample::render(render::Renderer&)
{
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    shader_.bind();
    shader_.set_mat4("uModel", transform_.model_matrix());
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix());

    // Send Grid Data to Shader
    shader_.set_int("uTParam", t_param_);
    shader_.set_int("uGridSize", static_cast<int>(grid_light_.size()));
    
    // Upload Zonal Harmonics of BRDF up to t_param_
    for (int l = 0; l < t_param_; ++l) {
        double ZH_l = std::sqrt((4.0 * std::numbers::pi_v<double>) / (2 * l + 1)) * base_brdf_(l, 0);
        shader_.set_float("uBrdfZH[" + std::to_string(l) + "]", static_cast<float>(ZH_l));
    }
    
    // For arrays, we just use glUniform since Shader class doesn't have vector arrays
    GLint locLight = glGetUniformLocation(shader_.id(), "uLightGrid");
    if (locLight != -1) glUniform3fv(locLight, grid_light_.size(), (const GLfloat*)grid_light_.data());
    
    GLint locVis = glGetUniformLocation(shader_.id(), "uVisGrid");
    if (locVis != -1) glUniform1fv(locVis, grid_vis_.size(), grid_vis_.data());
    
    GLint locWeights = glGetUniformLocation(shader_.id(), "uGridWeights");
    if (locWeights != -1) glUniform1fv(locWeights, grid_weights_.size(), grid_weights_.data());
    
    GLint locDirs = glGetUniformLocation(shader_.id(), "uGridDirs");
    if (locDirs != -1) glUniform3fv(locDirs, grid_dirs_.size(), (const GLfloat*)grid_dirs_.data());

    mesh_->draw();
}

} // namespace examples

