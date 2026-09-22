#include "house_example.hpp"
#include "render/obj_loader.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/gauss_legendre.hpp"
#include "sylph/spherical_grid.hpp"
#include "sylph/sh_grid_transform.hpp"
#include "sylph/spherical_grid.hpp"
#include "sylph/sh_grid_transform.hpp"
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

HouseExample::HouseExample(int width, int height)
    : shader_(load_shader_source("shaders/glossy.vert"), load_shader_source("shaders/glossy.frag")),
      lighting_(8)
{
    std::cout << "Loading full-scale house and display meshes...\n";
    
    auto load_mesh = [](const std::string& path, bool recenter) {
        render::MeshData mesh_data = render::ObjLoader::load(path);
        if (recenter) {
            Eigen::Vector3f min_b = mesh_data.vertices[0].position;
            Eigen::Vector3f max_b = mesh_data.vertices[0].position;
            for (const auto& v : mesh_data.vertices) {
                min_b = min_b.cwiseMin(v.position);
                max_b = max_b.cwiseMax(v.position);
            }
            Eigen::Vector3f center = (min_b + max_b) * 0.5f;
            float extent = (max_b - min_b).maxCoeff();
            for (auto& v : mesh_data.vertices) {
                v.position = (v.position - center) * (2.0f / extent);
            }
        }
        return std::make_shared<render::Mesh>(mesh_data.vertices, mesh_data.indices);
    };

    Eigen::Vector3f color_floor(0.4f, 0.25f, 0.15f);
    Eigen::Vector3f color_wall(0.9f, 0.9f, 0.85f);
    Eigen::Vector3f color_gold(1.0f, 0.8f, 0.2f);
    Eigen::Vector3f color_jade(0.3f, 0.8f, 0.5f);
    Eigen::Vector3f color_clay(0.7f, 0.6f, 0.5f);

    auto mesh_floor = load_mesh("assets/house/Floor_Square.obj", false);
    auto mesh_wall = load_mesh("assets/house/Wall_Straight_Main.obj", false);
    auto mesh_window = load_mesh("assets/house/Exterior_Window_Double.obj", false);
    auto mesh_doorway = load_mesh("assets/house/Wall_Doorway.obj", false);
    
    auto mesh_bunny = load_mesh("assets/meshes/bunny.obj", true);
    auto mesh_teapot = load_mesh("assets/meshes/teapot.obj", true);
    auto mesh_dragon = load_mesh("assets/meshes/dragon.obj", true);
    
    for (int x = -5; x <= 4; ++x) {
        for (int z = -5; z <= 4; ++z) {
            SceneObject obj;
            obj.mesh = mesh_floor;
            obj.transform.set_position(Eigen::Vector3f(float(x)*4.0f + 2.0f, -0.5f, float(z)*4.0f + 2.0f));
            obj.transform.set_scale(Eigen::Vector3f(4.0f, 4.0f, 4.0f));
            obj.roughness = 40.0f;
            obj.base_color = color_floor;
            objects_.push_back(std::move(obj));
        }
    }
    
    for (int x = -5; x <= 4; ++x) {
        SceneObject back_wall;
        back_wall.mesh = (x == 0) ? mesh_window : mesh_wall;
        back_wall.transform.set_position(Eigen::Vector3f(float(x)*4.0f + 2.0f, -0.5f, -18.0f));
        back_wall.transform.set_scale(Eigen::Vector3f(4.0f, 4.0f, 4.0f));
        back_wall.roughness = 10.0f;
        back_wall.base_color = color_wall;
        objects_.push_back(std::move(back_wall));
        
        SceneObject front_wall;
        front_wall.mesh = (x == 0) ? mesh_doorway : mesh_wall;
        front_wall.transform.set_position(Eigen::Vector3f(float(x)*4.0f + 2.0f, -0.5f, 22.0f));
        front_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(std::numbers::pi_v<float>, Eigen::Vector3f::UnitY())));
        front_wall.transform.set_scale(Eigen::Vector3f(4.0f, 4.0f, 4.0f));
        front_wall.roughness = 10.0f;
        front_wall.base_color = color_wall;
        objects_.push_back(std::move(front_wall));
    }
    for (int z = -4; z <= 4; ++z) {
        SceneObject left_wall;
        left_wall.mesh = mesh_wall;
        left_wall.transform.set_position(Eigen::Vector3f(-18.0f, -0.5f, float(z)*4.0f + 2.0f));
        left_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(-std::numbers::pi_v<float>/2.0f, Eigen::Vector3f::UnitY())));
        left_wall.transform.set_scale(Eigen::Vector3f(4.0f, 4.0f, 4.0f));
        left_wall.roughness = 10.0f;
        left_wall.base_color = color_wall;
        objects_.push_back(std::move(left_wall));
        
        SceneObject right_wall;
        right_wall.mesh = mesh_wall;
        right_wall.transform.set_position(Eigen::Vector3f(22.0f, -0.5f, float(z)*4.0f + 2.0f));
        right_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(std::numbers::pi_v<float>/2.0f, Eigen::Vector3f::UnitY())));
        right_wall.transform.set_scale(Eigen::Vector3f(4.0f, 4.0f, 4.0f));
        right_wall.roughness = 10.0f;
        right_wall.base_color = color_wall;
        objects_.push_back(std::move(right_wall));
    }
    
    SceneObject bunny;
    bunny.mesh = mesh_bunny;
    bunny.transform.set_position(Eigen::Vector3f(-4.0f, 0.5f, 0.0f));
    bunny.transform.set_scale(Eigen::Vector3f(2.0f, 2.0f, 2.0f));
    bunny.roughness = 10.0f;
    bunny.base_color = color_clay;
    objects_.push_back(std::move(bunny));
    
    SceneObject teapot;
    teapot.mesh = mesh_teapot;
    teapot.transform.set_position(Eigen::Vector3f(0.0f, 0.5f, 0.0f));
    teapot.transform.set_scale(Eigen::Vector3f(2.0f, 2.0f, 2.0f));
    teapot.roughness = 40.0f;
    teapot.base_color = color_gold;
    objects_.push_back(std::move(teapot));
    
    SceneObject dragon;
    dragon.mesh = mesh_dragon;
    dragon.transform.set_position(Eigen::Vector3f(4.0f, 0.5f, 0.0f));
    dragon.transform.set_scale(Eigen::Vector3f(2.0f, 2.0f, 2.0f));
    dragon.roughness = 150.0f;
    dragon.base_color = color_jade;
    objects_.push_back(std::move(dragon));

    camera_.set_position(Eigen::Vector3f(0.0f, 1.5f, 8.0f));
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);
    
    camera_controller_.set_move_speed(8.0f);
    camera_controller_.set_sprint_multiplier(3.0f);

    auto& light_r = lighting_.coefficients(render::LightChannel::Red);
    auto& light_g = lighting_.coefficients(render::LightChannel::Green);
    auto& light_b = lighting_.coefficients(render::LightChannel::Blue);
    light_r(0, 0) = 0.5; light_g(0, 0) = 0.5; light_b(0, 0) = 0.6;
    light_r(1, 0) = 1.0; light_g(1, 0) = 1.0; light_b(1, 0) = 0.9;
    
    t_param_ = 6;
    update_grid_data();
}

void HouseExample::update_grid_data() {
    int nTheta = (t_param_ + 1) / 2;
    int nPhi = t_param_;
    int n_points = nTheta * nPhi;
    
    grid_light_.assign(n_points, Eigen::Vector3f::Zero());
    grid_vis_.assign(n_points, 1.0f);
    grid_weights_.assign(n_points, 0.0f);
    grid_dirs_.assign(n_points, Eigen::Vector3f::Zero());
    
    sylph::gauss_legendre gl(nTheta);
    double dPhi = 2.0 * std::numbers::pi_v<double> / nPhi;
    
    sylph::SphericalHarmonics sh;
    
    for (int i = 0; i < nTheta; ++i) {
        double theta = std::acos(gl.nodes()[i]);
        double w_theta = gl.weights()[i];
        
        for (int j = 0; j < nPhi; ++j) {
            double phi = j * dPhi;
            int idx = i * nPhi + j;
            
            grid_weights_[idx] = static_cast<float>(w_theta * dPhi);
            
            float x = static_cast<float>(std::sin(theta) * std::cos(phi));
            float z = static_cast<float>(std::sin(theta) * std::sin(phi));
            float y = static_cast<float>(std::cos(theta));
            grid_dirs_[idx] = Eigen::Vector3f(x, y, z);
            
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
        }
    }
}
void HouseExample::process_event(const SDL_Event& event, SDL_Window* window) {
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
    
    if (event.type == SDL_EVENT_KEY_DOWN && event.key.key == SDLK_T) {
        t_param_ = std::clamp(t_param_ + 1, 1, 8);
        update_grid_data();
    }
    if (event.type == SDL_EVENT_KEY_DOWN && event.key.key == SDLK_G) {
        t_param_ = std::clamp(t_param_ - 1, 1, 8);
        update_grid_data();
    }
}

void HouseExample::update(float dt) {
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, dt);
    
    objects_.back().transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(dt * 0.5f, Eigen::Vector3f::UnitY())));
}

void HouseExample::resize(int width, int height) {
    glViewport(0, 0, width, height);
    camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void HouseExample::render(render::Renderer& renderer) {
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    
    shader_.bind();
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix());
    
    shader_.set_int("uTParam", t_param_);
    shader_.set_int("uGridSize", static_cast<int>(static_cast<GLsizei>(grid_light_.size())));
    
    GLint locLight = glGetUniformLocation(shader_.id(), "uLightGrid");
    if (locLight != -1) glUniform3fv(locLight, static_cast<GLsizei>(grid_light_.size()), (const GLfloat*)grid_light_.data());
    
    GLint locVis = glGetUniformLocation(shader_.id(), "uVisGrid");
    if (locVis != -1) glUniform1fv(locVis, static_cast<GLsizei>(grid_vis_.size()), grid_vis_.data());
    
    GLint locWeights = glGetUniformLocation(shader_.id(), "uGridWeights");
    if (locWeights != -1) glUniform1fv(locWeights, static_cast<GLsizei>(grid_weights_.size()), grid_weights_.data());
    
    GLint locDirs = glGetUniformLocation(shader_.id(), "uGridDirs");
    if (locDirs != -1) glUniform3fv(locDirs, static_cast<GLsizei>(grid_dirs_.size()), (const GLfloat*)grid_dirs_.data());

    for (const auto& obj : objects_) {
        shader_.set_mat4("uModel", obj.transform.model_matrix());
        
        // Pass base color to the shader!
        shader_.set_vec3("uBaseColor", obj.base_color);
        
        for (int l = 0; l < t_param_; ++l) {
            double ZH_l = std::sqrt((4.0 * std::numbers::pi_v<double>) / (2 * l + 1)) * obj.base_brdf(l, 0);
            shader_.set_float("uBrdfZH[" + std::to_string(l) + "]", static_cast<float>(ZH_l));
        }
        
        obj.mesh->draw();
    }
}
}







