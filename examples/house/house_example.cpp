/*
 * ==============================================================================
 * EDUCATIONAL GUIDE: REAL-TIME GLOSSY PRT USING SPHERICAL GRIDS
 * ==============================================================================
 * 
 * This file demonstrates the core algorithm of the paper "Fast and Accurate 
 * Spherical Harmonics Products" applied to Glossy Precomputed Radiance Transfer.
 * 
 * THE CHALLENGE:
 * To render a glossy object with self-shadows, we need to solve the Rendering Equation 
 * for every pixel on the screen. This requires computing the integral over the sphere 
 * of three spherical functions multiplied together:
 *    Integral( L(w) * B(w) * V(w) ) dw
 * 
 * - L(w): Incident Light (Changes globally over time)
 * - B(w): BRDF Material (Changes per-pixel based on camera reflection vector R)
 * - V(w): Visibility/Shadows (Static, changes per-vertex based on geometry)
 * 
 * If we represent L, B, and V as Spherical Harmonics (SH) coefficients, multiplying 
 * them in frequency space requires Clebsch-Gordan coefficients. This is incredibly 
 * slow (O(n^5) complexity) and impossible to do per-pixel in real-time.
 * 
 * THE SPHERICAL GRID SOLUTION:
 * The paper proposes a brilliant alternative: instead of multiplying the SH 
 * coefficients directly, we transform the SH functions back into the spatial domain, 
 * but only at very specific, carefully chosen angles—a "Spherical Grid".
 * 
 * At these grid points, we simply do normal multiplication (L * B * V). Then, we 
 * multiply each point by a "quadrature weight" and sum them all up. This magically 
 * gives us the exact integral we need in a fraction of the time!
 * ==============================================================================
 */

#include "house_example.hpp"
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

HouseExample::HouseExample(int width, int height)
    : shader_(load_shader_source("shaders/glossy.vert"), load_shader_source("shaders/glossy.frag")),
      lighting_(8)
{
    std::cout << "Loading full-scale house and display meshes...\n";
    
    // 1. Helper function to load a mesh with a specific color
    auto load_mesh = [](const std::string& path, bool recenter, const Eigen::Vector3f& color) {
        render::MeshData mesh_data = render::ObjLoader::load(path);
        std::vector<render::ColorVertex> vertices;
        vertices.reserve(mesh_data.vertices.size());
        
        Eigen::Vector3f min_b = mesh_data.vertices[0].position;
        Eigen::Vector3f max_b = mesh_data.vertices[0].position;
        for (const auto& v : mesh_data.vertices) {
            min_b = min_b.cwiseMin(v.position);
            max_b = max_b.cwiseMax(v.position);
            render::ColorVertex cv;
            cv.position = v.position;
            cv.color = color; // Use the provided base color instead of the normal!
            vertices.push_back(cv);
        }
        
        if (recenter) {
            Eigen::Vector3f center = (min_b + max_b) * 0.5f;
            float extent = (max_b - min_b).maxCoeff();
            for (auto& v : vertices) {
                v.position = (v.position - center) * (2.0f / extent);
            }
        }
        
        return std::make_shared<render::ColorMesh>(vertices, mesh_data.indices);
    };

    // 2. Load the base assets with realistic colors
    Eigen::Vector3f color_floor(0.4f, 0.25f, 0.15f); // Brown wood
    Eigen::Vector3f color_wall(0.9f, 0.9f, 0.85f);   // Off-white drywall
    Eigen::Vector3f color_gold(1.0f, 0.8f, 0.2f);    // Gold Teapot
    Eigen::Vector3f color_jade(0.3f, 0.8f, 0.5f);    // Jade Dragon
    Eigen::Vector3f color_clay(0.7f, 0.6f, 0.5f);    // Clay Bunny

    auto mesh_floor = load_mesh("assets/house/Floor_Square.obj", false, color_floor);
    auto mesh_wall = load_mesh("assets/house/Wall_Straight_Main.obj", false, color_wall);
    auto mesh_window = load_mesh("assets/house/Exterior_Window_Double.obj", false, color_wall);
    auto mesh_doorway = load_mesh("assets/house/Wall_Doorway.obj", false, color_wall);
    
    auto mesh_bunny = load_mesh("assets/meshes/bunny.obj", true, color_clay);
    auto mesh_teapot = load_mesh("assets/meshes/teapot.obj", true, color_gold);
    auto mesh_dragon = load_mesh("assets/meshes/dragon.obj", true, color_jade);
    
    // 3. Build the HUGE House! (10x10 units)
    for (int x = -5; x <= 4; ++x) {
        for (int z = -5; z <= 4; ++z) {
            SceneObject obj;
            obj.mesh = mesh_floor;
            obj.transform.set_position(Eigen::Vector3f(float(x) + 0.5f, -0.5f, float(z) + 0.5f));
            obj.roughness = 40.0f; // Semi-polished wood floor
            objects_.push_back(std::move(obj));
        }
    }
    
    // Add Walls around the 10x10 perimeter
    for (int x = -5; x <= 4; ++x) {
        // Back Wall
        SceneObject back_wall;
        back_wall.mesh = (x == 0) ? mesh_window : mesh_wall; // Add a window!
        back_wall.transform.set_position(Eigen::Vector3f(float(x) + 0.5f, -0.5f, -4.5f));
        back_wall.roughness = 10.0f;
        objects_.push_back(std::move(back_wall));
        
        // Front Wall
        SceneObject front_wall;
        front_wall.mesh = (x == 0) ? mesh_doorway : mesh_wall; // Add a doorway!
        front_wall.transform.set_position(Eigen::Vector3f(float(x) + 0.5f, -0.5f, 5.5f));
        front_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(std::numbers::pi_v<float>, Eigen::Vector3f::UnitY())));
        front_wall.roughness = 10.0f;
        objects_.push_back(std::move(front_wall));
    }
    for (int z = -4; z <= 4; ++z) {
        // Left Wall
        SceneObject left_wall;
        left_wall.mesh = mesh_wall;
        left_wall.transform.set_position(Eigen::Vector3f(-4.5f, -0.5f, float(z) + 0.5f));
        left_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(-std::numbers::pi_v<float>/2.0f, Eigen::Vector3f::UnitY())));
        left_wall.roughness = 10.0f;
        objects_.push_back(std::move(left_wall));
        
        // Right Wall
        SceneObject right_wall;
        right_wall.mesh = mesh_wall;
        right_wall.transform.set_position(Eigen::Vector3f(5.5f, -0.5f, float(z) + 0.5f));
        right_wall.transform.set_rotation(Eigen::Quaternionf(Eigen::AngleAxisf(std::numbers::pi_v<float>/2.0f, Eigen::Vector3f::UnitY())));
        right_wall.roughness = 10.0f;
        objects_.push_back(std::move(right_wall));
    }
    
    // 4. Place the display objects inside the grand hall!
    SceneObject bunny;
    bunny.mesh = mesh_bunny;
    bunny.transform.set_position(Eigen::Vector3f(-2.5f, 0.5f, 0.0f));
    bunny.transform.set_scale(Eigen::Vector3f::Constant(2.0f)); // Double size!
    bunny.roughness = 10.0f; // Matte
    objects_.push_back(std::move(bunny));
    
    SceneObject teapot;
    teapot.mesh = mesh_teapot;
    teapot.transform.set_position(Eigen::Vector3f(0.0f, 0.5f, 0.0f));
    teapot.transform.set_scale(Eigen::Vector3f::Constant(2.0f));
    teapot.roughness = 40.0f; // Glossy gold
    objects_.push_back(std::move(teapot));
    
    SceneObject dragon;
    dragon.mesh = mesh_dragon;
    dragon.transform.set_position(Eigen::Vector3f(2.5f, 0.5f, 0.0f));
    dragon.transform.set_scale(Eigen::Vector3f::Constant(2.0f));
    dragon.roughness = 150.0f; // Highly polished jade
    objects_.push_back(std::move(dragon));
    

    camera_.set_position(Eigen::Vector3f(0.0f, 1.5f, 4.0f)); // Move camera up slightly for the larger scale
    camera_.set_fov_degrees(45.0f);
    camera_.set_clip_planes(0.1f, 100.0f);
    resize(width, height);

    camera_controller_.set_move_speed(4.0f);
    camera_controller_.set_sprint_multiplier(3.0f);

    update_roughness();
    
    std::cout << "-------------------------------------------\n";
    std::cout << "TRUE Paper Algorithm (Spherical Grid Triple Product)\n";
    std::cout << "Scene: 3 Objects with distinct material roughnesses!\n";
    std::cout << "Controls:\n";
    std::cout << " UP/DOWN: Adjust Global Roughness Offset\n";
    std::cout << " 1-8: Adjust approximation parameter 't' (truncation)\n";
    std::cout << "-------------------------------------------\n";
}

void HouseExample::update_grid_data()
{
    /*
     * EDUCATIONAL: HOW TO CHOOSE THE GRID RESOLUTION
     * The paper dictates that to exactly represent a product of order 't', 
     * the Spherical Grid must have size: Theta = (t+1)/2, Phi = t.
     * This makes the grid incredibly small! (e.g. t=8 means only 32 grid points).
     */
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
            
            /*
             * EDUCATIONAL: QUADRATURE WEIGHTS
             * The integral over the sphere involves summing the points multiplied by their 
             * weights. The weight for a point is exactly the Gauss-Legendre weight for theta 
             * multiplied by the Riemann sum weight for phi (2pi / N_phi).
             */
            grid_weights_[idx] = static_cast<float>(w_theta * dPhi);
            
            float x = static_cast<float>(std::sin(theta) * std::cos(phi));
            float z = static_cast<float>(std::sin(theta) * std::sin(phi));
            float y = static_cast<float>(std::cos(theta));
            grid_dirs_[idx] = Eigen::Vector3f(x, y, z);
            
            /*
             * EDUCATIONAL: CPU-SIDE FORWARD TRANSFORM
             * L(w) and V(w) don't care about the camera position. To save the GPU from 
             * evaluating 64 complex SH basis functions per grid point per pixel, we evaluate 
             * L and V into the grid values here on the CPU just once per frame!
             */
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

void HouseExample::update_roughness()
{
    sylph::SphericalHarmonics sh;
    
    // Project Zonal Harmonics for each object based on its base roughness + global modifier
    for (auto& obj : objects_) {
        float s = std::max(1.0f, obj.roughness + roughness_);
        
        obj.base_brdf = sh.project(
            [s](double theta, double phi) {
                (void)phi;
                if (theta >= std::numbers::pi_v<double> / 2.0) return 0.0;
                return std::pow(std::cos(theta), s);
            },
            8, 128, 64
        );
    }
    
    update_grid_data(); // Rebuild grid if t changes (handled in process_event)
}

void HouseExample::process_event(const SDL_Event& event, SDL_Window* window)
{
    camera_controller_.process_event(event, window);
    camera_controller_.process_mouse_motion(camera_, event);
    
    if (event.type == SDL_EVENT_KEY_DOWN) {
        if (event.key.scancode == SDL_SCANCODE_UP) {
            roughness_ += 5.0f;
            update_roughness();
            std::cout << "Global Roughness Offset: " << roughness_ << "\n";
        }
        if (event.key.scancode == SDL_SCANCODE_DOWN) {
            roughness_ -= 5.0f;
            update_roughness();
            std::cout << "Global Roughness Offset: " << roughness_ << "\n";
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

void HouseExample::generate_lighting(float time)
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

void HouseExample::update(float dt)
{
    camera_controller_.apply_pending_actions(camera_);
    camera_controller_.update(camera_, dt);

    time_ += dt;
    generate_lighting(time_);
}

void HouseExample::resize(int width, int height)
{
    if (height > 0)
        camera_.set_aspect_ratio(static_cast<float>(width) / height);
}

void HouseExample::render(render::Renderer&)
{
    shader_.bind();
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix());

    // Send Grid Data to Shader
    shader_.set_int("uTParam", t_param_);
    shader_.set_int("uGridSize", static_cast<int>(grid_light_.size()));
    
    GLint locLight = glGetUniformLocation(shader_.id(), "uLightGrid");
    if (locLight != -1) glUniform3fv(locLight, grid_light_.size(), (const GLfloat*)grid_light_.data());
    
    GLint locVis = glGetUniformLocation(shader_.id(), "uVisGrid");
    if (locVis != -1) glUniform1fv(locVis, grid_vis_.size(), grid_vis_.data());
    
    GLint locWeights = glGetUniformLocation(shader_.id(), "uGridWeights");
    if (locWeights != -1) glUniform1fv(locWeights, grid_weights_.size(), grid_weights_.data());
    
    GLint locDirs = glGetUniformLocation(shader_.id(), "uGridDirs");
    if (locDirs != -1) glUniform3fv(locDirs, grid_dirs_.size(), (const GLfloat*)grid_dirs_.data());

    // Render each object with its specific transform and Zonal Harmonics
    for (const auto& obj : objects_) {
        shader_.set_mat4("uModel", obj.transform.model_matrix());
        
        for (int l = 0; l < t_param_; ++l) {
            double ZH_l = std::sqrt((4.0 * std::numbers::pi_v<double>) / (2 * l + 1)) * obj.base_brdf(l, 0);
            shader_.set_float("uBrdfZH[" + std::to_string(l) + "]", static_cast<float>(ZH_l));
        }
        
        obj.mesh->draw();
    }
}

} // namespace examples
