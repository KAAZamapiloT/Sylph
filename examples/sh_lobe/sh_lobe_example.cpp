#include "sh_lobe_example.hpp"

#include <iostream>

namespace examples
{

namespace
{

constexpr const char* VERTEX_SHADER = R"(
#version 330 core

layout(location = 0) in vec3 a_position;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec2 a_uv;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;

out vec3 v_normal;

void main()
{
    v_normal = mat3(uModel) * a_normal;

    gl_Position =
        uProjection *
        uView *
        uModel *
        vec4(a_position, 1.0);
}
)";

constexpr const char* FRAGMENT_SHADER = R"(
#version 330 core

in vec3 v_normal;

out vec4 frag_color;

void main()
{
    vec3 n = normalize(v_normal);

    vec3 light_dir =
        normalize(vec3(0.5, 0.8, 1.0));

    float diffuse =
        max(dot(n, light_dir), 0.0);

    vec3 base =
        vec3(0.15, 0.45, 0.95);

    vec3 color =
        base * (0.15 + 0.85 * diffuse);

    frag_color =
        vec4(color, 1.0);
}
)";

} // namespace


SHLobeExample::SHLobeExample(
    int width,
    int height
)
    /*
     * Order 4:
     *
     * l = 0
     * l = 1
     * l = 2
     * l = 3
     */
    :
    coefficients_(4),
    coeff_a_(4),
    coeff_b_(4),

    lobe_(
        coefficients_
    ),

    shader_(
        VERTEX_SHADER,
        FRAGMENT_SHADER
    )
{
    // --------------------------------------------------------
    // Test function: Y_1^0
    // --------------------------------------------------------

    coefficients_(1, 0) = 1.0;

    /*
     * IMPORTANT:
     *
     * lobe_ was constructed before the coefficient was changed.
     *
     * Rebuild after modifying coefficients.
     */
    lobe_.rebuild(coefficients_);

    // --------------------------------------------------------
    // Camera
    // --------------------------------------------------------

    camera_.set_position(
        Eigen::Vector3f(
            0.0f,
            0.0f,
            3.5f
        )
    );

    camera_.set_fov_degrees(45.0f);

    camera_.set_clip_planes(
        0.1f,
        100.0f
    );

    resize(width, height);

    // --------------------------------------------------------
    // Camera controls
    // --------------------------------------------------------

    camera_controller_.set_move_speed(4.0f);
    camera_controller_.set_sprint_multiplier(3.0f);
    camera_controller_.set_slow_multiplier(0.30f);
    camera_controller_.set_mouse_sensitivity(0.0025f);

    // --------------------------------------------------------
    // Object
    // --------------------------------------------------------

    transform_.set_position(
        Eigen::Vector3f::Zero()
    );

    transform_.set_scale(
        Eigen::Vector3f::Ones()
    );
}


void SHLobeExample::process_event(
    const SDL_Event& event,
    SDL_Window* window
)
{
    camera_controller_.process_event(
        event,
        window
    );

    camera_controller_.process_mouse_motion(
        camera_,
        event
    );

    if (event.type == SDL_EVENT_KEY_DOWN) {
        if (event.key.scancode == SDL_SCANCODE_1) current_view_ = 1;
        if (event.key.scancode == SDL_SCANCODE_2) current_view_ = 2;
        if (event.key.scancode == SDL_SCANCODE_3) current_view_ = 3;
    }
}


void SHLobeExample::update(float dt)
{
    camera_controller_.apply_pending_actions(
        camera_
    );

    camera_controller_.update(
        camera_,
        dt
    );

    // 1. Animate Shape A over time (Moving Directional Lobe)
    time_ += dt;
    coeff_a_(0, 0) = 0.8f; // Base spherical shape
    coeff_a_(1, 0) = std::cos(time_ * 2.0f); // Y-axis movement
    coeff_a_(1, 1) = std::sin(time_ * 2.0f); // X-axis movement

    // 2. Shape B is static (Static Quadrupole Lobe)
    coeff_b_(0, 0) = 0.5f; // Base spherical shape
    coeff_b_(2, 0) = 1.0f; // Pinching effect

    // 2. THE PAPER'S ALGORITHM: Multiply A and B
    sylph::SHCoefficients result = sh_transform_.product(coeff_a_, coeff_b_, 4);

    // 3. Rebuild the mesh based on what the user wants to see
    if (current_view_ == 1) {
        lobe_.rebuild(coeff_a_);
    } else if (current_view_ == 2) {
        lobe_.rebuild(coeff_b_);
    } else {
        lobe_.rebuild(result); // Watch the product morph dynamically!
    }
}


void SHLobeExample::resize(
    int width,
    int height
)
{
    if (height <= 0)
        return;

    camera_.set_aspect_ratio(
        static_cast<float>(width) /
        static_cast<float>(height)
    );
}


void SHLobeExample::render(render::Renderer& renderer)
{
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    renderer.draw(
        lobe_.mesh(),
        transform_,
        shader_,
        camera_
    );
}

} // namespace examples
