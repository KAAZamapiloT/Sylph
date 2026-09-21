#include "render/renderer.hpp"

#include "render/camera.hpp"
#include "render/mesh.hpp"
#include "render/shader.hpp"
#include "render/transform.hpp"

#include <glad/gl.h>

#include <stdexcept>

namespace render
{

void Renderer::initialize()
{
    if (initialized_) {
        return;
    }

    // --------------------------------------------------------
    // Basic OpenGL state
    // --------------------------------------------------------

    glEnable(GL_DEPTH_TEST);

    glDepthFunc(GL_LESS);

    glClearDepth(1.0);

    initialized_ = true;
}

void Renderer::resize(
    std::int32_t width,
    std::int32_t height)
{
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument(
            "Renderer dimensions must be greater than zero."
        );
    }

    width_ = width;
    height_ = height;

    glViewport(
        0,
        0,
        width_,
        height_
    );
}

void Renderer::begin_frame(
    const Eigen::Vector4f& clear_color)
{
    if (!initialized_) {
        throw std::logic_error(
            "Renderer must be initialized before begin_frame()."
        );
    }

    glClearColor(
        clear_color.x(),
        clear_color.y(),
        clear_color.z(),
        clear_color.w()
    );

    glClear(
        GL_COLOR_BUFFER_BIT |
        GL_DEPTH_BUFFER_BIT
    );
}

void Renderer::end_frame()
{
    // Nothing required here yet.
    //
    // SDL_GL_SwapWindow() belongs to the application/window
    // layer, not the renderer.
}

void Renderer::draw(
    const Mesh& mesh,
    const Transform& transform,
    const Shader& shader,
    const Camera& camera)
{
    if (!initialized_) {
        throw std::logic_error(
            "Renderer must be initialized before draw()."
        );
    }

    shader.bind();

    shader.set_mat4(
        "uModel",
        transform.model_matrix()
    );

    shader.set_mat4(
        "uView",
        camera.view_matrix()
    );

    shader.set_mat4(
        "uProjection",
        camera.projection_matrix()
    );

    mesh.draw();
}

} // namespace render