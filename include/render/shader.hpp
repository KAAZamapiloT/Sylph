/**
 * @file shader.hpp
 * @brief OpenGL Shader Program Manager.
 * 
 * Handles the compilation of Vertex and Fragment shader strings, linking them 
 * into a valid GPU program. Provides utility functions for safely uploading 
 * uniform variables (matrices, floats, vectors) to the GPU.
 */
#pragma once

#include <Eigen/Dense>

#include <cstdint>
#include <string_view>

#include <glad/gl.h>

namespace render
{

class Shader
{
public:
    Shader(
        std::string_view vertex_source,
        std::string_view fragment_source
    );

    ~Shader();

    // OpenGL program owns a GPU resource.
    Shader(const Shader&) = delete;
    Shader& operator=(const Shader&) = delete;

    Shader(Shader&& other) noexcept;
    Shader& operator=(Shader&& other) noexcept;

    // Bind this shader program.
    void bind() const;

    // Upload uniforms.
    void set_int(
        std::string_view name,
        std::int32_t value
    ) const;

    void set_float(
        std::string_view name,
        float value
    ) const;

    void set_vec3(
        std::string_view name,
        const Eigen::Vector3f& value
    ) const;

    void set_mat4(
        std::string_view name,
        const Eigen::Matrix4f& value
    ) const;

    // OpenGL program handle.
    [[nodiscard]]
    GLuint id() const;

private:
    static GLuint compile_shader(
        GLenum type,
        std::string_view source
    );

    static GLuint link_program(
        GLuint vertex_shader,
        GLuint fragment_shader
    );

    [[nodiscard]]
    GLint uniform_location(
        std::string_view name
    ) const;

private:
    GLuint program_{0};
};

} // namespace render