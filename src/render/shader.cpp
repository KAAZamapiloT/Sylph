#include "render/shader.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace render
{

Shader::Shader(
    std::string_view vertex_source,
    std::string_view fragment_source)
{
    const GLuint vertex_shader =
        compile_shader(GL_VERTEX_SHADER, vertex_source);

    GLuint fragment_shader = 0;

    try {
        fragment_shader =
            compile_shader(GL_FRAGMENT_SHADER, fragment_source);

        program_ =
            link_program(vertex_shader, fragment_shader);
    }
    catch (...) {
        if (fragment_shader != 0) {
            glDeleteShader(fragment_shader);
        }

        glDeleteShader(vertex_shader);
        throw;
    }

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
}

Shader::~Shader()
{
    if (program_ != 0) {
        glDeleteProgram(program_);
    }
}

// ------------------------------------------------------------
// Move constructor
// ------------------------------------------------------------

Shader::Shader(Shader&& other) noexcept
    : program_(other.program_)
{
    other.program_ = 0;
}

// ------------------------------------------------------------
// Move assignment
// ------------------------------------------------------------

Shader& Shader::operator=(Shader&& other) noexcept
{
    if (this == &other) {
        return *this;
    }

    if (program_ != 0) {
        glDeleteProgram(program_);
    }

    program_ = other.program_;
    other.program_ = 0;

    return *this;
}

// ------------------------------------------------------------
// Bind
// ------------------------------------------------------------

void Shader::bind() const
{
    glUseProgram(program_);
}

// ------------------------------------------------------------
// Uniform setters
// ------------------------------------------------------------

void Shader::set_int(
    std::string_view name,
    std::int32_t value) const
{
    const GLint location = uniform_location(name);

    if (location == -1) {
        return;
    }

    glUniform1i(location, value);
}

void Shader::set_float(
    std::string_view name,
    float value) const
{
    const GLint location = uniform_location(name);

    if (location == -1) {
        return;
    }

    glUniform1f(location, value);
}

void Shader::set_vec3(
    std::string_view name,
    const Eigen::Vector3f& value) const
{
    const GLint location = uniform_location(name);

    if (location == -1) {
        return;
    }

    glUniform3fv(
        location,
        1,
        value.data()
    );
}

void Shader::set_mat4(
    std::string_view name,
    const Eigen::Matrix4f& value) const
{
    const GLint location = uniform_location(name);

    if (location == -1) {
        return;
    }

    glUniformMatrix4fv(
        location,
        1,
        GL_FALSE,
        value.data()
    );
}

// ------------------------------------------------------------
// Program ID
// ------------------------------------------------------------

GLuint Shader::id() const
{
    return program_;
}

// ------------------------------------------------------------
// Shader compilation
// ------------------------------------------------------------

GLuint Shader::compile_shader(
    GLenum type,
    std::string_view source)
{
    const GLuint shader = glCreateShader(type);

    if (shader == 0) {
        throw std::runtime_error(
            "Failed to create OpenGL shader."
        );
    }

    const GLchar* source_ptr = source.data();
    const GLint source_length =
        static_cast<GLint>(source.size());

    glShaderSource(
        shader,
        1,
        &source_ptr,
        &source_length
    );

    glCompileShader(shader);

    GLint success = GL_FALSE;

    glGetShaderiv(
        shader,
        GL_COMPILE_STATUS,
        &success
    );

    if (success == GL_FALSE) {

        GLint log_length = 0;

        glGetShaderiv(
            shader,
            GL_INFO_LOG_LENGTH,
            &log_length
        );

        std::string log(
            static_cast<std::size_t>(log_length),
            '\0'
        );

        if (log_length > 0) {
            glGetShaderInfoLog(
                shader,
                log_length,
                nullptr,
                log.data()
            );
        }

        glDeleteShader(shader);

        const char* shader_type =
            type == GL_VERTEX_SHADER
                ? "vertex"
                : "fragment";

        throw std::runtime_error(
            std::string("Failed to compile ") +
            shader_type +
            " shader:\n" +
            log
        );
    }

    return shader;
}

// ------------------------------------------------------------
// Program linking
// ------------------------------------------------------------

GLuint Shader::link_program(
    GLuint vertex_shader,
    GLuint fragment_shader)
{
    const GLuint program =
        glCreateProgram();

    if (program == 0) {
        throw std::runtime_error(
            "Failed to create OpenGL shader program."
        );
    }

    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);

    glLinkProgram(program);

    GLint success = GL_FALSE;

    glGetProgramiv(
        program,
        GL_LINK_STATUS,
        &success
    );

    if (success == GL_FALSE) {

        GLint log_length = 0;

        glGetProgramiv(
            program,
            GL_INFO_LOG_LENGTH,
            &log_length
        );

        std::string log(
            static_cast<std::size_t>(log_length),
            '\0'
        );

        if (log_length > 0) {
            glGetProgramInfoLog(
                program,
                log_length,
                nullptr,
                log.data()
            );
        }

        glDeleteProgram(program);

        throw std::runtime_error(
            "Failed to link OpenGL shader program:\n" +
            log
        );
    }

    return program;
}

// ------------------------------------------------------------
// Uniform lookup
// ------------------------------------------------------------

GLint Shader::uniform_location(
    std::string_view name) const
{
    // string_view is not required to be null terminated.
    const std::string name_string(name);

    return glGetUniformLocation(
        program_,
        name_string.c_str()
    );
}

} // namespace render