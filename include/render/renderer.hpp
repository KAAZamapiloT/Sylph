#pragma once

#include <Eigen/Dense>

#include <cstdint>

namespace render
{

class Camera;
class Mesh;
class Shader;
class Transform;

class Renderer
{
public:
    Renderer() = default;
    ~Renderer() = default;

    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;

    // --------------------------------------------------------
    // Renderer state
    // --------------------------------------------------------

    void initialize();

    void resize(
        std::int32_t width,
        std::int32_t height
    );

    // --------------------------------------------------------
    // Frame
    // --------------------------------------------------------

    void begin_frame(
        const Eigen::Vector4f& clear_color
    );

    void end_frame();

    // --------------------------------------------------------
    // Drawing
    // --------------------------------------------------------

    void draw(
        const Mesh& mesh,
        const Transform& transform,
        const Shader& shader,
        const Camera& camera
    );

private:
    std::int32_t width_{0};
    std::int32_t height_{0};

    bool initialized_{false};
};

} // namespace render