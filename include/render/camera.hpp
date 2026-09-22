#pragma once

#include <Eigen/Dense>

namespace render
{

class Camera
{
public:
    Camera();

    // ------------------------------------------------------------
    // Position / orientation
    // ------------------------------------------------------------

    [[nodiscard]] const Eigen::Vector3f& position() const noexcept;
    [[nodiscard]] const Eigen::Vector3f& forward() const noexcept;
    [[nodiscard]] const Eigen::Vector3f& right() const noexcept;
    [[nodiscard]] const Eigen::Vector3f& up() const noexcept;

    void set_position(const Eigen::Vector3f& position);
    void set_aspect_ratio(float aspect_ratio);

    // Move in WORLD space.
    void translate(const Eigen::Vector3f& delta);

    // Move relative to the camera basis.
    // local.x = right, local.y = up, local.z = forward.
    void move_local(const Eigen::Vector3f& local_delta);

    // Angles are in radians.
    // pitch_delta -> look up/down
    // yaw_delta   -> look left/right
    // roll_delta  -> reserved for future camera roll support
    void rotate(float pitch_delta, float yaw_delta, float roll_delta = 0.0f);

    // Compatibility with the older camera example.
    void Translate(const Eigen::Vector3f& delta);
    void Rotate(float pitch_delta, float yaw_delta, float roll_delta = 0.0f);

    void reset();

    // ------------------------------------------------------------
    // Projection settings
    // ------------------------------------------------------------

    void set_fov_degrees(float fov_degrees);
    void set_clip_planes(float near_plane, float far_plane);

    [[nodiscard]] float fov_degrees() const noexcept;
    [[nodiscard]] float aspect_ratio() const noexcept;

    // ------------------------------------------------------------
    // Matrices consumed by Renderer
    // ------------------------------------------------------------

    [[nodiscard]] Eigen::Matrix4f view_matrix() const;
    [[nodiscard]] Eigen::Matrix4f projection_matrix() const;

private:
    void rebuild_basis();

private:
    Eigen::Vector3f position_;
    Eigen::Vector3f forward_;
    Eigen::Vector3f right_;
    Eigen::Vector3f up_;
    Eigen::Vector3f world_up_;

    // Camera orientation in radians.
    float yaw_;
    float pitch_;

    float fov_degrees_;
    float aspect_ratio_;
    float near_plane_;
    float far_plane_;
};

} // namespace render
