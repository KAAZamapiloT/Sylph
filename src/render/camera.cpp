#include "render/camera.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace render
{

namespace
{
constexpr float kPi = 3.14159265358979323846f;
constexpr float kHalfPi = kPi * 0.5f;
constexpr float kPitchLimit = kHalfPi - 0.01f;
}

Camera::Camera()
    : position_(0.0f, 0.0f, 3.0f),
      forward_(0.0f, 0.0f, -1.0f),
      right_(1.0f, 0.0f, 0.0f),
      up_(0.0f, 1.0f, 0.0f),
      world_up_(0.0f, 1.0f, 0.0f),
      yaw_(-kHalfPi),
      pitch_(0.0f),
      fov_degrees_(60.0f),
      aspect_ratio_(16.0f / 9.0f),
      near_plane_(0.05f),
      far_plane_(500.0f)
{
    rebuild_basis();
}

const Eigen::Vector3f& Camera::position() const noexcept
{
    return position_;
}

const Eigen::Vector3f& Camera::forward() const noexcept
{
    return forward_;
}

const Eigen::Vector3f& Camera::right() const noexcept
{
    return right_;
}

const Eigen::Vector3f& Camera::up() const noexcept
{
    return up_;
}

void Camera::set_position(const Eigen::Vector3f& position)
{
    position_ = position;
}

void Camera::set_aspect_ratio(float aspect_ratio)
{
    if (!(aspect_ratio > 0.0f))
        throw std::invalid_argument("Camera aspect ratio must be positive.");

    aspect_ratio_ = aspect_ratio;
}

void Camera::translate(const Eigen::Vector3f& delta)
{
    position_ += delta;
}

void Camera::move_local(const Eigen::Vector3f& local_delta)
{
    position_ +=
        right_ * local_delta.x() +
        up_ * local_delta.y() +
        forward_ * local_delta.z();
}

void Camera::rotate(float pitch_delta, float yaw_delta, float /*roll_delta*/)
{
    pitch_ = std::clamp(
        pitch_ + pitch_delta,
        -kPitchLimit,
        kPitchLimit);

    yaw_ += yaw_delta;

    // Keep yaw numerically bounded without changing its orientation.
    if (yaw_ > kPi || yaw_ < -kPi)
        yaw_ = std::remainder(yaw_, 2.0f * kPi);

    rebuild_basis();
}

void Camera::Translate(const Eigen::Vector3f& delta)
{
    translate(delta);
}

void Camera::Rotate(float pitch_delta, float yaw_delta, float roll_delta)
{
    rotate(pitch_delta, yaw_delta, roll_delta);
}

void Camera::reset()
{
    position_ = Eigen::Vector3f(0.0f, 0.0f, 3.0f);
    yaw_ = -kHalfPi;
    pitch_ = 0.0f;
    rebuild_basis();
}

void Camera::set_fov_degrees(float fov_degrees)
{
    if (!(fov_degrees > 1.0f && fov_degrees < 179.0f))
        throw std::invalid_argument("Camera FOV must be in (1, 179) degrees.");

    fov_degrees_ = fov_degrees;
}

void Camera::set_clip_planes(float near_plane, float far_plane)
{
    if (!(near_plane > 0.0f) || !(far_plane > near_plane))
        throw std::invalid_argument("Camera clip planes are invalid.");

    near_plane_ = near_plane;
    far_plane_ = far_plane;
}

float Camera::fov_degrees() const noexcept
{
    return fov_degrees_;
}

float Camera::aspect_ratio() const noexcept
{
    return aspect_ratio_;
}

Eigen::Matrix4f Camera::view_matrix() const
{
    // Right-handed OpenGL look-at matrix.
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    view(0, 0) = right_.x();
    view(0, 1) = right_.y();
    view(0, 2) = right_.z();
    view(0, 3) = -right_.dot(position_);

    view(1, 0) = up_.x();
    view(1, 1) = up_.y();
    view(1, 2) = up_.z();
    view(1, 3) = -up_.dot(position_);

    view(2, 0) = -forward_.x();
    view(2, 1) = -forward_.y();
    view(2, 2) = -forward_.z();
    view(2, 3) = forward_.dot(position_);

    return view;
}

Eigen::Matrix4f Camera::projection_matrix() const
{
    const float radians =
        fov_degrees_ * (kPi / 180.0f);

    const float tan_half_fov =
        std::tan(radians * 0.5f);

    const float focal = 1.0f / tan_half_fov;

    Eigen::Matrix4f projection = Eigen::Matrix4f::Zero();

    projection(0, 0) = focal / aspect_ratio_;
    projection(1, 1) = focal;
    projection(2, 2) =
        (far_plane_ + near_plane_) /
        (near_plane_ - far_plane_);
    projection(2, 3) =
        (2.0f * far_plane_ * near_plane_) /
        (near_plane_ - far_plane_);
    projection(3, 2) = -1.0f;

    return projection;
}

void Camera::rebuild_basis()
{
    const float cos_pitch = std::cos(pitch_);

    forward_ = Eigen::Vector3f(
        cos_pitch * std::cos(yaw_),
        std::sin(pitch_),
        cos_pitch * std::sin(yaw_));
    forward_.normalize();

    right_ = forward_.cross(world_up_);
    right_.normalize();

    up_ = right_.cross(forward_);
    up_.normalize();
}

} // namespace render
