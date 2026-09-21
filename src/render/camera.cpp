#include "render/camera.hpp"

#include <cmath>
#include <numbers>

namespace render
{

Camera::Camera()
    : position_(0.0f, 0.0f, 3.0f),
      orientation_(Eigen::Quaternionf::Identity()),
      fov_y_(60.0f * std::numbers::pi_v<float> / 180.0f),
      aspect_(16.0f / 9.0f),
      near_plane_(0.1f),
      far_plane_(100.0f),
      view_(Eigen::Matrix4f::Identity()),
      projection_(Eigen::Matrix4f::Identity()),
      view_projection_(Eigen::Matrix4f::Identity()),
      view_dirty_(true),
      projection_dirty_(true)
{
}



void Camera::set_position(const Eigen::Vector3f& position)
{
    position_ = position;
    view_dirty_ = true;
}

const Eigen::Vector3f& Camera::position() const
{
    return position_;
}



void Camera::set_orientation(
    const Eigen::Quaternionf& orientation)
{
    orientation_ = orientation.normalized();
    view_dirty_ = true;
}

const Eigen::Quaternionf& Camera::orientation() const
{
    return orientation_;
}



void Camera::set_fov_y(float radians)
{
    fov_y_ = radians;
    projection_dirty_ = true;
}

float Camera::fov_y() const
{
    return fov_y_;
}

void Camera::set_aspect(float aspect)
{
    aspect_ = aspect;
    projection_dirty_ = true;
}

float Camera::aspect() const
{
    return aspect_;
}

void Camera::set_near_plane(float near_plane)
{
    near_plane_ = near_plane;
    projection_dirty_ = true;
}

float Camera::near_plane() const
{
    return near_plane_;
}

void Camera::set_far_plane(float far_plane)
{
    far_plane_ = far_plane;
    projection_dirty_ = true;
}

float Camera::far_plane() const
{
    return far_plane_;
}



const Eigen::Matrix4f& Camera::view_matrix() const
{
    if (view_dirty_) {
        update_view();
    }

    return view_;
}

const Eigen::Matrix4f& Camera::projection_matrix() const
{
    if (projection_dirty_) {
        update_projection();
    }

    return projection_;
}

const Eigen::Matrix4f& Camera::view_projection_matrix() const
{
    if (view_dirty_ || projection_dirty_) {
        update_view_projection();
    }

    return view_projection_;
}

void Camera::update_view() const
{
    
    Eigen::Affine3f camera_transform =
        Eigen::Affine3f::Identity();

    camera_transform.translate(position_);
    camera_transform.rotate(orientation_);

    view_ = camera_transform.inverse().matrix();

    view_dirty_ = false;
}

void Camera::update_projection() const
{
   
 const float tan_half_fov =
        std::tan(fov_y_ * 0.5f);

    projection_.setZero();

    projection_(0, 0) =
        1.0f / (aspect_ * tan_half_fov);

    projection_(1, 1) =
        1.0f / tan_half_fov;

    projection_(2, 2) =
        -(far_plane_ + near_plane_) /
        (far_plane_ - near_plane_);

    projection_(2, 3) =
        -(2.0f * far_plane_ * near_plane_) /
        (far_plane_ - near_plane_);

    projection_(3, 2) = -1.0f;

    projection_dirty_ = false;
}

void Camera::update_view_projection() const
{
    if (view_dirty_) {
        update_view();
    }

    if (projection_dirty_) {
        update_projection();
    }

    view_projection_ =
        projection_ * view_;
}

} 