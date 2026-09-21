#include "render/transform.hpp"

namespace render
{

Transform::Transform()
    : position_(0.0f, 0.0f, 0.0f),
      rotation_(Eigen::Quaternionf::Identity()),
      scale_(1.0f, 1.0f, 1.0f),
      model_matrix_(Eigen::Matrix4f::Identity()),
      dirty_(true)
{
}

void Transform::set_position(
    const Eigen::Vector3f& position)
{
    position_ = position;
    dirty_ = true;
}

const Eigen::Vector3f& Transform::position() const
{
    return position_;
}


void Transform::set_rotation(
    const Eigen::Quaternionf& rotation)
{
    rotation_ = rotation.normalized();
    dirty_ = true;
}

const Eigen::Quaternionf& Transform::rotation() const
{
    return rotation_;
}



void Transform::set_scale(
    const Eigen::Vector3f& scale)
{
    scale_ = scale;
    dirty_ = true;
}

const Eigen::Vector3f& Transform::scale() const
{
    return scale_;
}


const Eigen::Matrix4f& Transform::model_matrix() const
{
    if (dirty_) {
        update_model_matrix();
    }

    return model_matrix_;
}



void Transform::update_model_matrix() const
{
    Eigen::Affine3f transform =
        Eigen::Affine3f::Identity();

    transform.translate(position_);
    transform.rotate(rotation_);
    transform.scale(scale_);

    model_matrix_ = transform.matrix();

    dirty_ = false;
}

} 