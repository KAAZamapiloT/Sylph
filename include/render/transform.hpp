#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace render
{

class Transform
{
public:
    Transform();

    
    void set_position(const Eigen::Vector3f& position);
    const Eigen::Vector3f& position() const;

    void set_rotation(const Eigen::Quaternionf& rotation);
    const Eigen::Quaternionf& rotation() const;

   
    void set_scale(const Eigen::Vector3f& scale);
    const Eigen::Vector3f& scale() const;

   
    const Eigen::Matrix4f& model_matrix() const;

private:
    void update_model_matrix() const;

private:
    Eigen::Vector3f position_;
    Eigen::Quaternionf rotation_;
    Eigen::Vector3f scale_;

    mutable Eigen::Matrix4f model_matrix_;
    mutable bool dirty_;
};

} 