#pragma once
#include<Eigen/Dense>
#include<Eigen/Geometry>

namespace render{

    class Camera{

        public:
        Camera();

        void set_position(const Eigen::Vector3f&position);

        const Eigen::Vector3f& position() const;

        void set_orientation(const Eigen::Quaternionf&orientation);
        const Eigen::Quaternionf&orientation() const;

        void set_fov_y(float radians);
        float fov_y() const;

        void set_aspect(float aspect);
        float aspect() const;

        void set_near_plane(float near_plane);
        float near_plane() const ;

        void set_far_plane(float far_plane);
        float far_plane() const;

        const Eigen::Matrix4f& view_matrix() const;
        const Eigen::Matrix4f& projection_matrix() const;
        const Eigen::Matrix4f& view_projection_matrix() const;

        private:
    void update_view() const;
    void update_projection() const;
    void update_view_projection() const;

    private:
    Eigen::Vector3f position_;
    Eigen::Quaternionf orientation_;

    float fov_y_;
    float aspect_;
    float near_plane_;
    float far_plane_;

    mutable Eigen::Matrix4f view_;
    mutable Eigen::Matrix4f projection_;
    mutable Eigen::Matrix4f view_projection_;

    mutable bool view_dirty_;
    mutable bool projection_dirty_;
    };
}