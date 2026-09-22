#include "render/oof.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace render {

OOF::OOF(
    int resolution_x,
    int resolution_y,
    int resolution_z,
    const Eigen::Vector3f& min,
    const Eigen::Vector3f& max,
    int visibility_order
)
    : resolution_x_(resolution_x)
    , resolution_y_(resolution_y)
    , resolution_z_(resolution_z)
    , min_(min)
    , max_(max)
    , visibility_order_(visibility_order)
{
    if (resolution_x <= 0 ||
        resolution_y <= 0 ||
        resolution_z <= 0)
    {
        throw std::invalid_argument(
            "OOF resolutions must be positive"
        );
    }

    if (visibility_order <= 0) {
        throw std::invalid_argument(
            "OOF visibility order must be positive"
        );
    }

    if ((max.array() <= min.array()).any()) {
        throw std::invalid_argument(
            "OOF max bounds must be greater than min bounds"
        );
    }

    const std::size_t cell_count =
        static_cast<std::size_t>(resolution_x_) *
        static_cast<std::size_t>(resolution_y_) *
        static_cast<std::size_t>(resolution_z_);

    cells_.reserve(cell_count);

    for (std::size_t i = 0; i < cell_count; ++i) {
        cells_.emplace_back(
            visibility_order_
        );
    }
}

int OOF::resolution_x() const noexcept
{
    return resolution_x_;
}

int OOF::resolution_y() const noexcept
{
    return resolution_y_;
}

int OOF::resolution_z() const noexcept
{
    return resolution_z_;
}

int OOF::visibility_order() const noexcept
{
    return visibility_order_;
}

const Eigen::Vector3f& OOF::min() const noexcept
{
    return min_;
}

const Eigen::Vector3f& OOF::max() const noexcept
{
    return max_;
}

Visibility& OOF::at(
    int x,
    int y,
    int z
)
{
    if (!valid_index(x, y, z)) {
        throw std::out_of_range(
            "OOF cell index out of range"
        );
    }

    return cells_[index(x, y, z)];
}

const Visibility& OOF::at(
    int x,
    int y,
    int z
) const
{
    if (!valid_index(x, y, z)) {
        throw std::out_of_range(
            "OOF cell index out of range"
        );
    }

    return cells_[index(x, y, z)];
}

const Visibility OOF::lookup(
    const Eigen::Vector3f& position
) const
{
    Eigen::Vector3f normalized = (position - min_).cwiseQuotient(max_ - min_);
    normalized = normalized.cwiseMax(0.0f).cwiseMin(1.0f);
    
    Eigen::Vector3f scaled = normalized.cwiseProduct(
        Eigen::Vector3f(
            static_cast<float>(resolution_x_ - 1),
            static_cast<float>(resolution_y_ - 1),
            static_cast<float>(resolution_z_ - 1)
        )
    );
    
    int x0 = static_cast<int>(std::floor(scaled.x()));
    int y0 = static_cast<int>(std::floor(scaled.y()));
    int z0 = static_cast<int>(std::floor(scaled.z()));
    
    int x1 = std::min(x0 + 1, resolution_x_ - 1);
    int y1 = std::min(y0 + 1, resolution_y_ - 1);
    int z1 = std::min(z0 + 1, resolution_z_ - 1);
    
    float tx = scaled.x() - x0;
    float ty = scaled.y() - y0;
    float tz = scaled.z() - z0;
    
    const Visibility& v000 = at(x0, y0, z0);
    const Visibility& v100 = at(x1, y0, z0);
    const Visibility& v010 = at(x0, y1, z0);
    const Visibility& v110 = at(x1, y1, z0);
    const Visibility& v001 = at(x0, y0, z1);
    const Visibility& v101 = at(x1, y0, z1);
    const Visibility& v011 = at(x0, y1, z1);
    const Visibility& v111 = at(x1, y1, z1);
    
    Visibility result(visibility_order_);
    
    // Trilinear interpolation of SH coefficients
    for (int l = 0; l < visibility_order_; ++l) {
        for (int m = -l; m <= l; ++m) {
            float c00 = v000.coefficients()(l, m) * (1 - tx) + v100.coefficients()(l, m) * tx;
            float c10 = v010.coefficients()(l, m) * (1 - tx) + v110.coefficients()(l, m) * tx;
            float c01 = v001.coefficients()(l, m) * (1 - tx) + v101.coefficients()(l, m) * tx;
            float c11 = v011.coefficients()(l, m) * (1 - tx) + v111.coefficients()(l, m) * tx;
            
            float c0 = c00 * (1 - ty) + c10 * ty;
            float c1 = c01 * (1 - ty) + c11 * ty;
            
            result.coefficients()(l, m) = c0 * (1 - tz) + c1 * tz;
        }
    }
    
    return result;
}

std::size_t OOF::index(
    int x,
    int y,
    int z
) const noexcept
{
    return
        static_cast<std::size_t>(z) *
        static_cast<std::size_t>(resolution_y_) *
        static_cast<std::size_t>(resolution_x_)
        +
        static_cast<std::size_t>(y) *
        static_cast<std::size_t>(resolution_x_)
        +
        static_cast<std::size_t>(x);
}

bool OOF::valid_index(
    int x,
    int y,
    int z
) const noexcept
{
    return
        x >= 0 &&
        x < resolution_x_ &&
        y >= 0 &&
        y < resolution_y_ &&
        z >= 0 &&
        z < resolution_z_;
}

Eigen::Vector3i OOF::position_to_cell(
    const Eigen::Vector3f& position
) const
{
    const Eigen::Vector3f normalized =
        (position - min_).cwiseQuotient(
            max_ - min_
        );

    const Eigen::Vector3f clamped =
        normalized.cwiseMax(
            Eigen::Vector3f::Zero()
        ).cwiseMin(
            Eigen::Vector3f::Ones()
        );

    const Eigen::Vector3f scaled =
        clamped.cwiseProduct(
            Eigen::Vector3f(
                static_cast<float>(resolution_x_),
                static_cast<float>(resolution_y_),
                static_cast<float>(resolution_z_)
            )
        );

    Eigen::Vector3i cell;

    cell.x() = std::min(
        static_cast<int>(std::floor(scaled.x())),
        resolution_x_ - 1
    );

    cell.y() = std::min(
        static_cast<int>(std::floor(scaled.y())),
        resolution_y_ - 1
    );

    cell.z() = std::min(
        static_cast<int>(std::floor(scaled.z())),
        resolution_z_ - 1
    );

    return cell;
}

} // namespace render