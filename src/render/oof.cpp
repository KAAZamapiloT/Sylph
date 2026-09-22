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

const Visibility& OOF::lookup(
    const Eigen::Vector3f& position
) const
{
    const Eigen::Vector3i cell =
        position_to_cell(position);

    return at(
        cell.x(),
        cell.y(),
        cell.z()
    );
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