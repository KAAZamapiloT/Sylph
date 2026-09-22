#pragma once

#include "render/visibility.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace render {

class OOF {
public:
    OOF(
        int resolution_x,
        int resolution_y,
        int resolution_z,
        const Eigen::Vector3f& min,
        const Eigen::Vector3f& max,
        int visibility_order
    );

    [[nodiscard]]
    int resolution_x() const noexcept;

    [[nodiscard]]
    int resolution_y() const noexcept;

    [[nodiscard]]
    int resolution_z() const noexcept;

    [[nodiscard]]
    int visibility_order() const noexcept;

    [[nodiscard]]
    const Eigen::Vector3f& min() const noexcept;

    [[nodiscard]]
    const Eigen::Vector3f& max() const noexcept;

    [[nodiscard]]
    Visibility& at(
        int x,
        int y,
        int z
    );

    [[nodiscard]]
    const Visibility& at(
        int x,
        int y,
        int z
    ) const;

    [[nodiscard]]
    const Visibility& lookup(
        const Eigen::Vector3f& position
    ) const;

private:
    [[nodiscard]]
    std::size_t index(
        int x,
        int y,
        int z
    ) const noexcept;

    [[nodiscard]]
    bool valid_index(
        int x,
        int y,
        int z
    ) const noexcept;

    [[nodiscard]]
    Eigen::Vector3i position_to_cell(
        const Eigen::Vector3f& position
    ) const;

private:
    int resolution_x_;
    int resolution_y_;
    int resolution_z_;

    Eigen::Vector3f min_;
    Eigen::Vector3f max_;

    int visibility_order_;

    std::vector<Visibility> cells_;
};

} // namespace render