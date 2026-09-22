#pragma once

#include "sylph/sh_coefficients.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace render
{

enum class LightChannel : std::uint8_t
{
    Red = 0,
    Green = 1,
    Blue = 2
};

class Lighting
{
public:
    explicit Lighting(int order);

    // SH order
    [[nodiscard]]
    int order() const;

    // Access one color channel's SH coefficients.
    [[nodiscard]]
    const sylph::SHCoefficients&
    coefficients(LightChannel channel) const;

    [[nodiscard]]
    sylph::SHCoefficients&
    coefficients(LightChannel channel);

    // Replace one channel.
    void set_coefficients(
        LightChannel channel,
        const sylph::SHCoefficients& coefficients
    );

    // Replace all channels.
    void set_coefficients(
        const std::array<
            sylph::SHCoefficients,
            3
        >& coefficients
    );

    // Reset all channels to zero.
    void clear();

private:
    static std::size_t
    channel_index(LightChannel channel);

private:
    std::array<
        sylph::SHCoefficients,
        3
    > coefficients_;
};

} // namespace render