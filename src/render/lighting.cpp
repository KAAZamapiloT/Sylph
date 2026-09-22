#include "render/lighting.hpp"

#include <stdexcept>

namespace render
{

Lighting::Lighting(int order)
    : coefficients_{
          sylph::SHCoefficients(order),
          sylph::SHCoefficients(order),
          sylph::SHCoefficients(order)
      }
{
    if (order <= 0) {
        throw std::invalid_argument(
            "Lighting order must be positive."
        );
    }

    clear();
}

// ------------------------------------------------------------
// Order
// ------------------------------------------------------------

int Lighting::order() const
{
    return coefficients_[0].order();
}

// ------------------------------------------------------------
// Channel access
// ------------------------------------------------------------

const sylph::SHCoefficients&
Lighting::coefficients(
    LightChannel channel) const
{
    return coefficients_[channel_index(channel)];
}

sylph::SHCoefficients&
Lighting::coefficients(
    LightChannel channel)
{
    return coefficients_[channel_index(channel)];
}

// ------------------------------------------------------------
// Set one channel
// ------------------------------------------------------------

void Lighting::set_coefficients(
    LightChannel channel,
    const sylph::SHCoefficients& coefficients)
{
    if (coefficients.order() != order()) {
        throw std::invalid_argument(
            "Lighting channel order does not match "
            "the lighting order."
        );
    }

    coefficients_[channel_index(channel)] =
        coefficients;
}

// ------------------------------------------------------------
// Set all channels
// ------------------------------------------------------------

void Lighting::set_coefficients(
    const std::array<
        sylph::SHCoefficients,
        3
    >& coefficients)
{
    const int expected_order =
        coefficients_[0].order();

    for (const auto& channel : coefficients) {
        if (channel.order() != expected_order) {
            throw std::invalid_argument(
                "All lighting channels must have "
                "the same SH order."
            );
        }
    }

    coefficients_ = coefficients;
}

// ------------------------------------------------------------
// Clear
// ------------------------------------------------------------

void Lighting::clear()
{
    for (auto& channel : coefficients_) {
        for (int l = 0; l < channel.order(); ++l) {
            for (int m = -l; m <= l; ++m) {
                channel(l, m) = 0.0;
            }
        }
    }
}

// ------------------------------------------------------------
// Channel → array index
// ------------------------------------------------------------

std::size_t Lighting::channel_index(
    LightChannel channel)
{
    const auto index =
        static_cast<std::size_t>(channel);

    if (index >= 3) {
        throw std::invalid_argument(
            "Invalid lighting channel."
        );
    }

    return index;
}

} // namespace render