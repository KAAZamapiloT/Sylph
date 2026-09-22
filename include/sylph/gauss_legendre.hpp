/**
 * @file gauss_legendre.hpp
 * @brief Gauss-Legendre Quadrature Generator.
 * 
 * Computes the optimal node angles and weights for integrating functions over a sphere. 
 * Used heavily by the Spherical Grid algorithm to accurately sum the triple/multiple 
 * products without losing energy.
 */
#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>
#include<types/constants.hpp>
namespace sylph {

class gauss_legendre {
public:

    explicit gauss_legendre(int n)
        : nodes_(n),
          weights_(n)
    {
        if (n <= 0) {
            throw std::invalid_argument(
                "Gauss-Legendre order must be positive"
            );
        }

        compute(n);
    }

    const std::vector<double>& nodes() const {
        return nodes_;
    }

    const std::vector<double>& weights() const {
        return weights_;
    }

private:

    std::vector<double> nodes_;
    std::vector<double> weights_;

    static std::pair<double, double>
    legendre_and_derivative(int n, double x)
    {
        double p0 = 1.0;
        double p1 = x;

        if (n == 0) {
            return {p0, 0.0};
        }

        if (n == 1) {
            return {p1, 1.0};
        }

        for (int k = 2; k <= n; ++k) {

            double p2 =
                (
                    (2.0 * k - 1.0) * x * p1
                    - (k - 1.0) * p0
                ) / k;

            p0 = p1;
            p1 = p2;
        }

        // P'_n(x)
        double derivative =
            n * (p0 - x * p1) /
            (1.0 - x * x);

        return {p1, derivative};
    }

    void compute(int n)
    {
        const int half = (n + 1) / 2;

        for (int i = 0; i < half; ++i) {

            // Initial Newton guess
            double x =
                std::cos(
                    pi * (i + 0.75) /
                    (n + 0.5)
                );

            for (int iteration = 0;
                 iteration < 100;
                 ++iteration)
            {
                auto [pn, dpn] =
                    legendre_and_derivative(n, x);

                double next =
                    x - pn / dpn;

                if (std::abs(next - x) < 1e-14) {
                    x = next;
                    break;
                }

                x = next;
            }

            auto [pn, dpn] =
                legendre_and_derivative(n, x);

            double weight =
                2.0 /
                (
                    (1.0 - x * x)
                    * dpn * dpn
                );

            // Positive root
            nodes_[i] = x;
            weights_[i] = weight;

            // Mirrored negative root
            int j = n - 1 - i;

            nodes_[j] = -x;
            weights_[j] = weight;
        }
    }
};

}