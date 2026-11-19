/**
 * @file gauss_legendre.hpp
 * @brief Cached Gauss-Legendre Quadrature Generator (Optimization #8).
 *
 * Computes optimal Gauss-Legendre nodes x_i = cos(theta_i) and weights w_i
 * for integrating functions over [-1, 1], cached thread-locally per order n
 * so repeated SH products incur zero Newton iterations and zero trig calls.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <types/constants.hpp>

namespace sylph {

class gauss_legendre {
public:
    struct Rule {
        int n = 0;
        std::vector<double> nodes;      // x_i = cos(theta_i)
        std::vector<double> weights;    // w_i
        std::vector<double> theta;      // acos(x_i)
        std::vector<double> sin_theta;  // sqrt(1 - x_i^2)
    };

    /**
     * @brief Returns a cached Gauss-Legendre rule of order n (O(1) after first call).
     */
    static const Rule& get(int n) {
        if (n <= 0) {
            throw std::invalid_argument(
                "Gauss-Legendre order must be positive"
            );
        }

        thread_local std::unordered_map<int, Rule> cache;
        auto it = cache.find(n);
        if (it != cache.end()) {
            return it->second;
        }

        auto [inserted, _] = cache.emplace(n, compute_rule(n));
        return inserted->second;
    }

    /**
     * @brief Uncached computation retained solely for before/after #8 benchmarking.
     */
    static Rule compute_uncached(int n) {
        if (n <= 0) {
            throw std::invalid_argument(
                "Gauss-Legendre order must be positive"
            );
        }
        return compute_rule(n);
    }

    explicit gauss_legendre(int n)
        : nodes_(get(n).nodes),
          weights_(get(n).weights)
    {}

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
                ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
            p0 = p1;
            p1 = p2;
        }

        // P'_n(x)
        double derivative =
            n * (p0 - x * p1) / (1.0 - x * x);

        return {p1, derivative};
    }

    static Rule compute_rule(int n)
    {
        Rule rule;
        rule.n = n;
        rule.nodes.resize(static_cast<std::size_t>(n));
        rule.weights.resize(static_cast<std::size_t>(n));
        rule.theta.resize(static_cast<std::size_t>(n));
        rule.sin_theta.resize(static_cast<std::size_t>(n));

        const int half = (n + 1) / 2;

        for (int i = 0; i < half; ++i) {
            // Initial Newton guess
            double x = std::cos(pi * (i + 0.75) / (n + 0.5));

            for (int iteration = 0; iteration < 100; ++iteration) {
                auto [pn, dpn] = legendre_and_derivative(n, x);
                double next = x - pn / dpn;
                if (std::abs(next - x) < 1e-14) {
                    x = next;
                    break;
                }
                x = next;
            }

            auto [pn, dpn] = legendre_and_derivative(n, x);
            double weight = 2.0 / ((1.0 - x * x) * dpn * dpn);

            const double th = std::acos(std::clamp(x, -1.0, 1.0));
            const double st = std::sqrt(std::max(0.0, 1.0 - x * x));

            // Positive root (northern hemisphere)
            rule.nodes[static_cast<std::size_t>(i)]     = x;
            rule.weights[static_cast<std::size_t>(i)]   = weight;
            rule.theta[static_cast<std::size_t>(i)]     = th;
            rule.sin_theta[static_cast<std::size_t>(i)] = st;

            // Mirrored negative root (southern hemisphere)
            const int j = n - 1 - i;
            if (j != i) {
                rule.nodes[static_cast<std::size_t>(j)]     = -x;
                rule.weights[static_cast<std::size_t>(j)]   = weight;
                rule.theta[static_cast<std::size_t>(j)]     = pi - th;
                rule.sin_theta[static_cast<std::size_t>(j)] = st;
            } else {
                // Exact equator node for odd n
                rule.nodes[static_cast<std::size_t>(i)]     = 0.0;
                rule.theta[static_cast<std::size_t>(i)]     = 0.5 * pi;
                rule.sin_theta[static_cast<std::size_t>(i)] = 1.0;
            }
        }

        return rule;
    }
};

} // namespace sylph
