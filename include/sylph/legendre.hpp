/**
 * @file legendre.hpp
 * @brief Normalized Associated Legendre Polynomial mathematics.
 *
 * Uses a directly normalized 3-term recurrence relation (zero std::tgamma calls,
 * zero factorial overflow even at high orders N >= 64) with precomputed
 * column-major recurrence coefficients cached per order.
 *
 * Convention: includes the Condon-Shortley phase (-1)^m and real SH
 * normalization factor:
 *   N_l^m = sqrt( (2l + 1) / (4*pi) * (l - m)! / (l + m)! )
 *   P~_l^m(x) = N_l^m * P_l^m(x)
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <types/constants.hpp>

namespace sylph {

class Legendre {
public:
    /**
     * @brief Precomputed recurrence coefficients for all 0 <= m <= l < order.
     *
     * Stored column-major (contiguous in l = m .. order - 1 for each m):
     *   - diag[m]         = sqrt((2m + 1) / (2m))           for m >= 1
     *   - a[col(m) + l-m] = sqrt((4l^2 - 1) / (l^2 - m^2))  for l >= m + 1
     *   - b[col(m) + l-m] = a(l, m) / a(l - 1, m)           for l >= m + 2
     */
    struct RecurrenceTable {
        int order = 0;
        std::vector<std::size_t> col_offset;
        std::vector<double> diag;
        std::vector<double> a;
        std::vector<double> b;
    };

    static double inv_sqrt_4pi() noexcept {
        static const double val = 1.0 / std::sqrt(4.0 * pi);
        return val;
    }

    static const RecurrenceTable& get_recurrence_table(int order) {
        if (order <= 0) {
            throw std::invalid_argument("Legendre recurrence order must be positive");
        }

        thread_local std::unordered_map<int, RecurrenceTable> cache;
        auto it = cache.find(order);
        if (it != cache.end()) {
            return it->second;
        }

        auto [inserted, _] = cache.emplace(order, build_recurrence_table(order));
        return inserted->second;
    }

    /**
     * @brief Real Spherical Harmonic Y_l^m(theta, phi) via normalized 3-term recurrence.
     */
    double sh(int l, int m, double theta, double phi) const {
        const int am = std::abs(m);
        const double x = std::cos(theta);
        const double s = std::max(0.0, std::sin(theta));
        const double p_norm = normalized_associated_cs(l, am, x, s);

        if (m == 0) {
            return p_norm;
        }
        if (m > 0) {
            return std::numbers::sqrt2_v<double> * p_norm * std::cos(static_cast<double>(m) * phi);
        }
        return std::numbers::sqrt2_v<double> * p_norm * std::sin(static_cast<double>(am) * phi);
    }

    /**
     * @brief Directly normalized associated Legendre P~_l^m(x) = N_l^m * P_l^m(x).
     *
     * Computed via stable 3-term recurrence (no factorials, no tgamma).
     * Requires 0 <= m <= l and x in [-1, 1].
     */
    double normalized_associated(int l, int m, double x) const {
        const double s = std::sqrt(std::max(0.0, 1.0 - x * x));
        return normalized_associated_cs(l, m, x, s);
    }

    /**
     * @brief Directly normalized P~_l^m(cos theta) given precomputed (x = cos theta, s = sin theta).
     */
    double normalized_associated_cs(int l, int m, double x, double s) const {
        if (m < 0 || l < m) {
            throw std::invalid_argument("Invalid associated Legendre indices (l, m)");
        }

        // 1. Diagonal seed: P~_0^0 = 1 / sqrt(4*pi),
        //    P~_k^k = -sqrt((2k + 1) / (2k)) * s * P~_{k-1}^{k-1}
        double pmm = inv_sqrt_4pi();
        for (int k = 1; k <= m; ++k) {
            const double dk = std::sqrt(static_cast<double>(2 * k + 1) / static_cast<double>(2 * k));
            pmm = -dk * s * pmm;
        }

        if (l == m) {
            return pmm;
        }

        // 2. Sub-diagonal: P~_{m+1}^m = sqrt(2m + 3) * x * P~_m^m
        double a_prev = std::sqrt(static_cast<double>(2 * m + 3));
        double p1 = a_prev * x * pmm;
        if (l == m + 1) {
            return p1;
        }

        // 3. Three-term recurrence for k = m + 2 .. l:
        //    P~_k^m = a_{k,m} * x * P~_{k-1}^m - (a_{k,m} / a_{k-1,m}) * P~_{k-2}^m
        double p0 = pmm;
        const double m2 = static_cast<double>(m) * static_cast<double>(m);
        for (int k = m + 2; k <= l; ++k) {
            const double kd = static_cast<double>(k);
            const double a_curr = std::sqrt((4.0 * kd * kd - 1.0) / (kd * kd - m2));
            const double b_curr = a_curr / a_prev;
            const double p2 = a_curr * x * p1 - b_curr * p0;
            p0 = p1;
            p1 = p2;
            a_prev = a_curr;
        }

        return p1;
    }

    /**
     * @brief Unnormalized associated Legendre P_l^m(x) with Condon-Shortley phase (-1)^m,
     *        computed via 3-term recurrence (zero tgamma calls).
     */
    double associated(int l, int m, double x) const {
        if (m < 0 || l < m) {
            throw std::invalid_argument("Invalid associated Legendre indices (l, m)");
        }

        const double s = std::sqrt(std::max(0.0, 1.0 - x * x));
        double pmm = 1.0;
        for (int k = 1; k <= m; ++k) {
            pmm *= -static_cast<double>(2 * k - 1) * s;
        }

        if (l == m) {
            return pmm;
        }

        double p1 = static_cast<double>(2 * m + 1) * x * pmm;
        if (l == m + 1) {
            return p1;
        }

        double p0 = pmm;
        for (int k = m + 2; k <= l; ++k) {
            const double p2 =
                (static_cast<double>(2 * k - 1) * x * p1 -
                 static_cast<double>(k + m - 1) * p0) /
                static_cast<double>(k - m);
            p0 = p1;
            p1 = p2;
        }

        return p1;
    }

    double normalization(int l, int m) const {
        return Nl(m, l);
    }

    /**
     * @brief Legacy closed-form tgamma evaluator retained solely for before/after benchmarking.
     */
    double normalized_tgamma_legacy(int l, int m, double x) const {
        return Nl(m, l) * Pl_tgamma_legacy(m, l, x);
    }

private:
    static RecurrenceTable build_recurrence_table(int order) {
        RecurrenceTable tbl;
        tbl.order = order;
        tbl.col_offset.resize(static_cast<std::size_t>(order));
        tbl.diag.resize(static_cast<std::size_t>(order), 1.0);

        const std::size_t total =
            static_cast<std::size_t>(order) * static_cast<std::size_t>(order + 1) / 2;
        tbl.a.assign(total, 0.0);
        tbl.b.assign(total, 0.0);

        for (int m = 1; m < order; ++m) {
            tbl.diag[static_cast<std::size_t>(m)] =
                std::sqrt(static_cast<double>(2 * m + 1) / static_cast<double>(2 * m));
        }

        std::size_t offset = 0;
        for (int m = 0; m < order; ++m) {
            tbl.col_offset[static_cast<std::size_t>(m)] = offset;
            const double m2 = static_cast<double>(m) * static_cast<double>(m);

            double a_prev = 0.0;
            for (int l = m + 1; l < order; ++l) {
                const double ld = static_cast<double>(l);
                const double a_curr = std::sqrt((4.0 * ld * ld - 1.0) / (ld * ld - m2));
                const std::size_t idx = offset + static_cast<std::size_t>(l - m);
                tbl.a[idx] = a_curr;
                if (l >= m + 2) {
                    tbl.b[idx] = a_curr / a_prev;
                }
                a_prev = a_curr;
            }

            offset += static_cast<std::size_t>(order - m);
        }

        return tbl;
    }

    double Nl(int m, int l) const {
        const double A = (2.0 * l + 1.0) / (4.0 * pi);
        const double B = factorial_ratio(l - m, l + m);
        return std::sqrt(A * B);
    }

    double Pl_tgamma_legacy(int m, int l, double x) const {
        const double factor = std::pow(1.0 - x * x, 0.5 * m);
        const int lim = (l - m) / 2;
        double result = 0.0;

        for (int k = 0; k <= lim; ++k) {
            const double power = std::pow(x, l - 2 * k - m);
            const double coefficient =
                std::tgamma(2 * l - 2 * k + 1.0) /
                (std::pow(2.0, l) *
                 std::tgamma(k + 1.0) *
                 std::tgamma(l - k + 1.0) *
                 std::tgamma(l - 2 * k - m + 1.0));

            result += (k % 2 == 0 ? 1.0 : -1.0) * power * coefficient;
        }

        return ((m % 2 == 0) ? 1.0 : -1.0) * factor * result;
    }

    double factorial_ratio(int lo, int hi) const {
        double result = 1.0;
        for (int i = lo + 1; i <= hi; ++i) {
            result /= static_cast<double>(i);
        }
        return result;
    }
};

} // namespace sylph
