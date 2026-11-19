/**
 * @file polar_basis_cache.hpp
 * @brief Cached Normalized Legendre Polar Basis Table (Optimization #8 - Level 2).
 *
 * Precomputes and caches P~_l^m(x_i) for all northern/equator Gauss-Legendre
 * nodes i in [0, (nTheta + 1) / 2) and all 0 <= m <= l < order in a compact,
 * column-major contiguous layout:
 *
 *   col(i, m)[k] = P~_{m + k}^m(x_i),   k = 0 .. order - m - 1
 *
 * Eliminates:
 *  - Redundant 3-term Legendre recurrence evaluations across multiple inputs
 *    and inverse passes in SHGridTransform::product().
 *  - Serial recurrence loop-carried FMA dependencies, exposing full CPU ILP.
 *  - All per-call acos / cos / sin evaluations.
 */
#pragma once

#include "sylph/gauss_legendre.hpp"
#include "sylph/legendre.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace sylph {

class PolarBasisCache {
public:
    struct Table {
        int nTheta = 0;
        int nHalf  = 0;
        int order  = 0;
        std::size_t tri_size = 0;

        std::vector<std::size_t> col_offset;
        std::vector<double> values; // [nHalf][tri_size], contiguous per (i, m, l-m)

        const double* col(int i, int m) const noexcept {
            return values.data()
                 + static_cast<std::size_t>(i) * tri_size
                 + col_offset[static_cast<std::size_t>(m)];
        }
    };

    static const Table& get(int nTheta, int order) {
        if (nTheta <= 0 || order <= 0) {
            throw std::invalid_argument("PolarBasisCache requires positive nTheta and order");
        }

        const std::uint64_t key =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(nTheta)) << 32) |
             static_cast<std::uint64_t>(static_cast<std::uint32_t>(order));

        thread_local std::unordered_map<std::uint64_t, Table> cache;
        auto it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }

        auto [inserted, _] = cache.emplace(key, build_table(nTheta, order));
        return inserted->second;
    }

private:
    static Table build_table(int nTheta, int order) {
        const auto& rule = gauss_legendre::get(nTheta);
        const auto& rec  = Legendre::get_recurrence_table(order);

        Table tbl;
        tbl.nTheta   = nTheta;
        tbl.nHalf    = (nTheta + 1) / 2;
        tbl.order    = order;
        tbl.tri_size = static_cast<std::size_t>(order) * static_cast<std::size_t>(order + 1) / 2;
        tbl.col_offset = rec.col_offset;
        tbl.values.resize(static_cast<std::size_t>(tbl.nHalf) * tbl.tri_size);

        const double* diag  = rec.diag.data();
        const double* a_tbl = rec.a.data();
        const double* b_tbl = rec.b.data();
        const double p00    = Legendre::inv_sqrt_4pi();

        for (int i = 0; i < tbl.nHalf; ++i) {
            const double x  = rule.nodes[static_cast<std::size_t>(i)];
            const double st = rule.sin_theta[static_cast<std::size_t>(i)];
            double* node_base = tbl.values.data() + static_cast<std::size_t>(i) * tbl.tri_size;

            double pmm = p00;
            for (int m = 0; m < order; ++m) {
                if (m > 0) {
                    pmm = -diag[static_cast<std::size_t>(m)] * st * pmm;
                }

                const std::size_t off = tbl.col_offset[static_cast<std::size_t>(m)];
                double* dst = node_base + off;
                const double* a_col = a_tbl + off;
                const double* b_col = b_tbl + off;
                const int len = order - m;

                dst[0] = pmm;
                if (len > 1) {
                    dst[1] = a_col[1] * x * pmm;
                    for (int k = 2; k < len; ++k) {
                        dst[k] = a_col[k] * x * dst[k - 1] - b_col[k] * dst[k - 2];
                    }
                }
            }
        }

        return tbl;
    }
};

} // namespace sylph
