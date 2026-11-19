/**
 * @file sh_coefficients.hpp
 * @brief Core Spherical Harmonic Coefficient Data Structure.
 * 
 * Stores a vector of SH coefficients up to an arbitrary band limit (order). 
 * Provides operators for scaling, adding, and managing the memory of the 
 * frequency-domain data.
 */
#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

#include<types/constants.hpp>

namespace sylph{


    class SHCoefficients  {
       public:
        explicit SHCoefficients(int order):order_(order),data_((order) * (order), 0.0){

        }
    double& operator()(int l, int m) {
        validate(l, m);
        return data_[index(l, m)];
    }

    double operator()(int l, int m) const {
        validate(l, m);
        return data_[index(l, m)];
    }

    int order() const {
        return order_;
    }

    static constexpr int index(int l, int m) {
        return l * l + l + m;
    }

    static constexpr int size_for(int order) {
        return (order) * (order);
    }

    const double* data() const noexcept {
        return data_.data();
    }

    double* data() noexcept {
        return data_.data();
    }

    /**
     * @brief Returns true if all non-zonal coefficients (m != 0) satisfy |c(l,m)| <= eps.
     *
     * Exits on the first non-zero m != 0 entry (O(1) for dense general SH functions).
     */
    bool is_zonal(double eps = 0.0) const noexcept {
        for (int l = 1; l < order_; ++l) {
            const int base = l * l + l;
            for (int m = -l; m < 0; ++m) {
                if (std::abs(data_[static_cast<std::size_t>(base + m)]) > eps) {
                    return false;
                }
            }
            for (int m = 1; m <= l; ++m) {
                if (std::abs(data_[static_cast<std::size_t>(base + m)]) > eps) {
                    return false;
                }
            }
        }
        return true;
    }

        private:
        int order_;
    std::vector<double> data_;

    void validate(int l, int m) const {
        if (l < 0 || l >= order_ || std::abs(m) > l) {
            throw std::out_of_range("Invalid (l,m)");
        }
    }
    
    };
}