/**
 * @file spherical_grid.hpp
 * @brief Spherical Grid mathematical layout.
 * 
 * Determines the exact number of Theta and Phi grid points required to accurately 
 * represent an SH product of order 't'. Generates the 3D directional vectors 
 * and quadrature weights for each point.
 */
#pragma once
#include <vector>
#include <cstddef>
#include <stdexcept>
#include<numbers>
#include <cmath>
#include<types/constants.hpp>
namespace sylph{

    class Spherical_grid{
        public:
        Spherical_grid(int nTheta,int nPhi, const std::vector<double>&nodes)
        :nTheta_(nTheta),nPhi_(nPhi),theta_(nTheta),
          phi_(nPhi),
          values_(static_cast<std::size_t>(nTheta) *
                  static_cast<std::size_t>(nPhi), 0.0)
                  {

                    if (nTheta <= 0 || nPhi <= 0) {
            throw std::invalid_argument(
                "Grid dimensions must be positive"
            );
        }

        if (static_cast<int>(nodes.size()) != nTheta) {
            throw std::invalid_argument(
                "Invalid number of Gauss-Legendre nodes"
            );
        }

           for (int i = 0; i < nTheta_; ++i) {
            theta_[i] = std::acos(nodes[i]);
        }

        for (int j = 0; j < nPhi_; ++j) {
    phi_[j] =
        2.0 * pi
        * static_cast<double>(j)
        / static_cast<double>(nPhi_);
}

        }


        int theta_size() const {
        return nTheta_;
    }

    int phi_size() const {
        return nPhi_;
    }

    double theta(int i) const {
        return theta_.at(i);
    }

    double phi(int j) const {
        return phi_.at(j);
    }

    double& operator()(int i, int j) {
        return values_[index(i, j)];
    }

    double operator()(int i, int j) const {
        return values_[index(i, j)];
    }
    
    void multiply_inplace(const Spherical_grid& rhs){
        if (theta_size() != rhs.theta_size() ||
        phi_size() != rhs.phi_size())
    {
        throw std::invalid_argument(
            "Spherical_grid dimensions must match"
        );
    }

    for (int i = 0; i < theta_size(); ++i)
    {
        for (int j = 0; j < phi_size(); ++j)
        {
            (*this)(i, j) *= rhs(i, j);
        }
    }
    }

    static Spherical_grid multiply(
    const Spherical_grid& a,
    const Spherical_grid& b
){
      Spherical_grid result = a;
    result.multiply_inplace(b);
    return result;
}
    private:

    int nTheta_;
    int nPhi_;

    std::vector<double> theta_;
    std::vector<double> phi_;
    std::vector<double> values_;

    std::size_t index(int i, int j) const{
        return static_cast<std::size_t>(i) * nPhi_ + j;
    }
  
    };

template<typename Function>
void sample(
    Spherical_grid& grid,
    Function&& f
) {
    for (int i = 0; i < grid.theta_size(); ++i) {
        for (int j = 0; j < grid.phi_size(); ++j) {

            grid(i, j) =
                f(
                    grid.theta(i),
                    grid.phi(j)
                );
        }
    }
}

}