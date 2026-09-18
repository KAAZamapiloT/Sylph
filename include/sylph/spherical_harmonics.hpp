#pragma once
#include"legendre.hpp"

#include<cmath>
#include <vector>
#include<stdexcept>

namespace sylph{
class SphericalHarmonics {
public:
    double evaluate(
            int l,
            int m,
            double theta,
            double phi
        ){
            validate(l, m);
            return legendre_.sh(l, m, theta, phi);
        }


        static constexpr int index(int l, int m) {
            return l * l + l + m;
        }

        static constexpr int coefficient_count(int lmax) {
            return (lmax + 1) * (lmax + 1);
        }

    std::vector<double> basis(
        int lmax,
        double theta,
        double phi
    ) {

        std::vector<double> result(coefficient_count(lmax));

        for(int l=0;l<=lmax;++l) {
            for(int m=-l;m<=l;++m){
                result[index(l, m)] = evaluate(l, m, theta, phi);
            }
        }
        return result;
    }

private:
    Legendre legendre_;
    static void validate(int l, int m) {
            if (l < 0 || std::abs(m) > l) {
                throw std::invalid_argument(
                    "Invalid spherical harmonic indices"
                );
            }
        }
};

}
