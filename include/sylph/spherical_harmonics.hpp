#pragma once
#include"legendre.hpp"

#include<cmath>
#include <vector>
#include<stdexcept>
#include<sylph/sh_coefficients.hpp>
namespace sylph{
class SphericalHarmonics {
public:
    double evaluate(
            int l,
            int m,
            double theta,
            double phi
        )const{
            validate(l, m);
            return legendre_.sh(l, m, theta, phi);
        }


        static constexpr int index(int l, int m) {
            return l * l + l + m;
        }

        static constexpr int coefficient_count(int order) {
            return (order) * (order);
        }

    std::vector<double> basis(
        int order,
        double theta,
        double phi
    ) {

        std::vector<double> result(coefficient_count(order));

        for(int l=0;l<order;++l) {
            for(int m=-l;m<=l;++m){
                result[index(l, m)] = evaluate(l, m, theta, phi);
            }
        }
        return result;
    }

    double reconstruct(
        const SHCoefficients& coefficients,
        double theta,
        double phi
    )const {
        int siz=coefficients.order();
        
        double result=0;
        for(int l=0;l<siz;++l){
            for(int m=-l;m<=l;++m){
                  result+=coefficients(l,m)*evaluate(l,m,theta,phi);
            }

        }
        return result;
    }
    template<typename Function>
    SHCoefficients project(Function&&function,int order,
        int thetaSamples,int phiSamples)
        {
        
            if (order <= 0 ||
            thetaSamples <= 0 ||
            phiSamples <= 0) {
            throw std::invalid_argument(
                "Invalid projection parameters"
            );
        }

            SHCoefficients result(order);

            const double dTheta =
    std::numbers::pi_v<double> / thetaSamples;

const double dPhi =
    2.0 * std::numbers::pi_v<double> / phiSamples;


            for(int l=0;l<order;++l){
                for(int m=-l;m<=l;++m){

                    double coefficient=0.0;

                    for(int i=0;i<thetaSamples;++i){
                        double theta=(i+0.5)*dTheta;
                        for(int j=0;j<phiSamples;++j){
                            double phi=(j+0.5)*dPhi;
                            coefficient+=function(theta,phi)*
                            evaluate(l,m,theta,phi)*std::sin(theta);
                        }
                    }
                    coefficient*=dTheta*dPhi;
                    result(l,m)=coefficient;
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
