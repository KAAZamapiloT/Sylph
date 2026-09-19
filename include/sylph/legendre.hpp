#include "iostream"

#include <cmath>
#include <numbers>
#include <stdexcept>


namespace sylph {
constexpr double pi = std::numbers::pi_v<double>;
    class Legendre {
    public:
        double sh(int l, int m, double theta, double phi){
            if(m==0){
                 return Nl(0, l) * Pl(0, l, std::cos(theta));
            }else if(m>0){
                return std::sqrt(2)*Nl(m,l)*Pl(m,l,cos(theta))*std::cos(m*phi);
            }else if(m<0){
                return std::sqrt(2)*Nl(std::abs(m),l)*Pl(std::abs(m),l,cos(theta))*std::sin(std::abs(m)*phi);
            }
             return 0.0;
        }
    private:

    double Nl(int m,int l){
        double A = (2.0 * l + 1.0) / (4.0 * pi);
        double B = factorial_ratio(l - m, l + m);
        double C=std::sqrt(A*B);
        return ((m%2==0)?1:-1)*C;
    }

    double Pl(int m, int l, double x){
        const double factor =
                    std::pow(1.0 - x * x, 0.5 * m);

                const int lim = (l - m) / 2;

                double result = 0.0;

                for (int k = 0; k <= lim; ++k) {

                    const double power =
                        std::pow(x, l - 2 * k - m);

                    const double coefficient =
                        std::tgamma(2 * l - 2 * k + 1.0) /
                        (
                            std::pow(2.0, l)
                            * std::tgamma(k + 1.0)
                            * std::tgamma(l - k + 1.0)
                            * std::tgamma(l - 2 * k - m + 1.0)
                        );

                    result += (k % 2 == 0 ? 1.0 : -1.0)
                            * power
                            * coefficient;
                }

                return ((m%2==0)?1:-1)*factor * result;
    }
    double factorial_ratio(int lo, int hi) {
        double result = 1.0;
        for (int i = lo + 1; i <= hi; ++i)
            result /= i;
        return result;
    }
    };



}
