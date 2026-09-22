#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include <vector>

double factorial(int n) {
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

int index(int l, int m) { return l * l + l + m; }

int main() {
    int order = 8;
    int count = order * order;
    std::cout << "const float uSHNorm[" << count << "] = float[](\n";
    for (int l = 0; l < order; ++l) {
        for (int m = -l; m <= l; ++m) {
            double A = (2.0 * l + 1.0) / (4.0 * std::numbers::pi_v<double>);
            double B = factorial(l - std::abs(m)) / factorial(l + std::abs(m));
            double Nl = std::sqrt(A * B);
            if (m != 0) Nl *= std::sqrt(2.0);
            
            std::cout << std::scientific << std::setprecision(7) << Nl << "f, ";
        }
        std::cout << "\n";
    }
    std::cout << ");\n";
    return 0;
}
