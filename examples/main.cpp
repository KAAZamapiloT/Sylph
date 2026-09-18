#include "sylph/spherical_harmonics.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

bool approx(double a, double b, double eps = 1e-10) {
    return std::abs(a - b) < eps;
}

int main() {

    sylph::SphericalHarmonics sh;

    double theta = 0.7;
    double phi = 1.1;

    auto values = sh.basis(2, theta, phi);

    assert(values.size() == 9);

    // Verify that basis() actually produces
    // the same values as evaluate().

    for (int l = 0; l <= 2; ++l) {
        for (int m = -l; m <= l; ++m) {

            double expected =
                sh.evaluate(l, m, theta, phi);

            double got =
                values[sylph::SphericalHarmonics::index(l, m)];

            assert(approx(got, expected));
        }
    }
    for (int l = 0; l <= 2; ++l) {
        for (int m = -l; m <= l; ++m) {

            int i = sylph::SphericalHarmonics::index(l, m);

            std::cout
                << "index=" << i
                << " (" << l << "," << m << ")"
                << " = " << values[i]
                << '\n';
        }
    }
    std::cout << "Basis test passed\n";
}
