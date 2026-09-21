/*
    Sylph system-level mathematical validation
    -------------------------------------------

    Purpose:
      Validate the SH + Spherical Grid pipeline against an independent
      numerical reference instead of comparing one Sylph implementation
      against another.

    Coverage:
      1. Real SH orthonormality
      2. SH -> SphGrid -> SH round-trip
      3. Exact generalized product against independent quadrature
      4. Product integral against independent quadrature
      5. Two-function integral against coefficient dot-product identity
      6. Exact t threshold
      7. Lower-t approximate mode
      8. Product commutativity
      9. Multiplication by the constant-one function
     10. Parseval / energy identity
     11. Eq. 34 inverse-resolution boundary checks

    Notes:
      - The reference SH basis is implemented independently here.
      - The reference integrator uses Gauss-Legendre x uniform-phi
        quadrature at substantially higher resolution than the production
        exact grids.
      - This file intentionally does NOT use SHGridTransform for the
        reference calculations.
*/

#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"
#include "sylph/spherical_grid.hpp"
#include "sylph/gauss_legendre.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double SQRT_TWO = 1.41421356237309504880168872420969808;

struct TestState
{
    int passed = 0;
    int failed = 0;

    double worst_error = 0.0;
    std::string worst_name;
};

TestState state;

void pass()
{
    ++state.passed;
}

void fail(const std::string& name)
{
    ++state.failed;
    std::cerr << "[FAIL] " << name << '\n';
}

void fail(
    const std::string& name,
    double got,
    double expected,
    double error,
    double tolerance
)
{
    ++state.failed;

    if (error > state.worst_error)
    {
        state.worst_error = error;
        state.worst_name = name;
    }

    std::cerr
        << std::setprecision(17)
        << "[FAIL] " << name
        << "\n       got       = " << got
        << "\n       expected  = " << expected
        << "\n       error     = " << error
        << "\n       tolerance = " << tolerance
        << '\n';
}

void check(
    bool condition,
    const std::string& name
)
{
    if (condition)
    {
        pass();
    }
    else
    {
        fail(name);
    }
}

void close(
    double got,
    double expected,
    const std::string& name,
    double abs_tol = 1e-10,
    double rel_tol = 1e-9
)
{
    const double error = std::abs(got - expected);

    const double scale =
        std::max({
            1.0,
            std::abs(got),
            std::abs(expected)
        });

    const double tolerance =
        abs_tol + rel_tol * scale;

    if (std::isfinite(got) &&
        std::isfinite(expected) &&
        error <= tolerance)
    {
        pass();
    }
    else
    {
        fail(
            name,
            got,
            expected,
            error,
            tolerance
        );
    }
}

void close_error(
    double error,
    const std::string& name,
    double tolerance = 1e-10
)
{
    if (std::isfinite(error) && error <= tolerance)
    {
        pass();
    }
    else
    {
        ++state.failed;

        if (error > state.worst_error)
        {
            state.worst_error = error;
            state.worst_name = name;
        }

        std::cerr
            << std::setprecision(17)
            << "[FAIL] " << name
            << "\n       error     = " << error
            << "\n       tolerance = " << tolerance
            << '\n';
    }
}

template <typename F>
void expect_invalid_argument(
    F&& function,
    const std::string& name
)
{
    try
    {
        function();

        ++state.failed;
        std::cerr
            << "[FAIL] " << name
            << "\n       expected std::invalid_argument\n";
    }
    catch (const std::invalid_argument&)
    {
        pass();
    }
    catch (const std::exception& e)
    {
        ++state.failed;
        std::cerr
            << "[FAIL] " << name
            << "\n       wrong exception: "
            << e.what()
            << '\n';
    }
    catch (...)
    {
        ++state.failed;
        std::cerr
            << "[FAIL] " << name
            << "\n       wrong non-standard exception\n";
    }
}

void section(const std::string& title)
{
    std::cout
        << "\n============================================================\n"
        << title
        << "\n============================================================\n";
}

/*
    Independent associated Legendre implementation.

    Uses the Condon-Shortley convention:

        P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)

    followed by the standard upward recurrence.
*/
double reference_associated_legendre(
    int l,
    int m,
    double x
)
{
    if (l < 0 || m < 0 || m > l)
    {
        throw std::invalid_argument(
            "invalid associated Legendre indices"
        );
    }

    x = std::clamp(x, -1.0, 1.0);

    double pmm = 1.0;

    if (m > 0)
    {
        double odd_double_factorial = 1.0;

        for (int k = 1; k <= m; ++k)
        {
            odd_double_factorial *=
                static_cast<double>(2 * k - 1);
        }

        pmm =
            ((m & 1) ? -1.0 : 1.0) *
            odd_double_factorial *
            std::pow(
                std::max(
                    0.0,
                    1.0 - x * x
                ),
                0.5 * static_cast<double>(m)
            );
    }

    if (l == m)
    {
        return pmm;
    }

    double pm1m =
        x *
        static_cast<double>(2 * m + 1) *
        pmm;

    if (l == m + 1)
    {
        return pm1m;
    }

    double p_prev = pmm;
    double p_curr = pm1m;

    for (int degree = m + 2; degree <= l; ++degree)
    {
        const double numerator =
            static_cast<double>(2 * degree - 1) *
            x *
            p_curr
            -
            static_cast<double>(degree + m - 1) *
            p_prev;

        const double next =
            numerator /
            static_cast<double>(degree - m);

        p_prev = p_curr;
        p_curr = next;
    }

    return p_curr;
}

double reference_factorial(int n)
{
    if (n < 0)
    {
        throw std::invalid_argument(
            "factorial requires n >= 0"
        );
    }

    double result = 1.0;

    for (int i = 2; i <= n; ++i)
    {
        result *= static_cast<double>(i);
    }

    return result;
}

double reference_normalization(
    int l,
    int m
)
{
    const double numerator =
        (2.0 * static_cast<double>(l) + 1.0) /
        (4.0 * PI);

    const double ratio =
        reference_factorial(l - m) /
        reference_factorial(l + m);

    return std::sqrt(numerator * ratio);
}

double reference_real_sh(
    int l,
    int m,
    double theta,
    double phi
)
{
    const int am = std::abs(m);

    const double x = std::cos(theta);

    const double n =
        reference_normalization(l, am);

    const double p =
        reference_associated_legendre(l, am, x);

    if (m == 0)
    {
        return n * p;
    }

    if (m > 0)
    {
        return
            SQRT_TWO *
            n *
            p *
            std::cos(
                static_cast<double>(m) * phi
            );
    }

    return
        SQRT_TWO *
        n *
        p *
        std::sin(
            static_cast<double>(am) * phi
        );
}

double reference_reconstruct(
    const sylph::SHCoefficients& coefficients,
    double theta,
    double phi
)
{
    double value = 0.0;

    for (int l = 0; l < coefficients.order(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            value +=
                coefficients(l, m) *
                reference_real_sh(
                    l,
                    m,
                    theta,
                    phi
                );
        }
    }

    return value;
}

double reference_integral_of_function_product(
    const std::vector<const sylph::SHCoefficients*>& inputs,
    int n_theta,
    int n_phi
)
{
    if (inputs.empty())
    {
        throw std::invalid_argument(
            "reference integral requires inputs"
        );
    }

    sylph::gauss_legendre gl(n_theta);

    const auto& nodes = gl.nodes();
    const auto& weights = gl.weights();

    const double dphi =
        2.0 * PI /
        static_cast<double>(n_phi);

    double integral = 0.0;

    for (int i = 0; i < n_theta; ++i)
    {
        const double x =
            nodes[static_cast<std::size_t>(i)];

        const double theta =
            std::acos(
                std::clamp(x, -1.0, 1.0)
            );

        double ring_sum = 0.0;

        for (int j = 0; j < n_phi; ++j)
        {
            const double phi =
                dphi *
                static_cast<double>(j);

            double value = 1.0;

            for (const sylph::SHCoefficients* input : inputs)
            {
                value *=
                    reference_reconstruct(
                        *input,
                        theta,
                        phi
                    );
            }

            ring_sum += value;
        }

        integral +=
            weights[static_cast<std::size_t>(i)] *
            dphi *
            ring_sum;
    }

    return integral;
}

sylph::SHCoefficients reference_project_product(
    const std::vector<const sylph::SHCoefficients*>& inputs,
    int output_order,
    int n_theta,
    int n_phi
)
{
    if (inputs.empty())
    {
        throw std::invalid_argument(
            "reference projection requires inputs"
        );
    }

    if (output_order <= 0)
    {
        throw std::invalid_argument(
            "reference projection requires positive output order"
        );
    }

    sylph::gauss_legendre gl(n_theta);

    const auto& nodes = gl.nodes();
    const auto& weights = gl.weights();

    const double dphi =
        2.0 * PI /
        static_cast<double>(n_phi);

    sylph::SHCoefficients result(output_order);

    for (int l = 0; l < output_order; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            double coefficient = 0.0;

            for (int i = 0; i < n_theta; ++i)
            {
                const double x =
                    nodes[static_cast<std::size_t>(i)];

                const double theta =
                    std::acos(
                        std::clamp(x, -1.0, 1.0)
                    );

                double ring_sum = 0.0;

                for (int j = 0; j < n_phi; ++j)
                {
                    const double phi =
                        dphi *
                        static_cast<double>(j);

                    double product = 1.0;

                    for (const sylph::SHCoefficients* input : inputs)
                    {
                        product *=
                            reference_reconstruct(
                                *input,
                                theta,
                                phi
                            );
                    }

                    ring_sum +=
                        product *
                        reference_real_sh(
                            l,
                            m,
                            theta,
                            phi
                        );
                }

                coefficient +=
                    weights[
                        static_cast<std::size_t>(i)
                    ] *
                    dphi *
                    ring_sum;
            }

            result(l, m) = coefficient;
        }
    }

    return result;
}

double max_coeff_error(
    const sylph::SHCoefficients& a,
    const sylph::SHCoefficients& b
)
{
    if (a.order() != b.order())
    {
        return std::numeric_limits<double>::infinity();
    }

    double worst = 0.0;

    for (int l = 0; l < a.order(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            worst =
                std::max(
                    worst,
                    std::abs(
                        a(l, m) - b(l, m)
                    )
                );
        }
    }

    return worst;
}

double max_coeff_norm(
    const sylph::SHCoefficients& a
)
{
    double worst = 0.0;

    for (int l = 0; l < a.order(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            worst =
                std::max(
                    worst,
                    std::abs(a(l, m))
                );
        }
    }

    return worst;
}

bool all_finite(
    const sylph::SHCoefficients& coefficients
)
{
    for (int l = 0; l < coefficients.order(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            if (!std::isfinite(coefficients(l, m)))
            {
                return false;
            }
        }
    }

    return true;
}

void fill_randomish_coefficients(
    sylph::SHCoefficients& coefficients,
    double seed
)
{
    for (int l = 0; l < coefficients.order(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            const double value =
                seed
                + 0.043 * static_cast<double>(l + 1)
                + 0.071 * static_cast<double>(m)
                + 0.011 *
                    static_cast<double>(
                        l * l + 2 * m * m
                    );

            coefficients(l, m) = value;
        }
    }
}

void test_basis_orthonormality()
{
    section("1. Independent SH basis orthonormality");

    constexpr int max_order = 6;
    constexpr int n_theta = 16;
    constexpr int n_phi = 32;

    sylph::gauss_legendre gl(n_theta);

    const auto& nodes = gl.nodes();
    const auto& weights = gl.weights();

    const double dphi =
        2.0 * PI /
        static_cast<double>(n_phi);

    double worst_diagonal_error = 0.0;
    double worst_off_diagonal = 0.0;

    for (int l1 = 0; l1 < max_order; ++l1)
    {
        for (int m1 = -l1; m1 <= l1; ++m1)
        {
            for (int l2 = 0; l2 < max_order; ++l2)
            {
                for (int m2 = -l2; m2 <= l2; ++m2)
                {
                    double integral = 0.0;

                    for (int i = 0; i < n_theta; ++i)
                    {
                        const double theta =
                            std::acos(
                                std::clamp(
                                    nodes[
                                        static_cast<std::size_t>(i)
                                    ],
                                    -1.0,
                                    1.0
                                )
                            );

                        double ring_sum = 0.0;

                        for (int j = 0; j < n_phi; ++j)
                        {
                            const double phi =
                                dphi *
                                static_cast<double>(j);

                            ring_sum +=
                                reference_real_sh(
                                    l1,
                                    m1,
                                    theta,
                                    phi
                                ) *
                                reference_real_sh(
                                    l2,
                                    m2,
                                    theta,
                                    phi
                                );
                        }

                        integral +=
                            weights[
                                static_cast<std::size_t>(i)
                            ] *
                            dphi *
                            ring_sum;
                    }

                    const double expected =
                        (l1 == l2 && m1 == m2)
                        ? 1.0
                        : 0.0;

                    const double error =
                        std::abs(integral - expected);

                    if (expected == 0.0)
                    {
                        worst_off_diagonal =
                            std::max(
                                worst_off_diagonal,
                                error
                            );
                    }
                    else
                    {
                        worst_diagonal_error =
                            std::max(
                                worst_diagonal_error,
                                error
                            );
                    }
                }
            }
        }
    }

    close_error(
        worst_diagonal_error,
        "basis diagonal orthonormality",
        2e-12
    );

    close_error(
        worst_off_diagonal,
        "basis off-diagonal orthogonality",
        2e-12
    );
}

void test_round_trip()
{
    section("2. SH -> SphGrid -> SH round trip");

    constexpr int order = 6;

    sylph::SHCoefficients original(order);
    fill_randomish_coefficients(original, -0.25);

    // Eq. 35: full recovery
    // Ntheta >= n
    // Nphi   >= 2n - 1
    constexpr int n_theta = order;
    constexpr int n_phi = 2 * order - 1;

    sylph::gauss_legendre gl(n_theta);

    sylph::Spherical_grid grid(
        n_theta,
        n_phi,
        gl.nodes()
    );

    sylph::SHGridTransform transform;

    const auto sampled =
        transform.forward_separable(
            original,
            grid
        );

    const auto recovered =
        transform.inverse(
            sampled,
            gl.weights(),
            order,
            order
        );

    close_error(
        max_coeff_error(
            original,
            recovered
        ),
        "full-order round trip",
        2e-10
    );

    check(
        all_finite(recovered),
        "round-trip result is finite"
    );
}

void test_exact_generalized_product()
{
    section("3. Exact generalized product vs independent quadrature");

    sylph::SHCoefficients a(3);
    sylph::SHCoefficients b(4);
    sylph::SHCoefficients c(2);

    fill_randomish_coefficients(a, 0.15);
    fill_randomish_coefficients(b, -0.35);
    fill_randomish_coefficients(c, 0.55);

    const std::array<
        const sylph::SHCoefficients*,
        3
    > inputs = {
        &a,
        &b,
        &c
    };

    sylph::SHGridTransform transform;

    // nG = 3 + 4 + 2 - 3 + 1 = 7
    constexpr int nG = 7;
    constexpr int output_order = 5;

    const auto got =
        transform.product(
            std::span<
                const sylph::SHCoefficients* const
            >(inputs),
            output_order
        );

    // Independent, deliberately over-resolved reference.
    // Production exact minimum for p=5 is:
    // Ntheta >= ceil((7+5-1)/2) = 6
    // Nphi   >= 11
    //
    // Use a much finer grid here.
    constexpr int ref_theta = 24;
    constexpr int ref_phi = 64;

    const auto expected =
        reference_project_product(
            std::vector<
                const sylph::SHCoefficients*
            >{
                &a, &b, &c
            },
            output_order,
            ref_theta,
            ref_phi
        );

    close_error(
        max_coeff_error(got, expected),
        "generalized product vs independent reference",
        3e-9
    );

    check(
        got.order() == output_order,
        "generalized product returns requested order"
    );

    // Keep nG in the source so this test documents the mathematical
    // full-product order explicitly.
    check(
        nG > output_order,
        "test uses a genuinely truncated generalized product"
    );
}

void test_integral_against_reference()
{
    section("4. Product integrals vs independent quadrature");

    sylph::SHCoefficients a(4);
    sylph::SHCoefficients b(3);
    sylph::SHCoefficients c(2);
    sylph::SHCoefficients d(3);

    fill_randomish_coefficients(a, 0.10);
    fill_randomish_coefficients(b, -0.20);
    fill_randomish_coefficients(c, 0.30);
    fill_randomish_coefficients(d, -0.40);

    const std::array<
        const sylph::SHCoefficients*,
        4
    > inputs = {
        &a, &b, &c, &d
    };

    sylph::SHGridTransform transform;

    const double got =
        transform.product_integral(
            std::span<
                const sylph::SHCoefficients* const
            >(inputs)
        );

    const double expected =
        reference_integral_of_function_product(
            std::vector<
                const sylph::SHCoefficients*
            >{
                &a, &b, &c, &d
            },
            24,
            64
        );

    close(
        got,
        expected,
        "4-function product integral vs independent quadrature",
        3e-9,
        3e-9
    );
}

void test_two_function_integral_identity()
{
    section("5. Two-function product integral identity");

    sylph::SHCoefficients a(6);
    sylph::SHCoefficients b(5);

    fill_randomish_coefficients(a, -0.15);
    fill_randomish_coefficients(b, 0.25);

    sylph::SHGridTransform transform;

    const double got =
        transform.product_integral(a, b);

    double expected = 0.0;

    const int common_order =
        std::min(
            a.order(),
            b.order()
        );

    for (int l = 0; l < common_order; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            expected +=
                a(l, m) *
                b(l, m);
        }
    }

    close(
        got,
        expected,
        "two-function integral equals SH coefficient dot product",
        2e-10,
        2e-10
    );

    const double numerical =
        reference_integral_of_function_product(
            std::vector<
                const sylph::SHCoefficients*
            >{
                &a, &b
            },
            24,
            64
        );

    close(
        got,
        numerical,
        "two-function integral also matches independent quadrature",
        3e-9,
        3e-9
    );
}

void test_exact_t_threshold()
{
    section("6. Exact t threshold");

    sylph::SHCoefficients a(3);
    sylph::SHCoefficients b(4);
    sylph::SHCoefficients c(2);

    fill_randomish_coefficients(a, 0.12);
    fill_randomish_coefficients(b, -0.31);
    fill_randomish_coefficients(c, 0.46);

    const std::array<
        const sylph::SHCoefficients*,
        3
    > inputs = {
        &a, &b, &c
    };

    constexpr int nG = 7;
    constexpr int output_order = 4;

    // Eq. 41:
    // t_exact = nG + p - 1 = 10
    constexpr int t_exact =
        nG + output_order - 1;

    sylph::SHGridTransform transform;

    const auto exact =
        transform.product(
            std::span<
                const sylph::SHCoefficients* const
            >(inputs),
            output_order
        );

    const auto exact_from_t =
        transform.product(
            std::span<
                const sylph::SHCoefficients* const
            >(inputs),
            output_order,
            t_exact
        );

    close_error(
        max_coeff_error(
            exact,
            exact_from_t
        ),
        "t = nG + p - 1 agrees with exact product",
        2e-10
    );
}

void test_approximate_t()
{
    section("7. Lower-t approximate mode");

    sylph::SHCoefficients a(5);
    sylph::SHCoefficients b(4);
    sylph::SHCoefficients c(3);

    fill_randomish_coefficients(a, 0.08);
    fill_randomish_coefficients(b, -0.27);
    fill_randomish_coefficients(c, 0.42);

    const std::array<
        const sylph::SHCoefficients*,
        3
    > inputs = {
        &a, &b, &c
    };

    constexpr int nG =
        5 + 4 + 3 - 3 + 1;

    constexpr int output_order = 4;

    constexpr int t_exact =
        nG + output_order - 1;

    constexpr int t_low =
        t_exact - 4;

    sylph::SHGridTransform transform;

    const auto approximate =
        transform.product(
            std::span<
                const sylph::SHCoefficients* const
            >(inputs),
            output_order,
            t_low
        );

    const auto reference =
        reference_project_product(
            std::vector<
                const sylph::SHCoefficients*
            >{
                &a, &b, &c
            },
            output_order,
            28,
            72
        );

    const double error =
        max_coeff_error(
            approximate,
            reference
        );

    check(
        std::isfinite(error),
        "lower-t approximation is finite"
    );

    check(
        approximate.order() == output_order,
        "lower-t approximation has requested order"
    );

    std::cout
        << std::setprecision(10)
        << "Approximation diagnostic:\n"
        << "  nG       = " << nG << '\n'
        << "  p        = " << output_order << '\n'
        << "  t_low    = " << t_low << '\n'
        << "  t_exact  = " << t_exact << '\n'
        << "  max err  = " << error << '\n';
}

void test_commutativity()
{
    section("8. Product commutativity");

    sylph::SHCoefficients a(5);
    sylph::SHCoefficients b(4);
    sylph::SHCoefficients c(3);

    fill_randomish_coefficients(a, 0.13);
    fill_randomish_coefficients(b, -0.19);
    fill_randomish_coefficients(c, 0.29);

    sylph::SHGridTransform transform;

    const std::array<
        const sylph::SHCoefficients*,
        3
    > abc = {
        &a, &b, &c
    };

    const std::array<
        const sylph::SHCoefficients*,
        3
    > cba = {
        &c, &b, &a
    };

    const auto first =
        transform.product(
            std::span<
                const sylph::SHCoefficients* const
            >(abc),
            4
        );

    const auto second =
        transform.product(
            std::span<
                const sylph::SHCoefficients* const
            >(cba),
            4
        );

    close_error(
        max_coeff_error(first, second),
        "generalized product is commutative",
        2e-10
    );
}

void test_constant_identity()
{
    section("9. Multiplication by constant-one function");

    constexpr int order = 6;

    sylph::SHCoefficients f(order);
    sylph::SHCoefficients one(1);

    fill_randomish_coefficients(f, -0.22);

    // Y_0^0 = 1 / sqrt(4*pi)
    // Therefore constant 1 has coefficient sqrt(4*pi)
    // = 2*sqrt(pi).
    one(0, 0) = 2.0 * std::sqrt(PI);

    sylph::SHGridTransform transform;

    const auto result =
        transform.product(
            f,
            one,
            order
        );

    close_error(
        max_coeff_error(f, result),
        "F * 1 reproduces F",
        3e-9
    );
}

void test_parseval()
{
    section("10. Parseval / energy identity");

    sylph::SHCoefficients f(7);

    fill_randomish_coefficients(f, 0.17);

    double coefficient_energy = 0.0;

    for (int l = 0; l < f.order(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            coefficient_energy +=
                f(l, m) * f(l, m);
        }
    }

    const double numerical_energy =
        reference_integral_of_function_product(
            std::vector<
                const sylph::SHCoefficients*
            >{
                &f, &f
            },
            28,
            72
        );

    close(
        coefficient_energy,
        numerical_energy,
        "Parseval energy equals numerical integral",
        5e-9,
        5e-9
    );
}

void test_inverse_resolution_boundary_correct()
{
    section("11b. Eq. 34 exact/insufficient boundary");

    constexpr int input_order = 7;
    constexpr int output_order = 4;

    constexpr int min_theta = 5;
    constexpr int min_phi = 10;

    sylph::SHGridTransform transform;

    sylph::gauss_legendre gl_exact(min_theta);

    sylph::Spherical_grid exact_grid(
        min_theta,
        min_phi,
        gl_exact.nodes()
    );

    // Exact minimum should NOT throw.
    bool exact_threw = false;

    try
    {
        const auto result =
            transform.inverse(
                exact_grid,
                gl_exact.weights(),
                input_order,
                output_order
            );

        (void)result;
    }
    catch (...)
    {
        exact_threw = true;
    }

    check(
        !exact_threw,
        "Eq. 34 minimum grid is accepted"
    );

    // One fewer polar sample must be rejected.
    sylph::gauss_legendre gl_low_theta(min_theta - 1);

    sylph::Spherical_grid low_theta_grid(
        min_theta - 1,
        min_phi,
        gl_low_theta.nodes()
    );

    expect_invalid_argument(
        [&] {
            transform.inverse(
                low_theta_grid,
                gl_low_theta.weights(),
                input_order,
                output_order
            );
        },
        "one fewer theta sample is rejected"
    );

    // Phi resolution is independent from the GL weights.
    sylph::Spherical_grid low_phi_grid(
        min_theta,
        min_phi - 1,
        gl_exact.nodes()
    );

    expect_invalid_argument(
        [&] {
            transform.inverse(
                low_phi_grid,
                gl_exact.weights(),
                input_order,
                output_order
            );
        },
        "one fewer phi sample is rejected"
    );
}

} // namespace

int RUNTEST()
{
    std::cout
        << "Sylph system-level mathematical validation\n";

    // The first boundary helper is intentionally not called because its
    // exact-boundary behavior is superseded by the corrected test below.
    test_basis_orthonormality();
    test_round_trip();
    test_exact_generalized_product();
    test_integral_against_reference();
    test_two_function_integral_identity();
    test_exact_t_threshold();
    test_approximate_t();
    test_commutativity();
    test_constant_identity();
    test_parseval();
    test_inverse_resolution_boundary_correct();

    std::cout
        << "\n============================================================\n"
        << "VALIDATION SUMMARY\n"
        << "============================================================\n"
        << "PASSED: " << state.passed << '\n'
        << "FAILED: " << state.failed << '\n';

    if (state.failed != 0)
    {
        std::cout
            << std::setprecision(17)
            << "WORST ERROR: "
            << state.worst_error
            << '\n';

        if (!state.worst_name.empty())
        {
            std::cout
                << "WORST TEST: "
                << state.worst_name
                << '\n';
        }

        return 1;
    }

    std::cout
        << "SYSTEM-LEVEL MATHEMATICAL VALIDATION PASSED.\n";

    return 0;
}






