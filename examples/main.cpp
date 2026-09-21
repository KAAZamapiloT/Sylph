#include "sylph/spherical_grid.hpp"
#include "sylph/gauss_legendre.hpp"
#include "types/constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = sylph::pi;

struct TestState
{
    int passed = 0;
    int failed = 0;

    double worst_error = 0.0;
    std::string worst_test;
};

TestState state;

void pass()
{
    ++state.passed;
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
        state.worst_test = name;
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
        return;
    }

    ++state.failed;
    std::cerr << "[FAIL] " << name << '\n';
}

void close(
    double got,
    double expected,
    const std::string& name,
    double abs_tol = 1e-12,
    double rel_tol = 1e-10
)
{
    const double error =
        std::abs(got - expected);

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
        return;
    }

    fail(
        name,
        got,
        expected,
        error,
        tolerance
    );
}

void expect_invalid_argument(
    auto&& function,
    const std::string& name
)
{
    try
    {
        function();

        ++state.failed;
        std::cerr
            << "[FAIL] " << name
            << "\n       expected std::invalid_argument"
               " but no exception was thrown\n";
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

void section(const std::string& name)
{
    std::cout
        << "\n============================================================\n"
        << name << '\n'
        << "============================================================\n";
}

sylph::Spherical_grid make_grid(
    int nTheta,
    int nPhi
)
{
    sylph::gauss_legendre gl(nTheta);

    return sylph::Spherical_grid(
        nTheta,
        nPhi,
        gl.nodes()
    );
}

void fill_formula(
    sylph::Spherical_grid& grid,
    double a,
    double b
)
{
    for (int i = 0;
         i < grid.theta_size();
         ++i)
    {
        for (int j = 0;
             j < grid.phi_size();
             ++j)
        {
            grid(i, j) =
                a * static_cast<double>(i + 1) +
                b * static_cast<double>(j + 1);
        }
    }
}

void fill_random(
    sylph::Spherical_grid& grid,
    std::mt19937_64& rng,
    double scale = 1.0
)
{
    std::uniform_real_distribution<double> dist(
        -scale,
        scale
    );

    for (int i = 0;
         i < grid.theta_size();
         ++i)
    {
        for (int j = 0;
             j < grid.phi_size();
             ++j)
        {
            grid(i, j) = dist(rng);
        }
    }
}

double max_grid_error(
    const sylph::Spherical_grid& got,
    const sylph::Spherical_grid& expected,
    int& worst_i,
    int& worst_j
)
{
    check(
        got.theta_size() == expected.theta_size() &&
        got.phi_size() == expected.phi_size(),
        "max_grid_error dimension agreement"
    );

    double maximum = 0.0;
    worst_i = -1;
    worst_j = -1;

    for (int i = 0;
         i < got.theta_size();
         ++i)
    {
        for (int j = 0;
             j < got.phi_size();
             ++j)
        {
            const double error =
                std::abs(
                    got(i, j) -
                    expected(i, j)
                );

            if (error > maximum)
            {
                maximum = error;
                worst_i = i;
                worst_j = j;
            }
        }
    }

    return maximum;
}

// ============================================================
// 1. Deterministic pointwise multiplication
//
// Verifies:
//     result(i,j) = A(i,j) * B(i,j)
// exactly at every grid point.
//
// This is the fundamental operation described by the paper's
// SphGrid product stage.
// ============================================================

void test_basic_pointwise_multiplication()
{
    section("1. BASIC POINTWISE MULTIPLICATION");

    for (int nTheta : {1, 2, 3, 5, 8})
    {
        for (int nPhi : {1, 2, 4, 7, 9})
        {
            auto a = make_grid(nTheta, nPhi);
            auto b = make_grid(nTheta, nPhi);

            fill_formula(
                a,
                1.25,
                -0.37
            );

            fill_formula(
                b,
                -0.81,
                0.23
            );

            const auto a_before = a;

            a.multiply_inplace(b);

            check(
                a.theta_size() == nTheta,
                "basic multiply preserves theta size"
            );

            check(
                a.phi_size() == nPhi,
                "basic multiply preserves phi size"
            );

            for (int i = 0;
                 i < nTheta;
                 ++i)
            {
                for (int j = 0;
                     j < nPhi;
                     ++j)
                {
                    const double expected =
                        a_before(i, j) *
                        b(i, j);

                    close(
                        a(i, j),
                        expected,
                        "pointwise product"
                        " theta=" +
                        std::to_string(i) +
                        " phi=" +
                        std::to_string(j),
                        1e-13,
                        1e-13
                    );
                }
            }
        }
    }
}

// ============================================================
// 2. Zero annihilator
//
//     A * 0 = 0
// ============================================================

void test_zero_annihilator()
{
    section("2. ZERO ANNIHILATOR");

    for (int nTheta : {1, 3, 5, 8})
    {
        for (int nPhi : {1, 4, 9})
        {
            auto a = make_grid(nTheta, nPhi);
            auto zero = make_grid(nTheta, nPhi);

            fill_formula(
                a,
                2.5,
                -1.7
            );

            for (int i = 0;
                 i < nTheta;
                 ++i)
            {
                for (int j = 0;
                     j < nPhi;
                     ++j)
                {
                    zero(i, j) = 0.0;
                }
            }

            a.multiply_inplace(zero);

            for (int i = 0;
                 i < nTheta;
                 ++i)
            {
                for (int j = 0;
                     j < nPhi;
                     ++j)
                {
                    close(
                        a(i, j),
                        0.0,
                        "A * 0"
                        " theta=" +
                        std::to_string(i) +
                        " phi=" +
                        std::to_string(j),
                        1e-13,
                        1e-13
                    );
                }
            }
        }
    }
}

// ============================================================
// 3. Identity
//
//     A * 1 = A
// ============================================================

void test_identity()
{
    section("3. MULTIPLICATIVE IDENTITY");

    for (int nTheta : {1, 2, 5, 8})
    {
        for (int nPhi : {1, 3, 7, 9})
        {
            auto a = make_grid(nTheta, nPhi);
            auto one = make_grid(nTheta, nPhi);

            fill_formula(
                a,
                -0.91,
                2.17
            );

            for (int i = 0;
                 i < nTheta;
                 ++i)
            {
                for (int j = 0;
                     j < nPhi;
                     ++j)
                {
                    one(i, j) = 1.0;
                }
            }

            const auto original = a;

            a.multiply_inplace(one);

            for (int i = 0;
                 i < nTheta;
                 ++i)
            {
                for (int j = 0;
                     j < nPhi;
                     ++j)
                {
                    close(
                        a(i, j),
                        original(i, j),
                        "A * 1 = A"
                        " theta=" +
                        std::to_string(i) +
                        " phi=" +
                        std::to_string(j),
                        1e-13,
                        1e-13
                    );
                }
            }
        }
    }
}

// ============================================================
// 4. Commutativity
//
//     A * B = B * A
//
// Multiplication is done pointwise, so this should hold exactly
// up to floating-point evaluation order.
// ============================================================

void test_commutativity()
{
    section("4. COMMUTATIVITY");

    std::mt19937_64 rng(
        0xA11CE55E12345678ULL
    );

    for (int trial = 0;
         trial < 100;
         ++trial)
    {
        const int nTheta =
            1 + static_cast<int>(rng() % 8);

        const int nPhi =
            1 + static_cast<int>(rng() % 10);

        auto a = make_grid(nTheta, nPhi);
        auto b = make_grid(nTheta, nPhi);

        fill_random(a, rng, 10.0);
        fill_random(b, rng, 10.0);

        auto ab = a;
        auto ba = b;

        ab.multiply_inplace(b);
        ba.multiply_inplace(a);

        int worst_i = -1;
        int worst_j = -1;

        const double error =
            max_grid_error(
                ab,
                ba,
                worst_i,
                worst_j
            );

        close(
            error,
            0.0,
            "A * B == B * A"
            " trial=" +
            std::to_string(trial),
            1e-14,
            1e-14
        );
    }
}

// ============================================================
// 5. Associativity
//
//     (A * B) * C = A * (B * C)
//
// Again, pointwise floating-point multiplication may differ by a
// few ulps, so this uses a numerical tolerance rather than exact
// equality.
// ============================================================

void test_associativity()
{
    section("5. ASSOCIATIVITY");

    std::mt19937_64 rng(
        0x51A551C8A55ULL
    );

    for (int trial = 0;
         trial < 100;
         ++trial)
    {
        const int nTheta =
            1 + static_cast<int>(rng() % 8);

        const int nPhi =
            1 + static_cast<int>(rng() % 10);

        auto a = make_grid(nTheta, nPhi);
        auto b = make_grid(nTheta, nPhi);
        auto c = make_grid(nTheta, nPhi);

        fill_random(a, rng, 2.0);
        fill_random(b, rng, 2.0);
        fill_random(c, rng, 2.0);

        auto left = a;
        left.multiply_inplace(b);
        left.multiply_inplace(c);

        auto right = b;
        right.multiply_inplace(c);

        right = a;

        // Rebuild the right side explicitly:
        auto bc = b;
        bc.multiply_inplace(c);

        right = a;
        right.multiply_inplace(bc);

        int worst_i = -1;
        int worst_j = -1;

        const double error =
            max_grid_error(
                left,
                right,
                worst_i,
                worst_j
            );

        close(
            error,
            0.0,
            "(A * B) * C == A * (B * C)"
            " trial=" +
            std::to_string(trial),
            1e-13,
            1e-13
        );
    }
}

// ============================================================
// 6. Self multiplication
//
//     A * A = A^2
//
// Explicitly tests aliasing of the lhs/rhs when the same object is
// passed to multiply_inplace().
// ============================================================

void test_self_multiplication()
{
    section("6. SELF MULTIPLICATION / ALIASING");

    std::mt19937_64 rng(
        0x5E1FA11ULL
    );

    for (int trial = 0;
         trial < 50;
         ++trial)
    {
        const int nTheta =
            1 + static_cast<int>(rng() % 7);

        const int nPhi =
            1 + static_cast<int>(rng() % 9);

        auto a = make_grid(nTheta, nPhi);
        fill_random(a, rng, 3.0);

        const auto original = a;

        a.multiply_inplace(a);

        for (int i = 0;
             i < nTheta;
             ++i)
        {
            for (int j = 0;
                 j < nPhi;
                 ++j)
            {
                const double expected =
                    original(i, j) *
                    original(i, j);

                close(
                    a(i, j),
                    expected,
                    "A * A"
                    " trial=" +
                    std::to_string(trial) +
                    " theta=" +
                    std::to_string(i) +
                    " phi=" +
                    std::to_string(j),
                    1e-13,
                    1e-13
                );
            }
        }
    }
}

// ============================================================
// 7. Random reference comparison
//
// Independent reference:
//
//     expected(i,j) = A(i,j) * B(i,j)
//
// This repeats the actual contract over many random grids.
// ============================================================

void test_random_reference()
{
    section("7. RANDOM REFERENCE COMPARISON");

    std::mt19937_64 rng(
        0xC0FFEE123456789ULL
    );

    double worst_error = 0.0;
    int worst_trial = -1;
    int worst_i = -1;
    int worst_j = -1;

    for (int trial = 0;
         trial < 200;
         ++trial)
    {
        const int nTheta =
            1 + static_cast<int>(rng() % 10);

        const int nPhi =
            1 + static_cast<int>(rng() % 12);

        auto a = make_grid(nTheta, nPhi);
        auto b = make_grid(nTheta, nPhi);
        auto expected = make_grid(nTheta, nPhi);

        fill_random(a, rng, 100.0);
        fill_random(b, rng, 100.0);

        for (int i = 0;
             i < nTheta;
             ++i)
        {
            for (int j = 0;
                 j < nPhi;
                 ++j)
            {
                expected(i, j) =
                    a(i, j) * b(i, j);
            }
        }

        a.multiply_inplace(b);

        int local_i = -1;
        int local_j = -1;

        const double error =
            max_grid_error(
                a,
                expected,
                local_i,
                local_j
            );

        if (error > worst_error)
        {
            worst_error = error;
            worst_trial = trial;
            worst_i = local_i;
            worst_j = local_j;
        }

        close(
            error,
            0.0,
            "random pointwise reference"
            " trial=" +
            std::to_string(trial),
            1e-12,
            1e-12
        );
    }

    std::cout
        << std::setprecision(17)
        << "\n[TRACE] random worst error = "
        << worst_error
        << " at trial="
        << worst_trial
        << " (" << worst_i
        << "," << worst_j
        << ")\n";
}

// ============================================================
// 8. Dimension mismatch
//
// Different SphGrid resolutions must NOT silently multiply.
//
// This prevents:
//     A(i,j) *= B(i,j)
//
// from becoming a partial or invalid operation.
// ============================================================

void test_dimension_mismatch()
{
    section("8. DIMENSION MISMATCH VALIDATION");

    {
        auto a = make_grid(4, 7);
        auto b = make_grid(5, 7);

        expect_invalid_argument(
            [&]
            {
                a.multiply_inplace(b);
            },
            "reject mismatched theta dimensions"
        );
    }

    {
        auto a = make_grid(4, 7);
        auto b = make_grid(4, 8);

        expect_invalid_argument(
            [&]
            {
                a.multiply_inplace(b);
            },
            "reject mismatched phi dimensions"
        );
    }

    {
        auto a = make_grid(1, 1);
        auto b = make_grid(2, 2);

        expect_invalid_argument(
            [&]
            {
                a.multiply_inplace(b);
            },
            "reject completely mismatched grids"
        );
    }
}

// ============================================================
// 9. Grid coordinate preservation
//
// Multiplication should modify ONLY values.
//
// theta(i) and phi(j) must remain unchanged.
// ============================================================

void test_coordinates_unchanged()
{
    section("9. COORDINATES UNCHANGED");

    constexpr int nTheta = 6;
    constexpr int nPhi = 11;

    auto a = make_grid(nTheta, nPhi);
    auto b = make_grid(nTheta, nPhi);

    fill_formula(a, 1.0, 2.0);
    fill_formula(b, -0.5, 0.25);

    std::vector<double> theta_before(
        static_cast<std::size_t>(nTheta)
    );

    std::vector<double> phi_before(
        static_cast<std::size_t>(nPhi)
    );

    for (int i = 0;
         i < nTheta;
         ++i)
    {
        theta_before[
            static_cast<std::size_t>(i)
        ] = a.theta(i);
    }

    for (int j = 0;
         j < nPhi;
         ++j)
    {
        phi_before[
            static_cast<std::size_t>(j)
        ] = a.phi(j);
    }

    a.multiply_inplace(b);

    for (int i = 0;
         i < nTheta;
         ++i)
    {
        close(
            a.theta(i),
            theta_before[
                static_cast<std::size_t>(i)
            ],
            "theta coordinate unchanged"
            " i=" +
            std::to_string(i),
            1e-15,
            1e-15
        );
    }

    for (int j = 0;
         j < nPhi;
         ++j)
    {
        close(
            a.phi(j),
            phi_before[
                static_cast<std::size_t>(j)
            ],
            "phi coordinate unchanged"
            " j=" +
            std::to_string(j),
            1e-15,
            1e-15
        );
    }
}

// ============================================================
// 10. Special floating-point values that remain finite
//
// This is not a mathematical identity test; it ensures normal
// finite inputs do not unexpectedly become NaN/Inf.
// ============================================================

void test_finite_values()
{
    section("10. FINITE VALUE PRESERVATION");

    auto a = make_grid(5, 9);
    auto b = make_grid(5, 9);

    for (int i = 0;
         i < a.theta_size();
         ++i)
    {
        for (int j = 0;
             j < a.phi_size();
             ++j)
        {
            a(i, j) =
                0.25 *
                static_cast<double>(i + 1);

            b(i, j) =
                -0.5 *
                static_cast<double>(j + 1);
        }
    }

    a.multiply_inplace(b);

    for (int i = 0;
         i < a.theta_size();
         ++i)
    {
        for (int j = 0;
             j < a.phi_size();
             ++j)
        {
            check(
                std::isfinite(a(i, j)),
                "product remains finite"
                " theta=" +
                std::to_string(i) +
                " phi=" +
                std::to_string(j)
            );
        }
    }
}

} // namespace

int main()
{
    test_basic_pointwise_multiplication();
    test_zero_annihilator();
    test_identity();
    test_commutativity();
    test_associativity();
    test_self_multiplication();
    test_random_reference();
    test_dimension_mismatch();
    test_coordinates_unchanged();
    test_finite_values();

    std::cerr
        << "\n\n============================================================\n"
        << "            SPHERICAL GRID MULTIPLICATION SUMMARY\n"
        << "============================================================\n"
        << "[INFO] Successful tests are intentionally silent.\n"
        << "[INFO] Only failures are printed above.\n"
        << "------------------------------------------------------------\n"
        << "TOTAL PASSED : " << state.passed << '\n'
        << "TOTAL FAILED : " << state.failed << '\n'
        << "------------------------------------------------------------\n";

    if (state.failed == 0)
    {
        std::cerr
            << "RESULT: ALL MULTIPLICATION TESTS PASSED\n";
        return 0;
    }

    std::cerr
        << "RESULT: MULTIPLICATION FAILURES DETECTED\n"
        << "WORST ERROR : "
        << std::setprecision(17)
        << state.worst_error
        << '\n'
        << "WORST TEST  : "
        << state.worst_test
        << '\n';

    return 1;
}
