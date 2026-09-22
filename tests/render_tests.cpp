#include "render_tests.hpp"

#include "render/brdf.hpp"
#include "render/glossy_relighting.hpp"
#include "render/lighting.hpp"
#include "render/oof.hpp"
#include "render/shadow_field.hpp"
#include "render/visibility.hpp"

#include "types/constants.hpp"

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <iostream>
#include <span>
#include <vector>
bool RenderTests::check_close(
    double got,
    double expected,
    double tolerance,
    const char* name
)
{
    const double error =
        std::abs(got - expected);

    const bool passed =
        error <= tolerance;

    std::cout
        << (passed ? "[PASS] " : "[FAIL] ")
        << name
        << " | got=" << got
        << " expected=" << expected
        << " error=" << error
        << '\n';

    return passed;
}

bool RenderTests::check_true(
    bool condition,
    const char* name
)
{
    std::cout
        << (condition ? "[PASS] " : "[FAIL] ")
        << name
        << '\n';

    return condition;
}


bool RenderTests::test_brdf()
{
    std::cout
        << "\n=== BRDF ===\n";

    bool passed = true;

    constexpr int order = 4;

    render::BRDF brdf(order);

    passed &= check_true(
        brdf.order() == order,
        "BRDF order"
    );

    brdf.coefficients()(2, -1) =
        3.14159;

    passed &= check_close(
        brdf.coefficients()(2, -1),
        3.14159,
        1e-12,
        "BRDF coefficient read/write"
    );

    brdf.clear();

    passed &= check_close(
        brdf.coefficients()(2, -1),
        0.0,
        1e-12,
        "BRDF clear"
    );

    return passed;
}


bool RenderTests::test_visibility()
{
    std::cout
        << "\n=== VISIBILITY ===\n";

    bool passed = true;

    constexpr int order = 5;

    render::Visibility visibility(order);

    passed &= check_true(
        visibility.order() == order,
        "Visibility order"
    );

    visibility.coefficients()(3, 2) =
        -2.75;

    passed &= check_close(
        visibility.coefficients()(3, 2),
        -2.75,
        1e-12,
        "Visibility coefficient read/write"
    );

    visibility.clear();

    passed &= check_close(
        visibility.coefficients()(3, 2),
        0.0,
        1e-12,
        "Visibility clear"
    );

    return passed;
}


bool RenderTests::test_oof()
{
    std::cout
        << "\n=== OOF ===\n";

    bool passed = true;

    constexpr int order = 3;

    const Eigen::Vector3f min(
        -1.0f,
        -2.0f,
        -3.0f
    );

    const Eigen::Vector3f max(
         1.0f,
         2.0f,
         3.0f
    );

    render::OOF oof(
        4,
        5,
        6,
        min,
        max,
        order
    );

    passed &= check_true(
        oof.resolution_x() == 4,
        "OOF X resolution"
    );

    passed &= check_true(
        oof.resolution_y() == 5,
        "OOF Y resolution"
    );

    passed &= check_true(
        oof.resolution_z() == 6,
        "OOF Z resolution"
    );

    passed &= check_true(
        oof.visibility_order() == order,
        "OOF visibility order"
    );

    oof.at(0, 0, 0)
        .coefficients()(0, 0) = 7.0;

    passed &= check_close(
        oof.at(0, 0, 0)
            .coefficients()(0, 0),
        7.0,
        1e-12,
        "OOF direct cell access"
    );

    // Exact minimum corner.
    const auto& min_visibility =
        oof.lookup(min);

    passed &= check_close(
        min_visibility.coefficients()(0, 0),
        7.0,
        1e-12,
        "OOF lookup at minimum bound"
    );

    // Out-of-bounds positions are clamped
    // into the volume by lookup().
    const auto& outside_visibility =
        oof.lookup(
            Eigen::Vector3f(
                -100.0f,
                -100.0f,
                -100.0f
            )
        );

    passed &= check_close(
        outside_visibility.coefficients()(0, 0),
        7.0,
        1e-12,
        "OOF clamped lookup"
    );

    return passed;
}


bool RenderTests::test_glossy_relighting()
{
    std::cout
        << "\n=== GLOSSY RELIGHTING ===\n";

    bool passed = true;

    /*
        Test constant spherical functions:

            L(ω) = 1
            B(ω) = 1
            V(ω) = 1

        Therefore:

            ∫ L B V dΩ = 4π

        The SH coefficient for constant 1 is:

            f00 = 2 * sqrt(pi)
    */

    constexpr int lighting_order = 2;
    constexpr int brdf_order = 3;
    constexpr int visibility_order = 4;

    const double constant_coefficient =
        2.0 * std::sqrt(sylph::pi);

    render::Lighting lighting(
        lighting_order
    );

    render::BRDF brdf(
        brdf_order
    );

    render::Visibility visibility(
        visibility_order
    );

    for (int channel = 0; channel < 3; ++channel) {

        lighting.coefficients(
            static_cast<render::LightChannel>(channel)
        )(0, 0) =
            constant_coefficient;
    }

    brdf.coefficients()(0, 0) =
        constant_coefficient;

    visibility.coefficients()(0, 0) =
        constant_coefficient;

    render::GlossyRelighting relighting;

    const Eigen::Vector3f result =
        relighting.evaluate(
            lighting,
            brdf,
            visibility
        );

    const double expected =
        4.0 * sylph::pi;

    passed &= check_close(
        result.x(),
        expected,
        1e-6,
        "Glossy relighting R"
    );

    passed &= check_close(
        result.y(),
        expected,
        1e-6,
        "Glossy relighting G"
    );

    passed &= check_close(
        result.z(),
        expected,
        1e-6,
        "Glossy relighting B"
    );

    return passed;
}


bool RenderTests::test_shadow_field()
{
    std::cout
        << "\n=== SHADOW FIELD ===\n";

    bool passed = true;

    /*
        Again use constant functions:

            L      = 1
            Vself  = 1
            Voof   = 1

        Therefore:

            ∫ L Vself Voof dΩ = 4π
    */

    constexpr int lighting_order = 2;
    constexpr int self_order = 3;
    constexpr int oof_order = 4;

    const double constant_coefficient =
        2.0 * std::sqrt(sylph::pi);

    render::Lighting lighting(
        lighting_order
    );

    render::Visibility self_visibility(
        self_order
    );

    const Eigen::Vector3f min(
        -1.0f,
        -1.0f,
        -1.0f
    );

    const Eigen::Vector3f max(
         1.0f,
         1.0f,
         1.0f
    );

    render::OOF oof(
        4,
        4,
        4,
        min,
        max,
        oof_order
    );

    for (int channel = 0; channel < 3; ++channel) {

        lighting.coefficients(
            static_cast<render::LightChannel>(channel)
        )(0, 0) =
            constant_coefficient;
    }

    self_visibility.coefficients()(0, 0) =
        constant_coefficient;

    // Make every OOF cell represent visibility = 1.
    for (int z = 0; z < oof.resolution_z(); ++z) {
        for (int y = 0; y < oof.resolution_y(); ++y) {
            for (int x = 0; x < oof.resolution_x(); ++x) {

                oof.at(x, y, z)
                    .coefficients()(0, 0) =
                    constant_coefficient;
            }
        }
    }

    const std::array<
        const render::OOF*,
        1
    > oofs = {
        &oof
    };

    render::ShadowField shadow_field;

    const Eigen::Vector3f position(
        0.0f,
        0.0f,
        0.0f
    );

    const Eigen::Vector3f result =
        shadow_field.evaluate(
            lighting,
            self_visibility,
            position,
            std::span<const render::OOF* const>(
                oofs
            )
        );

    const double expected =
        4.0 * sylph::pi;

    passed &= check_close(
        result.x(),
        expected,
        1e-6,
        "Shadow field R"
    );

    passed &= check_close(
        result.y(),
        expected,
        1e-6,
        "Shadow field G"
    );

    passed &= check_close(
        result.z(),
        expected,
        1e-6,
        "Shadow field B"
    );

    return passed;
}


bool RenderTests::run_all()
{
    std::cout
        << "\n"
        << "========================================\n"
        << "        SYLPH RENDER TESTS\n"
        << "========================================\n";

    bool passed = true;

    passed &= test_brdf();
    passed &= test_visibility();
    passed &= test_oof();
    passed &= test_glossy_relighting();
    passed &= test_shadow_field();

    std::cout
        << "\n========================================\n"
        << (passed
                ? "ALL RENDER TESTS PASSED"
                : "RENDER TESTS FAILED")
        << "\n========================================\n";

    return passed;
}