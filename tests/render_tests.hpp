#pragma once

class RenderTests {
public:
    static bool run_all();

private:
    static bool test_brdf();
    static bool test_visibility();
    static bool test_oof();
    static bool test_glossy_relighting();
    static bool test_shadow_field();

    static bool check_close(
        double got,
        double expected,
        double tolerance,
        const char* name
    );

    static bool check_true(
        bool condition,
        const char* name
    );
};