#include "glossy/glossy_example.hpp"
#include "sh_lobe/sh_lobe_example.hpp"
#include "shadow/shadow_example.hpp"
#include "oof_demo/oof_example.hpp"
#include "sh_basis/sh_basis_example.hpp"
#include "product_demo/product_demo_example.hpp"

#include "render/renderer.hpp"

#include <Eigen/Dense>
#include <SDL3/SDL.h>
#include <glad/gl.h>

#include <cstdlib>
#include <iostream>

namespace {
constexpr int WINDOW_WIDTH  = 1280;
constexpr int WINDOW_HEIGHT = 720;
}

int main()
{
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::cerr << "Failed to initialize SDL: " << SDL_GetError() << "\n";
        return EXIT_FAILURE;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    SDL_Window* window = SDL_CreateWindow(
        "Sylph Demo - SH Basis / Lobe / Glossy / Shadow PRT / OOF",
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE
    );

    if (!window) {
        std::cerr << "Failed to create window: " << SDL_GetError() << "\n";
        SDL_Quit();
        return EXIT_FAILURE;
    }

    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    if (!gl_context) {
        std::cerr << "Failed to create GL context: " << SDL_GetError() << "\n";
        SDL_DestroyWindow(window);
        SDL_Quit();
        return EXIT_FAILURE;
    }

    if (!gladLoadGL((GLADloadfunc)SDL_GL_GetProcAddress)) {
        std::cerr << "Failed to initialize GLAD\n";
        SDL_GL_DestroyContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return EXIT_FAILURE;
    }

    if (!SDL_GL_SetSwapInterval(1)) {
        std::cerr << "Warning: Unable to set VSync: " << SDL_GetError() << "\n";
    }

    {
        render::Renderer renderer;

        // ==========================================================
        // DEMO SELECTOR
        // Change the class name below to run a different demo!
        // 1. examples::SHBasisExample
        // 2. examples::SHLobeExample
        // 3. examples::ProductDemoExample
        // 4. examples::GlossyExample
        // 5. examples::ShadowExample
        // 6. examples::OofExample
        // ==========================================================
        examples::GlossyExample example(WINDOW_WIDTH, WINDOW_HEIGHT);

        bool running = true;
        Uint64 last_time = SDL_GetPerformanceCounter();

        while (running) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_EVENT_QUIT) {
                    running = false;
                }
                if (event.type == SDL_EVENT_WINDOW_RESIZED) {
                    example.resize(event.window.data1, event.window.data2);
                }
                example.process_event(event, window);
            }

            Uint64 current_time = SDL_GetPerformanceCounter();
            float dt = static_cast<float>(current_time - last_time) / static_cast<float>(SDL_GetPerformanceFrequency());
            last_time = current_time;

            example.update(dt);
            example.render(renderer);

            SDL_GL_SwapWindow(window);
        }
    }

    SDL_GL_DestroyContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return EXIT_SUCCESS;
}



