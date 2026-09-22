#include "glossy/glossy_example.hpp"
#include "gallery/gallery_example.hpp"
#include "gallery/gallery_example.hpp"
#include "house/house_example.hpp"
#include "oof_demo/oof_example.hpp"

#include "render/renderer.hpp"

#include <Eigen/Dense>

#include <SDL3/SDL.h>
#include <glad/gl.h>

#include <cstdlib>
#include <iostream>

namespace
{

constexpr int WINDOW_WIDTH  = 1280;
constexpr int WINDOW_HEIGHT = 720;

}

int main()
{
    // --------------------------------------------------------
    // SDL
    // --------------------------------------------------------

    if (SDL_Init(SDL_INIT_VIDEO) == 0)
    {
        std::cerr << "Failed to initialize SDL: "
                  << SDL_GetError()
                  << '\n';
        return EXIT_FAILURE;
    }

    // --------------------------------------------------------
    // OpenGL context
    // --------------------------------------------------------

    SDL_GL_SetAttribute(
        SDL_GL_CONTEXT_MAJOR_VERSION,
        3
    );

    SDL_GL_SetAttribute(
        SDL_GL_CONTEXT_MINOR_VERSION,
        3
    );

    SDL_GL_SetAttribute(
        SDL_GL_CONTEXT_PROFILE_MASK,
        SDL_GL_CONTEXT_PROFILE_CORE
    );

    SDL_GL_SetAttribute(
        SDL_GL_DOUBLEBUFFER,
        1
    );

    SDL_Window* window = SDL_CreateWindow(
        "Sylph Rendering Engine",
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE
    );

    if (window == nullptr)
    {
        std::cerr << "Failed to create SDL window: "
                  << SDL_GetError()
                  << '\n';
        SDL_Quit();
        return EXIT_FAILURE;
    }

    SDL_GLContext gl_context =
        SDL_GL_CreateContext(window);

    if (gl_context == nullptr)
    {
        std::cerr << "Failed to create OpenGL context: "
                  << SDL_GetError()
                  << '\n';
        SDL_DestroyWindow(window);
        SDL_Quit();

        return EXIT_FAILURE;
    }

    // --------------------------------------------------------
    // GLAD
    // --------------------------------------------------------

    if (!gladLoadGL((GLADloadfunc)SDL_GL_GetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD\n";

        SDL_GL_DestroyContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();

        return EXIT_FAILURE;
    }

    if (SDL_GL_SetSwapInterval(1) == 0) {
        // success
    } else {
        std::cerr << "Warning: Unable to set VSync: " << SDL_GetError() << "\n";
    }

    std::cout
        << "OpenGL: "
        << reinterpret_cast<const char*>(
               glGetString(GL_VERSION)
           )
        << '\n';

    // --------------------------------------------------------
    // Renderer
    // --------------------------------------------------------

    render::Renderer renderer;

    renderer.initialize();

    renderer.resize(
        WINDOW_WIDTH,
        WINDOW_HEIGHT
    );

    // --------------------------------------------------------
    // Example
    // --------------------------------------------------------

    examples::HouseExample example(
        WINDOW_WIDTH,
        WINDOW_HEIGHT
    );

    // --------------------------------------------------------
    // Timing
    // --------------------------------------------------------

    Uint64 previous =
        SDL_GetPerformanceCounter();

    const double frequency =
        static_cast<double>(
            SDL_GetPerformanceFrequency()
        );

    // --------------------------------------------------------
    // Main loop
    // --------------------------------------------------------

    bool running = true;

    while (running)
    {
        const Uint64 current =
            SDL_GetPerformanceCounter();

        const float dt =
            static_cast<float>(
                static_cast<double>(
                    current - previous
                ) / frequency
            );

        previous = current;

        // ----------------------------------------------------
        // Events
        // ----------------------------------------------------

        SDL_Event event;

        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_EVENT_QUIT)
            {
                running = false;
            }

            example.process_event(
                event,
                window
            );

            if (event.type ==
                SDL_EVENT_WINDOW_RESIZED)
            {
                const int width =
                    event.window.data1;

                const int height =
                    event.window.data2;

                if (height > 0)
                {
                    renderer.resize(
                        width,
                        height
                    );

                    example.resize(
                        width,
                        height
                    );
                }
            }
        }

        // ----------------------------------------------------
        // Update
        // ----------------------------------------------------

        example.update(dt);

        // ----------------------------------------------------
        // Render
        // ----------------------------------------------------

        renderer.begin_frame(
            Eigen::Vector4f(
                0.02f,
                0.02f,
                0.03f,
                1.0f
            )
        );

        example.render(renderer);

        renderer.end_frame();

        SDL_GL_SwapWindow(window);
    }

    // --------------------------------------------------------
    // Cleanup
    // --------------------------------------------------------

    SDL_GL_DestroyContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return EXIT_SUCCESS;
}