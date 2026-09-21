#include <iostream>

#include <SDL3/SDL.h>
#include <glad/gl.h>

#include <render/renderer.hpp>

int main()
{
    // --------------------------------------------------------
    // SDL initialization
    // --------------------------------------------------------

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::cerr
            << "SDL_Init failed: "
            << SDL_GetError()
            << '\n';

        return 1;
    }

    // --------------------------------------------------------
    // OpenGL 4.6 Core
    // --------------------------------------------------------

    SDL_GL_SetAttribute(
        SDL_GL_CONTEXT_MAJOR_VERSION,
        4
    );

    SDL_GL_SetAttribute(
        SDL_GL_CONTEXT_MINOR_VERSION,
        6
    );

    SDL_GL_SetAttribute(
        SDL_GL_CONTEXT_PROFILE_MASK,
        SDL_GL_CONTEXT_PROFILE_CORE
    );

    // --------------------------------------------------------
    // Window
    // --------------------------------------------------------

    SDL_Window* window = SDL_CreateWindow(
        "Sylph Renderer",
        1280,
        720,
        SDL_WINDOW_OPENGL
    );

    if (!window) {
        std::cerr
            << "SDL_CreateWindow failed: "
            << SDL_GetError()
            << '\n';

        SDL_Quit();
        return 1;
    }

    // --------------------------------------------------------
    // OpenGL context
    // --------------------------------------------------------

    SDL_GLContext context =
        SDL_GL_CreateContext(window);

    if (!context) {
        std::cerr
            << "SDL_GL_CreateContext failed: "
            << SDL_GetError()
            << '\n';

        SDL_DestroyWindow(window);
        SDL_Quit();

        return 1;
    }

    // --------------------------------------------------------
    // GLAD
    // --------------------------------------------------------

    if (!gladLoadGL(
            (GLADloadfunc)SDL_GL_GetProcAddress))
    {
        std::cerr
            << "Failed to initialize GLAD\n";

        SDL_GL_DestroyContext(context);
        SDL_DestroyWindow(window);
        SDL_Quit();

        return 1;
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

    renderer.resize(1280, 720);

    // --------------------------------------------------------
    // Main loop
    // --------------------------------------------------------

    bool running = true;
   

    float counter=0;
    float SI=0;
    while (running) {

        SDL_Event event;

        while (SDL_PollEvent(&event)) {

            if (event.type == SDL_EVENT_QUIT) {
                running = false;
            }
        }

        // Clear the screen.
        //
        // glClearColor uses [0, 1], not [0, 255].
        renderer.begin_frame({
            0.1f,
            SI,
            std::tan(SI),
            1.0f
        });
        counter+=0.001f;
        renderer.end_frame();
        SI=std::sin(counter);
        // Renderer does not swap the window.
        SDL_GL_SwapWindow(window);
    }

    // --------------------------------------------------------
    // Shutdown
    // --------------------------------------------------------

    SDL_GL_DestroyContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}