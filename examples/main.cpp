#include "glossy/glossy_example.hpp"
#include "sh_lobe/sh_lobe_example.hpp"
#include "shadow/shadow_example.hpp"
#include "oof_demo/oof_example.hpp"
#include "verification/verification_example.hpp"
#include "sh_basis/sh_basis_example.hpp"
#include "product_demo/product_demo_example.hpp"

#include "render/renderer.hpp"

#include <Eigen/Dense>
#include <SDL3/SDL.h>
#include <glad/gl.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

namespace {
constexpr int WINDOW_WIDTH  = 1280;
constexpr int WINDOW_HEIGHT = 720;

// Universal demo interface wrapper
struct IDemo {
    virtual ~IDemo() = default;
    virtual void process_event(const SDL_Event& event, SDL_Window* window) = 0;
    virtual void update(float dt) = 0;
    virtual void resize(int width, int height) = 0;
    virtual void render(render::Renderer& renderer) = 0;
};

template <typename T>
struct DemoWrapper final : IDemo {
    T demo;
    DemoWrapper(int w, int h) : demo(w, h) {}
    void process_event(const SDL_Event& event, SDL_Window* window) override { demo.process_event(event, window); }
    void update(float dt) override { demo.update(dt); }
    void resize(int width, int height) override { demo.resize(width, height); }
    void render(render::Renderer& renderer) override { demo.render(renderer); }
};
}

int main(int argc, char* argv[])
{
    try {
        int demo_choice = 0;
        if (argc > 1) {
            demo_choice = std::atoi(argv[1]);
        }

        if (demo_choice < 1 || demo_choice > 7) {
            std::cout << "\n======================================================\n";
            std::cout << "                 Sylph Demo Selector\n";
            std::cout << "======================================================\n";
            std::cout << " [1] SH Basis Visualizer (Spherical Harmonics Lobe Viewer)\n";
            std::cout << " [2] Interactive SH Product Lobe Demo\n";
            std::cout << " [3] Multi-Lobe Product Visualizer (Shape A x Shape B)\n";
            std::cout << " [4] Glossy PRT Real-Time Relighting (Section 7.3.1)\n";
            std::cout << " [5] Diffuse Shadowed PRT Demo\n";
            std::cout << " [6] Object Occlusion Field (OOF) Dynamic Soft Shadows\n";
            std::cout << " [7] Mathematical Verification & Benchmark Suite (19 Invariants)\n";
            std::cout << "======================================================\n";
            std::cout << "Enter choice [1-7] (default: 7): ";
            std::string input;
            if (std::getline(std::cin, input) && !input.empty()) {
                demo_choice = std::atoi(input.c_str());
            } else {
                demo_choice = 7;
            }
            if (demo_choice < 1 || demo_choice > 7) {
                demo_choice = 7;
            }
        }

        const char* demo_name = "";
        switch (demo_choice) {
            case 1: demo_name = "SH Basis Visualizer"; break;
            case 2: demo_name = "Interactive SH Product Lobe Demo"; break;
            case 3: demo_name = "Multi-Lobe Product Visualizer"; break;
            case 4: demo_name = "Glossy PRT Real-Time Relighting"; break;
            case 5: demo_name = "Diffuse Shadowed PRT Demo"; break;
            case 6: demo_name = "Object Occlusion Field (OOF) Dynamic Soft Shadows"; break;
            case 7:
            default: demo_name = "Mathematical Verification & Benchmark Suite"; break;
        }

        std::cout << "\n[Sylph] Launching: [" << demo_choice << "] " << demo_name << "\n\n";

        if (!SDL_Init(SDL_INIT_VIDEO)) {
            std::cerr << "Failed to initialize SDL: " << SDL_GetError() << "\n";
            return EXIT_FAILURE;
        }

        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

        std::string window_title = std::string("Sylph - [") + std::to_string(demo_choice) + "] " + demo_name;
        SDL_Window* window = SDL_CreateWindow(
            window_title.c_str(),
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
            renderer.initialize();

            std::unique_ptr<IDemo> demo;
            switch (demo_choice) {
                case 1:
                    demo = std::make_unique<DemoWrapper<examples::SHBasisExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
                case 2:
                    demo = std::make_unique<DemoWrapper<examples::SHLobeExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
                case 3:
                    demo = std::make_unique<DemoWrapper<examples::ProductDemoExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
                case 4:
                    demo = std::make_unique<DemoWrapper<examples::GlossyExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
                case 5:
                    demo = std::make_unique<DemoWrapper<examples::ShadowExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
                case 6:
                    demo = std::make_unique<DemoWrapper<examples::OofExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
                case 7:
                default:
                    demo = std::make_unique<DemoWrapper<examples::VerificationExample>>(WINDOW_WIDTH, WINDOW_HEIGHT);
                    break;
            }

            bool running = true;
            Uint64 last_time = SDL_GetPerformanceCounter();

            while (running) {
                SDL_Event event;
                while (SDL_PollEvent(&event)) {
                    if (event.type == SDL_EVENT_QUIT) {
                        running = false;
                    }
                    if (event.type == SDL_EVENT_WINDOW_RESIZED) {
                        demo->resize(event.window.data1, event.window.data2);
                    }
                    demo->process_event(event, window);
                }

                Uint64 current_time = SDL_GetPerformanceCounter();
                float dt = static_cast<float>(current_time - last_time) / static_cast<float>(SDL_GetPerformanceFrequency());
                last_time = current_time;

                demo->update(dt);
                demo->render(renderer);

                SDL_GL_SwapWindow(window);
            }
        }

        SDL_GL_DestroyContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();

        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
