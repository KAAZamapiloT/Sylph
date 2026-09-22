// examples/model_renderer_demo.cpp
//
// Sylph renderer demo:
//
//     OBJ -> ObjLoader -> MeshData -> render::Mesh -> Renderer
//
// Usage:
//     syplh_demo.exe
//     syplh_demo.exe assets/meshes/bunny.obj
//
// This example owns the SDL window/OpenGL context. The existing Renderer,
// Mesh and Shader classes only consume the already-created OpenGL context.

#include "render/camera.hpp"
#include "render/camera_controller.hpp"
#include "render/mesh.hpp"
#include "render/obj_loader.hpp"
#include "render/renderer.hpp"
#include "render/shader.hpp"
#include "render/transform.hpp"

#include <Eigen/Dense>

#include <SDL3/SDL.h>
#include <glad/gl.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr int kWidth = 1280;
constexpr int kHeight = 720;

constexpr char kVertexShader[] = R"GLSL(
#version 450 core

layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aUV;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;

out vec3 vWorldNormal;
out vec2 vUV;

void main()
{
    vec4 world = uModel * vec4(aPosition, 1.0);

    vWorldNormal = normalize(
        mat3(transpose(inverse(uModel))) * aNormal
    );

    vUV = aUV;
    gl_Position = uProjection * uView * world;
}
)GLSL";

constexpr char kFragmentShader[] = R"GLSL(
#version 450 core

in vec3 vWorldNormal;
in vec2 vUV;

out vec4 FragColor;

void main()
{
    const vec3 L = normalize(vec3(0.65, 1.0, 0.40));
    const vec3 base_color = vec3(0.68, 0.74, 0.82);
    const float ambient = 0.08;

    const vec3 N = normalize(vWorldNormal);
    const float NdotL = max(dot(N, L), 0.0);

    FragColor = vec4(base_color * (ambient + NdotL), 1.0);
}
)GLSL";

void check_sdl(bool ok, const char* operation)
{
    if (ok)
        return;

    const char* error = SDL_GetError();
    throw std::runtime_error(
        std::string(operation) + ": " +
        (error != nullptr && error[0] != '\0' ? error : "unknown SDL error"));
}

std::filesystem::path resolve_model_path(int argc, char** argv)
{
    namespace fs = std::filesystem;

    const fs::path requested =
        argc > 1 ? fs::path(argv[1]) : fs::path("assets/meshes/girlOBJ.obj");

    // Explicit paths are respected as-is when they already exist.
    if (requested.is_absolute() && fs::exists(requested))
        return fs::weakly_canonical(requested);

    // Search from the process working directory first, then a few common
    // build-directory layouts used by this project.
    const fs::path cwd = fs::current_path();
    const std::vector<fs::path> candidates = {
        requested,
        cwd / requested,
        cwd / ".." / requested,
        cwd / "../.." / requested,
        cwd / "../../.." / requested,
    };

    for (const fs::path& candidate : candidates)
    {
        std::error_code ec;
        if (fs::exists(candidate, ec) && fs::is_regular_file(candidate, ec))
            return fs::weakly_canonical(candidate, ec);
    }

    std::string message = "Model file not found: " + requested.string();
    message += "\nCurrent working directory: " + cwd.string();
    message += "\nPass the OBJ path explicitly, for example:";
    message += "\n  syplh_demo.exe D:/PersonalProj/Sylph/assets/meshes/girlOBJ.obj";
    throw std::runtime_error(message);
}

} // namespace

int main(int argc, char** argv)
{
    SDL_Window* window = nullptr;
    SDL_GLContext context = nullptr;

    try
    {
        const std::filesystem::path model_path =
            resolve_model_path(argc, argv);

        std::cout
            << "Loading model: "
            << model_path.string()
            << '\n';

        // -------------------------------------------------------------
        // 1. SDL + OpenGL context
        // -------------------------------------------------------------

        check_sdl(
            SDL_Init(SDL_INIT_VIDEO),
            "SDL_Init"
        );

        check_sdl(
            SDL_GL_SetAttribute(
                SDL_GL_CONTEXT_MAJOR_VERSION,
                4
            ),
            "SDL_GL_SetAttribute(GL major)"
        );

        check_sdl(
            SDL_GL_SetAttribute(
                SDL_GL_CONTEXT_MINOR_VERSION,
                5
            ),
            "SDL_GL_SetAttribute(GL minor)"
        );

        check_sdl(
            SDL_GL_SetAttribute(
                SDL_GL_CONTEXT_PROFILE_MASK,
                SDL_GL_CONTEXT_PROFILE_CORE
            ),
            "SDL_GL_SetAttribute(GL profile)"
        );

        check_sdl(
            SDL_GL_SetAttribute(
                SDL_GL_DOUBLEBUFFER,
                1
            ),
            "SDL_GL_SetAttribute(double buffer)"
        );

        window = SDL_CreateWindow(
            "Sylph - External OBJ + Renderer",
            kWidth,
            kHeight,
            SDL_WINDOW_OPENGL |
            SDL_WINDOW_RESIZABLE
        );

        if (window == nullptr)
        {
            throw std::runtime_error(
                std::string("SDL_CreateWindow: ") +
                SDL_GetError()
            );
        }

        context = SDL_GL_CreateContext(window);

        if (context == nullptr)
        {
            throw std::runtime_error(
                std::string("SDL_GL_CreateContext: ") +
                SDL_GetError()
            );
        }

        // VSync is optional.
        if (SDL_GL_SetSwapInterval(1) == 0)
        {
            std::cout << "VSync: enabled\n";
        }
        else
        {
            std::cout
                << "VSync: unavailable (continuing)\n";
        }

        // -------------------------------------------------------------
        // 2. Load OpenGL entry points
        // -------------------------------------------------------------

        if (
            gladLoadGL(
                reinterpret_cast<GLADloadfunc>(
                    SDL_GL_GetProcAddress
                )
            ) == 0
        )
        {
            throw std::runtime_error(
                "gladLoadGL failed"
            );
        }

        std::cout
            << "OpenGL vendor   : "
            << reinterpret_cast<const char*>(
                glGetString(GL_VENDOR)
            )
            << '\n';

        std::cout
            << "OpenGL renderer : "
            << reinterpret_cast<const char*>(
                glGetString(GL_RENDERER)
            )
            << '\n';

        std::cout
            << "OpenGL version  : "
            << reinterpret_cast<const char*>(
                glGetString(GL_VERSION)
            )
            << '\n';

        // -------------------------------------------------------------
        // 3. Everything owning GL resources lives in this scope.
        // -------------------------------------------------------------

        {
            const render::MeshData mesh_data =
                render::ObjLoader::load(model_path);

            std::cout
                << "Vertices : "
                << mesh_data.vertices.size()
                << '\n';

            std::cout
                << "Indices  : "
                << mesh_data.indices.size()
                << '\n';

            std::cout
                << "Triangles: "
                << mesh_data.indices.size() / 3
                << '\n';

            if (
                mesh_data.vertices.empty() ||
                mesh_data.indices.empty()
            )
            {
                throw std::runtime_error(
                    "OBJ loader returned an empty mesh"
                );
            }

            // Compute bounding box
            Eigen::Vector3f min_bound = mesh_data.vertices[0].position;
            Eigen::Vector3f max_bound = mesh_data.vertices[0].position;
            for (const auto& v : mesh_data.vertices)
            {
                min_bound = min_bound.cwiseMin(v.position);
                max_bound = max_bound.cwiseMax(v.position);
            }
            const Eigen::Vector3f center = (min_bound + max_bound) * 0.5f;
            const Eigen::Vector3f size = max_bound - min_bound;
            const float max_extent = std::max({size.x(), size.y(), size.z()});

            render::Mesh mesh(
                mesh_data.vertices,
                mesh_data.indices
            );

            render::Shader shader(
                kVertexShader,
                kFragmentShader
            );

            render::Transform transform;
            // Center the model at the origin and scale it to a reasonable size
            transform.set_position(-center * (2.0f / max_extent));
            transform.set_scale(Eigen::Vector3f::Constant(2.0f / max_extent));

            render::Camera camera;

            // Place camera at a reasonable distance
            camera.set_position(
                Eigen::Vector3f(
                    0.0f,
                    0.0f,
                    4.0f
                )
            );

            camera.set_aspect_ratio(
                static_cast<float>(kWidth) /
                static_cast<float>(kHeight)
            );

            camera.set_fov_degrees(60.0f);

            render::Renderer renderer;

            renderer.initialize();
            renderer.resize(
                kWidth,
                kHeight
            );

            // =============================================================
            // CAMERA CONTROLLER
            // =============================================================

            render::CameraController camera_controller;

            camera_controller.set_move_speed(4.0f);
            camera_controller.set_sprint_multiplier(3.0f);
            camera_controller.set_slow_multiplier(0.30f);
            camera_controller.set_mouse_sensitivity(0.0025f);

            // Mouse wheel controls FOV in the example.
            float fov_degrees = camera.fov_degrees();

            bool running = true;

            auto previous_time =
                std::chrono::steady_clock::now();

            // -------------------------------------------------------------
            // Main loop
            // -------------------------------------------------------------

            std::cout
                << "\nControls:\n"
                << "  WASD       Move\n"
                << "  Space/Ctrl Up / Down\n"
                << "  Shift      Sprint\n"
                << "  Alt        Slow\n"
                << "  RMB        Mouse look\n"
                << "  Wheel      FOV / zoom\n"
                << "  R          Reset camera\n"
                << "  Esc        Quit\n\n";

            while (running)
            {
                SDL_Event event{};

                // =========================================================
                // EVENT PROCESSING
                // =========================================================

                while (SDL_PollEvent(&event))
                {
                    // CameraController owns:
                    //   - RMB mouse-look activation
                    //   - relative mouse mode
                    //   - mouse-look motion
                    //   - R reset
                    camera_controller.process_event(
                        event,
                        window
                    );

                    camera_controller.process_mouse_motion(
                        camera,
                        event
                    );

                    switch (event.type)
                    {
                    case SDL_EVENT_QUIT:
                        running = false;
                        break;

                    case SDL_EVENT_KEY_DOWN:
                        if (
                            !event.key.repeat &&
                            event.key.scancode ==
                            SDL_SCANCODE_ESCAPE
                        )
                        {
                            running = false;
                        }
                        break;

                    case SDL_EVENT_MOUSE_WHEEL:
                        // Keep FOV as an example-level camera feature.
                        // Positive wheel -> zoom in -> smaller FOV.
                        fov_degrees -=
                            static_cast<float>(
                                event.wheel.y
                            ) * 3.0f;

                        fov_degrees =
                            std::clamp(
                                fov_degrees,
                                20.0f,
                                100.0f
                            );

                        camera.set_fov_degrees(
                            fov_degrees
                        );
                        break;

                    case SDL_EVENT_WINDOW_RESIZED:
                        renderer.resize(
                            event.window.data1,
                            event.window.data2
                        );

                        if (event.window.data2 > 0)
                        {
                            camera.set_aspect_ratio(
                                static_cast<float>(
                                    event.window.data1
                                ) /
                                static_cast<float>(
                                    event.window.data2
                                )
                            );
                        }
                        break;

                    default:
                        break;
                    }
                }

                // Apply controller-owned edge-triggered actions such as R.
                camera_controller.apply_pending_actions(camera);

                // =========================================================
                // FRAME TIME
                // =========================================================

                const auto now =
                    std::chrono::steady_clock::now();

                float dt =
                    std::chrono::duration<float>(
                        now - previous_time
                    ).count();

                previous_time = now;

                // Prevent giant jumps after a breakpoint/window stall.
                dt = std::clamp(
                    dt,
                    0.0f,
                    0.05f
                );

                // =========================================================
                // CAMERA UPDATE
                // =========================================================

                camera_controller.update(
                    camera,
                    dt
                );

                // =========================================================
                // RENDER
                // =========================================================

                renderer.begin_frame(
                    Eigen::Vector4f(
                        0.035f,
                        0.045f,
                        0.060f,
                        1.0f
                    )
                );

                renderer.draw(
                    mesh,
                    transform,
                    shader,
                    camera
                );

                renderer.end_frame();

                SDL_GL_SwapWindow(window);
            }
        }

        // Mesh + Shader were destroyed while the GL context was alive.
        SDL_GL_DestroyContext(context);
        context = nullptr;

        SDL_DestroyWindow(window);
        window = nullptr;

        SDL_Quit();
        return 0;
    }
    catch (const std::exception& e)
    {
        if (context != nullptr)
        {
            SDL_GL_DestroyContext(context);
            context = nullptr;
        }

        if (window != nullptr)
        {
            SDL_DestroyWindow(window);
            window = nullptr;
        }

        SDL_Quit();

        std::cerr
            << "\nSylph renderer demo failed:\n"
            << e.what()
            << '\n'
            << "\nPress Enter to close...\n";

        std::cin.get();
        return 1;
    }
}
