/**
 * @file camera_controller.hpp
 * @brief Interactive First-Person Camera Controller.
 * 
 * Captures user input (mouse for look, keyboard for movement) and applies it to 
 * the Camera object. Supports sprinting, pitch clamping, and smooth movement vectors
 * for navigating the 3D environment.
 */
#pragma once

#include "render/camera.hpp"

#include <SDL3/SDL.h>

#include <algorithm>

namespace render
{

class CameraController
{
public:
    void process_event(const SDL_Event& event, SDL_Window* window)
    {
        if (window == nullptr)
            return;

        if (event.type == SDL_EVENT_MOUSE_BUTTON_DOWN &&
            event.button.button == SDL_BUTTON_RIGHT)
        {
            mouse_look_active_ = true;
            SDL_SetWindowRelativeMouseMode(window, true);
            return;
        }

        if (event.type == SDL_EVENT_MOUSE_BUTTON_UP &&
            event.button.button == SDL_BUTTON_RIGHT)
        {
            mouse_look_active_ = false;
            SDL_SetWindowRelativeMouseMode(window, false);
            return;
        }

        if (event.type == SDL_EVENT_KEY_DOWN && !event.key.repeat)
        {
            if (event.key.scancode == SDL_SCANCODE_R)
                reset_requested_ = true;
        }
    }

    void update(Camera& camera, float dt) const
    {
        const bool* keys = SDL_GetKeyboardState(nullptr);

        Eigen::Vector3f movement = Eigen::Vector3f::Zero();

        if (keys[SDL_SCANCODE_W]) movement.z() += 1.0f;
        if (keys[SDL_SCANCODE_S]) movement.z() -= 1.0f;
        if (keys[SDL_SCANCODE_D]) movement.x() += 1.0f;
        if (keys[SDL_SCANCODE_A]) movement.x() -= 1.0f;
        if (keys[SDL_SCANCODE_SPACE]) movement.y() += 1.0f;
        if (keys[SDL_SCANCODE_LCTRL]) movement.y() -= 1.0f;

        if (movement.squaredNorm() > 0.0f)
            movement.normalize();

        float speed = move_speed_;

        if (keys[SDL_SCANCODE_LSHIFT] || keys[SDL_SCANCODE_RSHIFT])
            speed *= sprint_multiplier_;

        if (keys[SDL_SCANCODE_LALT] || keys[SDL_SCANCODE_RALT])
            speed *= slow_multiplier_;

        camera.move_local(movement * speed * dt);
    }

    void process_mouse_motion(Camera& camera, const SDL_Event& event) const
    {
        if (!mouse_look_active_ ||
            event.type != SDL_EVENT_MOUSE_MOTION)
            return;

        camera.rotate(
            -event.motion.yrel * mouse_sensitivity_,
            event.motion.xrel * mouse_sensitivity_);
    }

    void apply_pending_actions(Camera& camera)
    {
        if (!reset_requested_)
            return;

        camera.reset();
        reset_requested_ = false;
    }

    void set_move_speed(float speed)
    {
        move_speed_ = std::max(0.0f, speed);
    }

    void set_sprint_multiplier(float multiplier)
    {
        sprint_multiplier_ = std::max(1.0f, multiplier);
    }

    void set_slow_multiplier(float multiplier)
    {
        slow_multiplier_ = std::clamp(multiplier, 0.01f, 1.0f);
    }

    void set_mouse_sensitivity(float sensitivity)
    {
        mouse_sensitivity_ = std::max(0.0f, sensitivity);
    }

    [[nodiscard]] bool mouse_look_active() const noexcept
    {
        return mouse_look_active_;
    }

private:
    bool mouse_look_active_{false};
    bool reset_requested_{false};

    float move_speed_{4.0f};
    float sprint_multiplier_{3.0f};
    float slow_multiplier_{0.30f};
    float mouse_sensitivity_{0.0025f};
};

} // namespace render
