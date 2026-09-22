#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;

out vec3 vColor;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;

void main()
{
    gl_Position = uProjection * uView * uModel * vec4(aPos, 1.0);
    // For the ground plane, the CPU already mathematically solved the 
    // Spherical Harmonic multi-product integral (Light * V_self * V_oof)!
    vColor = aColor;
}
