#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec4 aSH0;
layout(location = 3) in vec4 aSH1;
layout(location = 4) in vec4 aSH2;
layout(location = 5) in vec4 aSH3;

out vec3 vColor;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;

uniform vec3 uBaseColor;

// Environment lighting SH coefficients (order 4 = 16 coefficients)
uniform vec3 uLightSH[16];

void main()
{
    gl_Position = uProjection * uView * uModel * vec4(aPos, 1.0);
    
    // The precomputed transfer function (aSH) integrates visibility and the cosine lobe.
    // So the final diffuse color is just the dot product of Light SH and Transfer SH!
    vec3 diffuse = vec3(0.0);
    
    float sh_arr[16] = float[](
        aSH0.x, aSH0.y, aSH0.z, aSH0.w,
        aSH1.x, aSH1.y, aSH1.z, aSH1.w,
        aSH2.x, aSH2.y, aSH2.z, aSH2.w,
        aSH3.x, aSH3.y, aSH3.z, aSH3.w
    );
    
    for(int i = 0; i < 16; ++i) {
        diffuse += uLightSH[i] * sh_arr[i];
    }
    
    // Prevent negative light
    diffuse = max(diffuse, vec3(0.0));
    
    vColor = uBaseColor * diffuse;
}
