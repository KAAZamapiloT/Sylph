#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormalOrColor; // For floor this is color, for walls this is normal!

out vec3 vColor;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;

uniform bool uIsFloor;
uniform vec3 uBaseColor;

// SH Lighting coefficients for diffuse lighting
uniform vec3 uLightSH[9];

// Basic spherical harmonic basis evaluation for l<=2
float eval_sh(int l, int m, vec3 n) {
    if (l == 0 && m == 0) return 0.28209479177;
    if (l == 1 && m == -1) return 0.4886025119 * n.y;
    if (l == 1 && m == 0) return 0.4886025119 * n.z;
    if (l == 1 && m == 1) return 0.4886025119 * n.x;
    if (l == 2 && m == -2) return 1.09254843059 * n.x * n.y;
    if (l == 2 && m == -1) return 1.09254843059 * n.y * n.z;
    if (l == 2 && m == 0) return 0.31539156525 * (3.0 * n.z * n.z - 1.0);
    if (l == 2 && m == 1) return 1.09254843059 * n.x * n.z;
    if (l == 2 && m == 2) return 0.54627421529 * (n.x * n.x - n.y * n.y);
    return 0.0;
}

void main()
{
    gl_Position = uProjection * uView * uModel * vec4(aPos, 1.0);
    
    if (uIsFloor) {
        // Floor already has CPU-evaluated SH Multiple Product Integral (Color)
        vColor = aNormalOrColor;
    } else {
        // For rigid objects (Bunny/Walls), evaluate simple SH diffuse lighting
        vec3 worldNormal = normalize(mat3(uModel) * aNormalOrColor);
        vec3 diffuse = vec3(0.0);
        
        // SH projection of cosine lobe (cl = sqrt(4pi/(2l+1)) * A_l)
        float A0 = 1.0;
        float A1 = 2.0 / 3.0;
        float A2 = 0.25;
        float coeffs[9] = float[](A0, A1, A1, A1, A2, A2, A2, A2, A2);
        
        int idx = 0;
        for (int l = 0; l <= 2; ++l) {
            for (int m = -l; m <= l; ++m) {
                float sh_val = eval_sh(l, m, worldNormal) * coeffs[idx];
                diffuse += uLightSH[idx] * sh_val;
                idx++;
            }
        }
        
        vColor = uBaseColor * diffuse;
    }
}
