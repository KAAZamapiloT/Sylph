#version 330 core

in vec3 vNormal;
in vec3 vWorldPos;
in vec3 vCameraPos;

// Spherical Grid Data for maximum t=8 (grid size 8x16 = 128)
uniform vec3 uLightGrid[128];
uniform float uVisGrid[128];
uniform float uGridWeights[128];
uniform vec3 uGridDirs[128];

// BRDF Zonal Harmonics
uniform float uBrdfZH[8];

uniform int uTParam; // 1 to 8
uniform int uGridSize; // current active grid size

out vec4 FragColor;

/*
 * ==============================================================================
 * EDUCATIONAL: EVALUATING THE BRDF ON THE GPU
 * ==============================================================================
 * The BRDF (B) depends on the Reflection vector (R), which changes for every pixel.
 * We could pass the full SH coefficients of B rotated to R, but that would mean 
 * evaluating all 64 SH basis functions per grid point in this shader!
 * 
 * INSTEAD, we use a brilliant math shortcut:
 * A Phong BRDF is circularly symmetric around R. Its SH coefficients are just 
 * Zonal Harmonics (which we pass in 'uBrdfZH'). To evaluate a Zonal Harmonic 
 * lobe at any direction (like our Grid points), we only need to calculate the 
 * unassociated Legendre polynomials P_l(cos(theta)), where cos(theta) is simply 
 * dot(R, GridDirection). 
 * 
 * This reduces evaluating 64 complex 3D functions down to an ultra-fast 1D 
 * polynomial recurrence loop!
 */
// Evaluates the BRDF Zonal Harmonics at the angle cosTheta using Legendre recurrence
float evaluate_brdf(float x) {
    float p0 = 1.0;
    float p1 = x;
    
    float sum = uBrdfZH[0] * p0;
    if (uTParam > 1) sum += uBrdfZH[1] * p1;
    
    float p_prev = p1;
    float p_prev2 = p0;
    
    for (int l = 2; l < uTParam; ++l) {
        float p = ((2.0 * float(l) - 1.0) * x * p_prev - (float(l) - 1.0) * p_prev2) / float(l);
        sum += uBrdfZH[l] * p;
        p_prev2 = p_prev;
        p_prev = p;
    }
    
    // The paper's grid product evaluates functions directly. 
    return max(sum, 0.0);
}

void main()
{
    vec3 N = normalize(vNormal);
    vec3 V = normalize(vCameraPos - vWorldPos);
    vec3 R = normalize(reflect(-V, N));
    
    // THE TRUE PAPER ALGORITHM:
    // Compute the Triple Product Integral (L x B x V) using Spherical Grids!
    vec3 color = vec3(0.0);
    
    /*
     * ==============================================================================
     * EDUCATIONAL: THE TRIPLE PRODUCT INTEGRAL SUM
     * ==============================================================================
     * This is where the magic happens!
     * 1. L_grid[i] = Incident light evaluated at this grid point (from CPU)
     * 2. V_grid[i] = Visibility/Shadow at this grid point (from CPU)
     * 3. evaluate_brdf = The view-dependent Material evaluated dynamically
     * 
     * We just multiply (L * B * V) exactly as if they were simple arrays.
     * Finally, we multiply by the Quadrature Weight and add to our sum.
     * 
     * The sum of these points perfectly equals the DC (average) component of the 
     * full Spherical Harmonic convolution, extracting the exact integral in one line!
     */
    for(int i = 0; i < uGridSize; i++) {
        vec3 dir = uGridDirs[i];
        
        // Only evaluate upper hemisphere for the BRDF
        float cosTheta = dot(R, dir);
        if (cosTheta > 0.0) {
            float b_val = evaluate_brdf(cosTheta);
            color += uLightGrid[i] * b_val * uVisGrid[i] * uGridWeights[i];
        }
    }
    // The integral over the sphere provides the final color (DC component extraction)
    
    // Simple ACES tonemapping
    float a = 2.51; float b = 0.03; float c = 2.43; float d = 0.59; float e = 0.14;
    color = clamp((color*(a*color+b))/(color*(c*color+d)+e), 0.0, 1.0);
    color = pow(color, vec3(1.0 / 2.2));
    
    FragColor = vec4(color, 1.0);
}
