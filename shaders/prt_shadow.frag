#version 330 core
in vec3 vColor;

out vec4 FragColor;

void main()
{
    vec3 color = vColor;
    
    // Simple ACES tonemapping
    float a = 2.51; float b = 0.03; float c = 2.43; float d = 0.59; float e = 0.14;
    color = clamp((color*(a*color+b))/(color*(c*color+d)+e), 0.0, 1.0);
    
    // Gamma correction
    color = pow(color, vec3(1.0 / 2.2));
    
    FragColor = vec4(color, 1.0);
}
