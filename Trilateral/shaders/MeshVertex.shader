#version 330 core
layout(location = 0) in vec3 coordinate; // 1 - x - 2 - y 3 - z
layout(location = 1) in vec3 color; // 1- r 2 -g 3 -b 
uniform mat4 u_MVP;
out vec3 o_color;
void main()
{
    gl_Position = u_MVP *  vec4(coordinate.x, coordinate.y, coordinate.z, 1.0);
    o_color = color;
};

