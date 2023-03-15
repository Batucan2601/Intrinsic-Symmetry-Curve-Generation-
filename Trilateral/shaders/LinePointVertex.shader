#version 330 core
layout(location = 0) in vec4 coordinate; // 1 - x - 2 - y 3 - z
uniform mat4 view;
uniform mat4 projection;
uniform mat4 model_1;
uniform mat4 model_2;

//if 4'th points is 1 this means it belongs to the model_1
// if 4't point is 2 this means it belongs to the model_1
out vec3 color; 
void main()
{
    mat4 MVP; 
    if (coordinate[3] < 1 )
    {
        MVP = projection * view * model_1;
        color = vec3(1.0f, 0.0f, 0.0f);
    }
    else if (coordinate[3] >  0)
    {
        MVP = projection * view * model_2;
        color = vec3(0.0f, 1.0f, 0.0f);

    }
    gl_Position = MVP * vec4(coordinate[0] , coordinate[1] , coordinate[2] , 1.0f);
};



