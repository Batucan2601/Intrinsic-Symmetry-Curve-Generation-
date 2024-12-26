#include "../Include/Camera.h"
#include <raymath.h>

void ModifiedCamera::set_camera()
{

    //camera 
    this->camera.position = { 0,0,-1 };
    this->camera.projection = CAMERA_PERSPECTIVE;
    this->camera.target = { 0,0,0 };
    this->camera.fovy = 90;
    this->camera.up = { 0 , 1 ,0 };
}
void ModifiedCamera::set_shader(Shader& shader )
{
    this->shader = shader; 
    shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(shader, "viewPos");
}
void ModifiedCamera::set_light()
{
    //light
    this->light = CreateLight(LIGHT_DIRECTIONAL, this->camera.position, Vector3Zero(), WHITE, shader);
}
void ModifiedCamera::update()
{
    float camera_pos[3] = { this->camera.position.x  , this->camera.position.y  ,this->camera.position.z};
    //update camera 
    UpdateCamera(&this->camera, CAMERA_FREE);
    //update light 
    SetShaderValue(shader, shader.locs[SHADER_LOC_VECTOR_VIEW], camera_pos, SHADER_UNIFORM_VEC3);
    light.position = camera.position;
    UpdateLightValues(shader, light);
}
