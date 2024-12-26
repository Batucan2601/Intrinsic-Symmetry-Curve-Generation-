/*#include <GL/glew.h>
#include <GLFW/glfw3.h>*/
#include <iostream>
#include "raylib.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "imgui/imgui_impl_opengl3.h"
#include "imgui/imgui_stdlib.h"
#include "imgui/imgui.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "Include/Sampling.h"
#include "Include/Mesh_imgui.h"
#include "Include/ImguiMenuBar.h"
#include "Include/Ray.h"
#include "Include/Camera.h"
#include "rlImGui/rlImGui.h"

#define RLIGHTS_IMPLEMENTATION
#include "raylib/examples/shaders/rlights.h"

//include prototypes
#include "Include/Prototypes.h"
//#include "TrilateralMesh.h"
#include "Include/MeshFactory.h"
#include "Include/Shader.h"
#include "Include/TrilateralMap.h"
#include <raymath.h>
#define GLSL_VERSION            330

static ModifiedCamera camera; 
static void imgui_display_camera(Camera3D& camera);
int main(void) 
{
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(1024, 768, " Trialteral");

    // build and compile our shader program
    //TrilateralMesh m1((char*)"C:\\Users\\Batuhan\\Desktop\\master\\Trilateral\\Trilateral\\Trilateral\\Mesh\\off\\0001.isometry.12.off");
    TrilateralMesh m1((char*)"C:\\Users\\Batuhan\\Desktop\\master\\Trilateral\\Trilateral\\Trilateral\\Mesh\\SCB\\DATA\\SCAPE\\Meshes\\mesh000.off");
    // Load basic lighting shader
    std::string vs_path = RAYLIB_PATH"/examples/shaders/resources/shaders/glsl330/lighting.vs";
    std::string fs_path = RAYLIB_PATH"/examples/shaders/resources/shaders/glsl330/lighting.fs";
    Shader shader = LoadShader(vs_path.c_str(), fs_path.c_str());

    camera.set_camera();
    camera.set_shader(shader);
    camera.set_light();
    int ambientLoc = GetShaderLocation(shader, "ambient");
    float ambient[4] = { 0.1f, 0.1f, 0.1f, 1.0f };
    //SetShaderValue(shader, ambientLoc, ambient, SHADER_UNIFORM_VEC4);
    rlImGuiSetup(true);
    while (!WindowShouldClose())
    {
        BeginDrawing();
            ClearBackground(LIGHTGRAY);
            camera.update();
            BeginMode3D(camera.camera);
                BeginShaderMode(shader);
                    SetShaderValue(shader, ambientLoc, ambient, SHADER_UNIFORM_VEC4);
                    draw_all(&m1 , shader );
                EndShaderMode();
            EndMode3D();
            rlImGuiBegin();
            imgui_menu_bar(&m1);
            imgui_display_camera(camera.camera);
            rlImGuiEnd();

        EndDrawing();
        
  
    }
    rlImGuiShutdown();
    CloseWindow_();
    return 0;
}
static Camera3D temp;
static void imgui_display_camera(Camera3D& camera)
{
    ImGui::Begin("Camera");
    ImGui::Text(" camera position X = %.6f", camera.position.x);
    ImGui::Text(" camera position Y = %.6f", camera.position.y);
    ImGui::Text(" camera position Z = %.6f", camera.position.z);
    ImGui::Text(" Teleport");
    ImGui::InputFloat("X : %f", &temp.position.x);
    ImGui::InputFloat("Y : %f", &temp.position.y);
    ImGui::InputFloat("Z : %f", &temp.position.z);
    if (ImGui::Button("teleport button"))
    {
        camera.position = temp.position;
    }
    ImGui::End();
}