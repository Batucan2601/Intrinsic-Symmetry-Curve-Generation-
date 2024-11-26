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
#include "rlImGui/rlImGui.h"


//include prototypes
#include "Include/Prototypes.h"
//#include "TrilateralMesh.h"
#include "Include/MeshFactory.h"
#include "Include/Shader.h"
#include "Include/TrilateralMap.h"

static void imgui_display_camera(Camera3D& camera);
int main(void) 
{
    InitWindow(1024, 768, " Trialteral");

    // build and compile our shader program
    //TrilateralMesh m1((char*)"C:\\Users\\Batuhan\\Desktop\\master\\Trilateral\\Trilateral\\Trilateral\\Mesh\\off\\0002.isometry.1.off");
    TrilateralMesh m1((char*)"C:\\Users\\Batuhan\\Desktop\\master\\Trilateral\\Trilateral\\Trilateral\\Mesh\\SCB\\DATA\\SCAPE\\Meshes\\mesh000.off");


    Camera camera;
    camera.position = { 0,0,-1 };
    camera.projection= CAMERA_PERSPECTIVE;
    camera.target = {0,0,0};
    camera.fovy = 90;
    camera.up = { 0 , 1 ,0 };

    
    rlImGuiSetup(true);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(LIGHTGRAY);
        UpdateCamera(&camera, CAMERA_FREE);
        BeginMode3D(camera);

        draw_all(&m1);
        EndMode3D();

        rlImGuiBegin();
        // show ImGui Content
        bool open = true;
        //ImGui::ShowDemoWindow(&open);
        /*imgui_mesh_window(selected_mesh, mesh_fac);
        imgui_selected_mesh_properties_window(selected_mesh, mesh_fac);
        imgui_KIDS_skeleton(selected_mesh, mesh_fac);
        imgui_N_Lateral_Parameters(selected_mesh, mesh_fac);
        imgui_debug_layer(selected_mesh, mesh_fac, cameraPos, cameraFront, cameraUp); */
        imgui_menu_bar(&m1);
        imgui_display_camera(camera);
        // end ImGui Content
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