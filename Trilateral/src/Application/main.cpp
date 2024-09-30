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

int main(void) 
{
    InitWindow(1024, 768, " Trialteral");

    // build and compile our shader program


    TrilateralMesh m1((char*)"C:\\Users\\Batuhan\\Desktop\\master\\Trilateral\\Trilateral\\Trilateral\\Mesh\\off\\0001.isometry.8.off");


    Camera camera;
    camera.position = { 0,0,-40 };
    camera.projection= CAMERA_PERSPECTIVE;
    camera.target = {0,0,0};
    camera.fovy = 90;
    camera.up = { 0 , 1 ,0 };

    Matrix identity = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    
    rlImGuiSetup(true);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(GREEN);
        UpdateCamera(&camera, CAMERA_FREE);
        BeginMode3D(camera);
        DrawMesh( m1.raylib_mesh, LoadMaterialDefault(), identity);
        draw_all();
        
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
        // end ImGui Content
        rlImGuiEnd();

        EndDrawing();
        
  
    }
    rlImGuiShutdown();
    CloseWindow_();
    return 0;
}

