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
#include "Include/Geodesic.h"
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
#include "Include/NLateralMapping.h"
#include "Include/Voronoi.h"
#define GLSL_VERSION            330

static ModifiedCamera camera; 
static void imgui_display_camera(Camera3D& camera, TrilateralMesh* m);

#define CONSOLE_MODE
#ifdef CONSOLE_MODE
int main(int argc, char* argv[]) {
    std::string inputFile;
    std::string outputFileCors;
    std::string outputFileAxis;
    ofstream correspondenceFile;
    ofstream axisPointFile;
    int sampleNo = 3;
    float biggest_dijkstra = 0;
    float fuziness = 0.1;
    float distance_to_mid_param = 0.85;
    float hks_dif_param = 0.08;
    float closeness_param = 0.2;
    int hist_no = 5;
    float voronoi_dif_param = 0.1;
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(1024, 768, " Trilateral");
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-I" && i + 1 < argc) {
            inputFile = argv[++i]; // Get the next argument
        }
        else if (arg == "-OCors" && i + 1 < argc) {
            outputFileCors = argv[++i];
        }
        else if (arg == "-OAxis" && i + 1 < argc) {
            outputFileAxis = argv[++i];
        }
        else if (arg == "-SampleNo" && i + 1 < argc) {
            sampleNo = std::stoi(argv[++i]);
        }
        else {
            std::cerr << "Unknown or incomplete argument: " << arg << std::endl;
        }
    }
    //generate mesh 
    TrilateralMesh m((char*)inputFile.c_str());

    inputFile = inputFile.substr(0, inputFile.size() - m.file_name.size());
    HKS_read_kernel_signature(&m, inputFile);
    // do a single Average Geodesic Pass 
    std::vector<unsigned int> sampled_agd_points = Geodesic_avg_dijkstra_modified(&m, 0.08, 2, false, biggest_dijkstra);
    std::vector<unsigned int> sampled_mgd_points = sampled_agd_points;
    // do mutiple Minimum Geodesic Pass
    for (size_t i = 0; i < sampleNo; i++)
    {
        sampled_mgd_points = Geodesic_min_dijkstra(&m, sampled_agd_points, 0.08, 0.7, false);
    }

    // the main function which does the matching. 
    NLateralMapping_generate_via_voronoi_midpoints(&m, sampled_mgd_points, 0.08, 0.7
        , fuziness, distance_to_mid_param, hks_dif_param, closeness_param, 0.2, hist_no, 0, biggest_dijkstra,
        sampled_agd_points, voronoi_dif_param);

    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output file: " << outputFileCors << std::endl;
    std::cout << "OutputCors file: " << outputFileCors << std::endl;

    // generate the voronoi region. 
    Voronoi voronoi = Voronoi_get_closest_voronoi(&m, voronoi_dif_param);
#ifdef PRUNING_MODE
    //recommended for SCAPE meshes where spill out happens.
    Voronoi_prune_voronoi(&m, voronoi, voronoi_dif_param);
#else
    //a single iteration
    // kill correspondences on same path
    voronoi = Voronoi_destroy_wrong_matches_and_recalculate(&m, voronoi_dif_param, voronoi);
    // generate new pairs. 
    voronoi = Voronoi_check_pair_closeness_and_recalculate(&m, voronoi_dif_param, distance_to_mid_param, voronoi);
#endif 
    correspondenceFile.open(outputFileCors);
    for (size_t i = 0; i < m.calculated_symmetry_pairs.size(); i++)
    {
        correspondenceFile << m.calculated_symmetry_pairs[i].first << " " << m.calculated_symmetry_pairs[i].second << std::endl;
    }

    axisPointFile.open(outputFileAxis);
   
    for (size_t i = 0; i < voronoi.indices.size(); i++)
    {
        axisPointFile << voronoi.indices[i] << std::endl;
    }
    return 0;
}
#else
int main(void) 
{
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(1024, 768, " Trilateral");
    // build and compile our shader program
    TrilateralMesh m1((char*)(SOURCE_PATH + std::string("\\Mesh\\SCB\\DATA\\SCAPE\\Meshes\\mesh000.off")).c_str());
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
            ClearBackground(WHITE);
            camera.update();
            BeginMode3D(camera.camera);
                draw_all(&m1);
                BeginShaderMode(shader);
                    SetShaderValue(shader, ambientLoc, ambient, SHADER_UNIFORM_VEC4);
                    draw_all_shader(&m1, camera.shader);
                EndShaderMode();
            EndMode3D();
            rlImGuiBegin();
            imgui_menu_bar(&m1);
            imgui_display_camera(camera.camera , &m1 );
            rlImGuiEnd();

        EndDrawing();
        
  
    }
    rlImGuiShutdown();
    CloseWindow_();
    return 0;
}
static Camera3D temp;
static int point_no;
static void imgui_display_camera(Camera3D& camera , TrilateralMesh* m )
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
    if (ImGui::InputInt("teleport to point index" , &point_no))
    {
        Vector3 p = { m->vertices[point_no].x , m->vertices[point_no].y ,m->vertices[point_no].z };
        camera.position = p;
    }
    if (ImGui::Button("Give Closest to Camera "))
    {
        glm::vec3 camPos = CoreType_conv_raylib_glm_vec3(camera.position);
        float closest = INFINITY; 
        int closest_index = -1;
        for (size_t i = 0; i < m->vertices.size(); i++)
        {
            float distance = glm::distance(camPos, m->vertices[i]);
            if (closest > distance)
            {
                closest = distance;
                closest_index = i;
            }
        }
        std::cout << " closest index " << closest_index << std::endl; 
    }
    ImGui::End();
}
#endif 
