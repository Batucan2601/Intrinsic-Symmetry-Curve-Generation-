#include "imgui/imgui_impl_opengl3.h"
#include "imgui/imgui_stdlib.h"
#include "../Include/ImguiMenuBar.h"
#include "../Include/Sampling.h"
#include "../Include/Laplace-Beltrami.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/CoreTypeDefs.h"
#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#include "../Include/Skeleton.h"
#include "../Include/NLateralDescriptor.h"
#include "../Include/ShapeDiameter.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/HeatKernelSignature.h"

static void imgui_menubar_save_mesh(TrilateralMesh* m );
static void dominant_symmetry_plane(TrilateralMesh* m);
static void skeleton_generation(TrilateralMesh* m);

static bool is_draw_plane = false;
static bool is_draw_skeleton = false;
static Plane plane;
static Skeleton skeleton;

static int dvorak_no_of_significant_points = 0;
void imgui_menu_bar(TrilateralMesh* m)
{
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("Open")) { /* Open action */ }
            if (ImGui::MenuItem("Save")) { /* Save action */ }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Edit"))
        {
            if (ImGui::MenuItem("Undo")) { /* Undo action */ }
            if (ImGui::MenuItem("Redo")) { /* Redo action */ }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Mesh"))
        {
            if (ImGui::MenuItem("Save TrilateralMesh")) { imgui_menubar_save_mesh(m); }
            if (ImGui::BeginMenu("Skeleton "))
            {
                skeleton_generation(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Trilateral "))
            {
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Dvorak"))
            {
                if (ImGui::BeginMenu("Dvorak Extraction of significant points "))
                {
                    if (ImGui::InputInt("no of points ", &dvorak_no_of_significant_points));
                    if (ImGui::MenuItem("Show"))
                    {
                        dvorak_show_signifcant_points(m, dvorak_no_of_significant_points);
                    }

                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            };
            if (ImGui::BeginMenu("Dominant Symmetry Plane"))
            {
                dominant_symmetry_plane(m);
                ImGui::EndMenu();

            }
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
}


//manaul save for now
// do this for spectral
static void imgui_menubar_save_mesh(TrilateralMesh* m )
{
    std::string file_name = "C:\\Users\\BATU\\Desktop\\Master\\Trilateral\\Trilateral\\TrilateralMesh\\off\\spectral_mesh.off";
    std::ofstream outFile(file_name);

    if (!outFile) {
        std::cerr << "Failed to open file for writing!" << std::endl;
        exit(1);
    }
    // Write OFF
    outFile << "OFF" << std::endl;
    outFile << m->vertices.size() << " " << m->triangles.size() / 3 << " " << 0 << std::endl;
    for (size_t i = 0; i < m->vertices.size(); i++)
    {
        outFile << m->vertices[i].x << " " << m->vertices[i].y << " " << m->vertices[i].z << std::endl;
    }
    for (size_t i = 0; i < m->triangles.size(); i += 3)
    {
        outFile << "3 " << m->triangles[i] << " " << m->triangles[i + 1] << " " << m->triangles[i + 2] << std::endl;
    }
    // Close the file
    outFile.close();
}


static void skeleton_generation(TrilateralMesh* m)
{
    if (ImGui::MenuItem("Generate Skeleton"))
    {
        skeleton = skeleton_read_swc_file(m, "0001.isometry.8.swc");
        is_draw_skeleton = true;
    }

}

static float convergence_ratio;
static void dominant_symmetry_plane(TrilateralMesh* m )
{
    ImGui::InputFloat("convergence ratio", &convergence_ratio);
    if (ImGui::MenuItem("Generate Dominant Symmetry Plane "))
    {
        plane = generate_dominant_symmetry_plane(m, convergence_ratio);
        is_draw_plane = true;
    }
    if (ImGui::MenuItem("Save dominant Symmetry Plane"))
    {
        dom_sym_save_plane(plane, m);
    }
    if (ImGui::MenuItem("read dominant Symmetry Plane"))
    {
        dom_sym_read_plane(m, plane);
    }
}





#pragma region drawsection
static void draw_dom_sym();
static void draw_skeleton();

void draw_all()
{
    draw_dom_sym();
    draw_skeleton();
}


static void draw_dom_sym()
{
    if (is_draw_plane)
    {

        DrawTriangle3D(CoreType_conv_glm_raylib_vec3(plane.p1), CoreType_conv_glm_raylib_vec3(plane.p2), CoreType_conv_glm_raylib_vec3(plane.p3), RED);
        DrawTriangle3D(CoreType_conv_glm_raylib_vec3(plane.p1), CoreType_conv_glm_raylib_vec3(plane.p3), CoreType_conv_glm_raylib_vec3(plane.p4), RED);
    }
}

static void draw_skeleton()
{
    if (is_draw_skeleton)
    {
        for (size_t i = 0; i < skeleton.adjacencies.size(); i++)
        {
            std::vector<unsigned int> skeleton_adj_i = skeleton.adjacencies[i];
            for (size_t j = 0; j < skeleton_adj_i.size(); j++)
            {
                DrawLine3D(CoreType_conv_glm_raylib_vec3(skeleton.skeletonFormat[i].point),
                    CoreType_conv_glm_raylib_vec3(skeleton.skeletonFormat[skeleton_adj_i[j]].point), RED);
            }
        }
    }
    
}