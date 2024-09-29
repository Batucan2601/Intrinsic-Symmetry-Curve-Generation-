#include "imgui/imgui_impl_glfw.h"
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

static void imgui_menubar_save_mesh(int& selected_mesh, MeshFactory& mesh_fac);
static void dominant_symmetry_plane(MeshFactory& mesh_fac, int selected_mesh);

static int dvorak_no_of_significant_points = 0;
static Plane plane;
void imgui_menu_bar(int& selected_mesh, MeshFactory& mesh_fac)
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
        if (ImGui::BeginMenu("TrilateralMesh"))
        {
            if (ImGui::MenuItem("Save TrilateralMesh")) { imgui_menubar_save_mesh(selected_mesh, mesh_fac); }
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
                        TrilateralMesh* m = &mesh_fac.mesh_vec[selected_mesh];
                        dvorak_show_signifcant_points(m, dvorak_no_of_significant_points);
                        mesh_fac.remove_all();
                        mesh_fac.add_all();
                    }

                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            };
            if (ImGui::BeginMenu("Dominant Symmetry Plane"))
            {
                dominant_symmetry_plane(mesh_fac, selected_mesh);
                ImGui::EndMenu();

            }
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
}


//manaul save for now
// do this for spectral
static void imgui_menubar_save_mesh(int& selected_mesh, MeshFactory& mesh_fac)
{
    TrilateralMesh* m = &mesh_fac.mesh_vec[selected_mesh];
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


static float convergence_ratio;
static void dominant_symmetry_plane(MeshFactory& mesh_fac , int selected_mesh)
{
    ImGui::InputFloat("convergence ratio", &convergence_ratio);
    if (ImGui::MenuItem("Generate Dominant Symmetry Plane "))
    {
        plane = generate_dominant_symmetry_plane(selected_mesh, mesh_fac, convergence_ratio);
    }
    if (ImGui::MenuItem("Save dominant Symmetry Plane"))
    {
        dom_sym_save_plane(plane, &mesh_fac.mesh_vec[selected_mesh]);
    }
    if (ImGui::MenuItem("read dominant Symmetry Plane"))
    {
        dom_sym_read_plane(mesh_fac, selected_mesh, plane);
    }
}
