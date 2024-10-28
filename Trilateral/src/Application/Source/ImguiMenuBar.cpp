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
#include "../Include/SpinImage.h"
#include "../Include/KIDS.h"
#include "ImGuiFileDialog.h"
#include "raymath.h"
#include <fstream>  // Include for file stream handling


static void imgui_menubar_save_mesh(TrilateralMesh* m );
static void dominant_symmetry_plane(TrilateralMesh* m);
static void skeleton_generation(TrilateralMesh* m);
static void trilateral_functions(TrilateralMesh* m);
static void dvorak_functions(TrilateralMesh* m);
static void distribution_functions(TrilateralMesh* m);
static void KIDS_dataset(TrilateralMesh* m);
static void mesh_drawing();
static void drawFileDialog(std::string& file_path, std::string& file_path_name, std::string file_type, bool& is_open);
static void display_file_dialogs(TrilateralMesh* m);
static void  display_descriptor(TrilateralMesh* m);
static void  display_descriptor_all(TrilateralMesh* m);
std::string file_path;
std::string file_path_name = "";
static bool is_searching_swc = false; 
static bool is_writing_dsc = false; 
static bool is_reading_dsc = false; 
static bool is_writing_pln = false;
static bool is_reading_pln = false;
static bool is_mesh_wires = false;
static bool is_draw_plane = false;
static bool is_draw_skeleton = false;
static bool is_draw_mesh = true;
static int descriptor_no = 0;
static Plane plane;
static Skeleton skeleton;
static int dvorak_no_of_significant_points = 0;
static int no_of_dist_points = 0;
static float dvorak_geodesic_dist_param = 0;
static int mesh_index = 0;
std::vector<TrilateralDescriptor> positive_desc; 
std::vector<TrilateralDescriptor> negative_desc; 
std::vector<unsigned int> distributed_indices;

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
            if (ImGui::BeginMenu("Mesh properties "))
            {
                if (!is_mesh_wires)
                {
                    if (ImGui::Button("Draw mesh with wireframe "))
                    {
                        is_mesh_wires = !is_mesh_wires;
                    }
                }
                else
                {
                    if (ImGui::Button("Draw mesh with polygon "))
                    {
                        is_mesh_wires = !is_mesh_wires;

                    }
                }
          
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Skeleton "))
            {
                skeleton_generation(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Trilateral "))
            {
                trilateral_functions(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Dvorak"))
            {
                dvorak_functions(m);
               
                ImGui::EndMenu();
            };
            if (ImGui::BeginMenu("Dominant Symmetry Plane"))
            {
                dominant_symmetry_plane(m);
                ImGui::EndMenu();

            }
            if (ImGui::BeginMenu("Distributions"))
            {
                distribution_functions(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Mesh Drawing"))
            {
                mesh_drawing();
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("dataset"))
            {
                KIDS_dataset(m);
                ImGui::EndMenu();

            }
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
    display_file_dialogs(m);
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
        is_searching_swc = true; 
    }
    if (ImGui::MenuItem("Show Skeleton End points on Mesh "))
    {
        std::vector<unsigned int> vertices;
        skeleton_get_end_points_update_mesh(m, skeleton, vertices);
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
        is_writing_pln = true; 
    }
    if (ImGui::MenuItem("read dominant Symmetry Plane"))
    {
        is_reading_pln = true;
    }
}

static void trilateral_functions(TrilateralMesh* m)
{
    if (ImGui::MenuItem("Setup mesh for operations  "))
    {
        std::cout << "reading symmetries " << std::endl;
        read_symmetry_format((char*)"../../Trilateral/Mesh/off/sym.txt", m);
        std::cout << "reading symmetries DONE " << std::endl;
        std::cout << "reading respective heat kernel signature" << std::endl;
        HKS_read_kernel_signature(m);
        std::cout << "reading respective heat kernel signature DONE " << std::endl;
        std::cout << "select skeleton " << std::endl;
        is_searching_swc = true;
    }
    if (ImGui::BeginMenu("Trilateral Descriptor"))
    {
        if (ImGui::MenuItem("Save trilateral desscriptors"))
        {
            is_writing_dsc = true;
        }
        if (ImGui::MenuItem("Load trilateral descriptors"))
        {
            is_reading_dsc = true;
        }
        if (ImGui::InputInt("descriptor no ", &descriptor_no));
        if (ImGui::MenuItem("Display descriptor"))
        {
            display_descriptor(m);
        }
        if (ImGui::MenuItem("Display descriptor all"))
        {
            display_descriptor_all(m);
        }
        ImGui::EndMenu();

    }
    if (ImGui::MenuItem("End point matching with Dvorak significant poins triangle area "))
    {
        trilateral_point_matching_with_dvorak_endpoints(m, positive_desc, negative_desc, plane , dvorak_no_of_significant_points , convergence_ratio);
        is_draw_plane = true;
    }
    if (ImGui::MenuItem("End point matching with Dvorak significant poins Optimal transform"))
    {
        trilateral_point_matching_with_gaussian_endpoints_and_OT(m, positive_desc, negative_desc, plane, dvorak_no_of_significant_points, convergence_ratio);
        is_draw_plane = true;
    }
    if (ImGui::MenuItem("End point matching with Dvorak significant poins Optimal transform w CDF"))
    {
        trilateral_point_matching_with_gaussian_endpoints_and_OT_w_CDF(m, positive_desc, negative_desc, plane, dvorak_no_of_significant_points, convergence_ratio);
        is_draw_plane = true;
    }
    if (ImGui::MenuItem("Generate trilaterals using endpoints "))
    {
        trilateral_display_trilateral_from_skeleton_endpoints(m, positive_desc , negative_desc , skeleton , plane );
    }
}

static void KIDS_dataset(TrilateralMesh* m)
{
    if (ImGui::MenuItem("read KIDS"))
    {
        KIDS_read_meshes();
    }
    ImGui::InputInt("Mesh Index", &mesh_index);
    if (ImGui::MenuItem("display mesh "))
    {
        KIDS_select_mesh(*m , mesh_index);
    }
    if (ImGui::MenuItem("generate KIDS planes"))
    {
        KIDS_dom_sym_generate_or_read_planes(convergence_ratio);
    }
    if (ImGui::MenuItem(" select N curvature points "))
    {
        KIDS_generate_gaussians(dvorak_no_of_significant_points, dvorak_geodesic_dist_param);
    }
    if (ImGui::MenuItem("use gaussian for endpoint matching "))
    {
        KIDS_endpoint_matching_w_gaussian(dvorak_no_of_significant_points, convergence_ratio);
    }
   
}
#pragma region draw section
static void draw_dom_sym();
static void draw_skeleton();
static void draw_resemblance_pairs(TrilateralMesh* m );
static void draw_mesh(TrilateralMesh* m);

void draw_all(TrilateralMesh* m)
{
    draw_dom_sym();
    draw_skeleton();
    draw_mesh(m);
    draw_resemblance_pairs(m);
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
                    CoreType_conv_glm_raylib_vec3(skeleton.skeletonFormat[skeleton_adj_i[j]].point), skeleton.skeletonFormat[i].color);
            }
        }
    }
    
}

static void draw_resemblance_pairs(TrilateralMesh* m )
{
    for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
    {
        int index1 = m->calculated_symmetry_pairs[i].first;
        int index2 = m->calculated_symmetry_pairs[i].second;
        DrawLine3D(CoreType_conv_glm_raylib_vec3(m->vertices[index1]),
            CoreType_conv_glm_raylib_vec3(m->vertices[index2]), BLUE);
    }
}
static void draw_mesh(TrilateralMesh* m)
{
    if (is_draw_mesh)
    {
        if (!is_mesh_wires)
        {
            DrawMesh(m->raylib_mesh, LoadMaterialDefault(), MatrixIdentity());
        }
        else
        {
            DrawMeshWires(m->raylib_mesh, Vector3{ 0,0,0 }, 1, BLACK);
        }
    }
    
}
static void drawFileDialog(std::string& file_path , std::string& file_path_name , std::string file_type , bool& is_open) {
    if (!is_open)
    {
        return; 
    }
    // open Dialog Simple
    IGFD::FileDialogConfig config;
    config.path = ".";
    ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose File", file_type.c_str(), config);
    // display
    if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey")) {
        if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
            file_path_name = ImGuiFileDialog::Instance()->GetFilePathName();
            file_path = ImGuiFileDialog::Instance()->GetCurrentPath();
            // action
        }

        // close
        ImGuiFileDialog::Instance()->Close();
        is_open = false;
    }
}

static void display_file_dialogs(TrilateralMesh* m )
{
    drawFileDialog(file_path, file_path_name, ".swc", is_searching_swc);
    if (file_path_name != "")
    {
        skeleton = skeleton_read_swc_file(m, file_path_name);
        is_draw_skeleton = true;
        file_path_name = "";
    }
    drawFileDialog(file_path, file_path_name, ".dsc", is_writing_dsc);
    if (file_path_name != "")
    {
        TrilateralDescriptor_write(file_path_name , positive_desc , negative_desc);
        file_path_name = "";
    }
    drawFileDialog(file_path, file_path_name, ".dsc", is_reading_dsc);
    if (file_path_name != "")
    {
        TrilateralDescriptor_read(file_path_name , positive_desc , negative_desc);
        file_path_name = "";
    }
    drawFileDialog(file_path, file_path_name, ".pln", is_writing_pln);
    if (file_path_name != "")
    {
        dom_sym_write_plane(m, plane,  file_path_name);
        is_draw_plane = true;
        file_path_name = "";
    }
    drawFileDialog(file_path, file_path_name, ".pln", is_reading_pln);
    if (file_path_name != "")
    {
        dom_sym_read_plane(m, plane, file_path_name);
        is_draw_plane = true;
        file_path_name = "";
    }
}

static void mesh_drawing()
{
    if (is_draw_plane)
    {
        if (ImGui::MenuItem("Disable Plane Drawing"))
        {
            is_draw_plane = !is_draw_plane;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable Plane Drawing "))
        {
            is_draw_plane = !is_draw_plane;
        }
    }
    if (is_draw_skeleton)
    {
        if (ImGui::MenuItem("Disable Skeleton Drawing"))
        {
            is_draw_skeleton = !is_draw_skeleton;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable Skeleton Drawing "))
        {
            is_draw_skeleton = !is_draw_skeleton;
        }
    }

    if (is_draw_mesh)
    {
        if (ImGui::MenuItem("Disable Mesh Drawing"))
        {
            is_draw_mesh = !is_draw_mesh;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable mesh Drawing "))
        {
            is_draw_mesh = !is_draw_mesh;
        }
    }

    
}
static void dvorak_functions(TrilateralMesh* m)
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
    if (ImGui::BeginMenu("Dvorak Sweep points with geodesic dist "))
    {
        if (ImGui::InputFloat("geodesic dist ", &dvorak_geodesic_dist_param));
        if (ImGui::MenuItem("start sweep"))
        {
            std::vector<DvorakPairs> dvorak_pairs = dvorak_extraction_of_significant_points(m, dvorak_no_of_significant_points);
            dvorak_distance_sweep(m, dvorak_pairs , dvorak_geodesic_dist_param);

        }
        ImGui::EndMenu();

    }
}
static void distribution_functions(TrilateralMesh* m)
{ 
    if (ImGui::BeginMenu("Furthest Point Sampling"))
    {
        if (ImGui::InputInt("no of points ", &no_of_dist_points));
        if (ImGui::MenuItem("Furthest point Sampling "))
        {
           distributed_indices = furthest_point_sampling(m, no_of_dist_points, true);
        }

        ImGui::EndMenu();
    }

}
#pragma region trilateral save load 



static void  display_descriptor(TrilateralMesh* m )
{
    TrilateralDescriptor desc; 
    if (descriptor_no < positive_desc.size())
    {
        desc = positive_desc[descriptor_no];
    }
    else
    {
        desc = negative_desc[descriptor_no - positive_desc.size()];
    }

    for (size_t i = 0; i < m->vertices.size(); i++)
    {
        m->raylib_mesh.colors[i * 4] = 0;
        m->raylib_mesh.colors[i * 4+1] = 0;
        m->raylib_mesh.colors[i * 4+2] = 0;
        m->raylib_mesh.colors[i * 4+3] = 255;
    }

    for (size_t i = 0; i < desc.visited_indices.size() ; i++)
    {
        int index = desc.visited_indices[i];

        m->raylib_mesh.colors[index * 4] = 0;
        m->raylib_mesh.colors[index * 4 + 1] = 255;
        m->raylib_mesh.colors[index * 4 + 2] = 0;
        m->raylib_mesh.colors[index * 4 + 3] = 255;
    }
    for (size_t i = 0; i < desc.path_1_2.size(); i++)
    {
        int index = desc.path_1_2[i];
        m->raylib_mesh.colors[index * 4] = 255;
        m->raylib_mesh.colors[index * 4 + 1] = 0;
        m->raylib_mesh.colors[index * 4 + 2] = 0;
        m->raylib_mesh.colors[index * 4 + 3] = 255;
    }
    for (size_t i = 0; i < desc.path_1_3.size(); i++)
    {
        int index = desc.path_1_3[i];
        m->raylib_mesh.colors[index * 4] = 255;
        m->raylib_mesh.colors[index * 4 + 1] = 0;
        m->raylib_mesh.colors[index * 4 + 2] = 0;
        m->raylib_mesh.colors[index * 4 + 3] = 255;
    }
    for (size_t i = 0; i < desc.path_2_3.size(); i++)
    {
        int index = desc.path_2_3[i];
        m->raylib_mesh.colors[index * 4] = 255;
        m->raylib_mesh.colors[index * 4 + 1] = 0;
        m->raylib_mesh.colors[index * 4 + 2] = 0;
        m->raylib_mesh.colors[index * 4 + 3] = 255;
    }
    m->raylib_mesh.colors[desc.p1 * 4] = 0;
    m->raylib_mesh.colors[desc.p1 * 4 + 1] = 0;
    m->raylib_mesh.colors[desc.p1 * 4 + 2] = 255;
    m->raylib_mesh.colors[desc.p1 * 4 + 3] = 255;
    m->update_raylib_mesh();
}

static void  display_descriptor_all(TrilateralMesh* m)
{
    int N_pos = positive_desc.size();
    int N_neg = negative_desc.size();
    //first color inside
    for (size_t i = 0; i < N_pos; i++)
    {
        for (size_t j = 0; j < positive_desc[i].visited_indices.size(); j++)
        {
            int index = positive_desc[i].visited_indices[j];
            m->raylib_mesh.colors[index * 4] = 0;
            m->raylib_mesh.colors[index * 4 + 1] = 255;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
    }
    for (size_t i = 0; i < N_neg; i++)
    {
        for (size_t j = 0; j < negative_desc[i].visited_indices.size(); j++)
        {
            int index = negative_desc[i].visited_indices[j];
            m->raylib_mesh.colors[index * 4] = 0;
            m->raylib_mesh.colors[index * 4 + 1] = 255;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
    }
    for (size_t i = 0; i < N_pos; i++)
    {
        for (size_t j = 0; j < positive_desc[i].path_1_2.size(); j++)
        {
            int index = positive_desc[i].path_1_2[j];
            m->raylib_mesh.colors[index * 4] = 255;
            m->raylib_mesh.colors[index * 4 + 1] = 0;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
        for (size_t j = 0; j < positive_desc[i].path_1_3.size(); j++)
        {
            int index = positive_desc[i].path_1_3[j];
            m->raylib_mesh.colors[index * 4] = 255;
            m->raylib_mesh.colors[index * 4 + 1] = 0;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
        for (size_t j = 0; j < positive_desc[i].path_2_3.size(); j++)
        {
            int index = positive_desc[i].path_2_3[j];
            m->raylib_mesh.colors[index * 4] = 255;
            m->raylib_mesh.colors[index * 4 + 1] = 0;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
    }
    for (size_t i = 0; i < N_neg; i++)
    {
        for (size_t j = 0; j < negative_desc[i].path_1_2.size(); j++)
        {
            int index = negative_desc[i].path_1_2[j];
            m->raylib_mesh.colors[index * 4] = 255;
            m->raylib_mesh.colors[index * 4 + 1] = 0;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
        for (size_t j = 0; j < negative_desc[i].path_1_3.size(); j++)
        {
            int index = negative_desc[i].path_1_3[j];
            m->raylib_mesh.colors[index * 4] = 255;
            m->raylib_mesh.colors[index * 4 + 1] = 0;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
        for (size_t j = 0; j < negative_desc[i].path_2_3.size(); j++)
        {
            int index = negative_desc[i].path_2_3[j];
            m->raylib_mesh.colors[index * 4] = 255;
            m->raylib_mesh.colors[index * 4 + 1] = 0;
            m->raylib_mesh.colors[index * 4 + 2] = 0;
            m->raylib_mesh.colors[index * 4 + 3] = 255;
        }
    }


    m->update_raylib_mesh();

}