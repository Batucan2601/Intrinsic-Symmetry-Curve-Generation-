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
#include "ImGuiFileDialog.h"
#include "raymath.h"
#include <fstream>  // Include for file stream handling


static void imgui_menubar_save_mesh(TrilateralMesh* m );
static void dominant_symmetry_plane(TrilateralMesh* m);
static void skeleton_generation(TrilateralMesh* m);
static void trilateral_functions(TrilateralMesh* m);
static void mesh_drawing();
static void drawFileDialog(std::string& file_path, std::string& file_path_name, std::string file_type, bool& is_open);
static void display_file_dialogs(TrilateralMesh* m);
static void save_trilateral_descriptors();
static void load_trilateral_descriptors();
static void write_trilateral_descriptors(std::string filename);
static void read_trilateral_descriptors(std::string filename);
static void  display_descriptor(TrilateralMesh* m);
std::string file_path;
std::string file_path_name = "";
static bool is_searching_swc = false; 
static bool is_writing_dsc = false; 
static bool is_reading_dsc = false; 
static bool is_mesh_wires = false;
static bool is_draw_plane = false;
static bool is_draw_skeleton = false;
static bool is_draw_mesh = true;
static int descriptor_no = 0;
static Plane plane;
static Skeleton skeleton;
static int dvorak_no_of_significant_points = 0;
std::vector<TrilateralDescriptor> positive_desc; 
std::vector<TrilateralDescriptor> negative_desc; 
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
            if (ImGui::BeginMenu("Mesh Drawing"))
            {
                mesh_drawing();
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
            save_trilateral_descriptors();
        }
        if (ImGui::MenuItem("Load trilateral descriptors"))
        {
            load_trilateral_descriptors();
        }
        if (ImGui::InputInt("descriptor no ", &descriptor_no));
        if (ImGui::MenuItem("Display descriptor"))
        {
            display_descriptor(m);
        }
        ImGui::EndMenu();

    }
    if (ImGui::MenuItem("Point matching with skeleton endpoints triangle area "))
    {
        trilateral_point_matching_with_skeleton_endpoints_w_HKS(m, skeleton, positive_desc, negative_desc, plane);
        is_draw_plane = true; 
    }
    if (ImGui::MenuItem("Point matching with skeleton endpoints with spin image  "))
    {
        trilateral_point_matching_with_skeleton_endpoints_SpinImage(m, skeleton, positive_desc,
            negative_desc, plane);
    }
    if (ImGui::MenuItem("Generate trilaterals using endpoints "))
    {
        trilateral_display_trilateral_from_skeleton_endpoints(m, positive_desc , negative_desc , skeleton , plane );
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
                    CoreType_conv_glm_raylib_vec3(skeleton.skeletonFormat[skeleton_adj_i[j]].point), RED);
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
        write_trilateral_descriptors(file_path_name);
        file_path_name = "";
    }
    drawFileDialog(file_path, file_path_name, ".dsc", is_reading_dsc);
    if (file_path_name != "")
    {
        read_trilateral_descriptors(file_path_name);

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


#pragma region trilateral save load 
static void save_trilateral_descriptors()
{
    is_writing_dsc = true; 
    
}
static void load_trilateral_descriptors()
{
    is_reading_dsc = true;
}
static void write_trilateral_descriptors(std::string filename)
{
    std::ofstream file;                // Create an ofstream object for file output

    // Open the file in write mode
    file.open(filename);

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return ;
    }

    // Write some data to the file
    //desc format 1 - 
    for (size_t i = 0; i < positive_desc.size(); i++)
    {
        file << " desc ";
        file << positive_desc[i].p1 << " " << positive_desc[i].p2 << " " << positive_desc[i].p3 << std::endl;
        for (size_t j = 0; j < positive_desc[i].visited_indices.size(); j++)
        {
            file << positive_desc[i].visited_indices[j] << " ";
        }
        file << std::endl;
    }
    file << "negative" << std::endl;
    for (size_t i = 0; i < negative_desc.size(); i++)
    {
        file << " desc ";
        file << negative_desc[i].p1 << " " << negative_desc[i].p2 << " " << negative_desc[i].p3 << std::endl;
        for (size_t j = 0; j < negative_desc[i].visited_indices.size(); j++)
        {
            file << negative_desc[i].visited_indices[j] << " ";
        }
        file << std::endl;
    }
    // Close the file
    file.close();

}

static void read_trilateral_descriptors(std::string filename)
{
    std::ifstream file(filename);                // Create an ofstream object for file output

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write some data to the file
    //desc format 1 - 
    int pos_size = positive_desc.size();
    std::vector<TrilateralDescriptor> descriptors;
    int negative_start_index = -1;
    bool is_pos_desc = true; 
    std::string line;

    TrilateralDescriptor desc; 
    while (std::getline(file, line))
    {

        if (line.find("negative") != std::string::npos)
        {
            is_pos_desc = false;
            negative_start_index = descriptors.size();
        }
        if (line.find("desc") != std::string::npos)
        {
            line = line.substr(5);
            std::stringstream ss(line);
            std::vector<int> nums;
            int num;
            while (ss >> num)
            {
                nums.push_back(num);
            }
            desc.p1 = nums[0];
            desc.p2 = nums[1];
            desc.p3 = nums[2];
        }
        else
        {
            std::stringstream ss(line);
            std::vector<unsigned int> visited;
            unsigned int num;
            while (ss >> num)
            {
                visited.push_back(num);
            }
            desc.visited_indices = visited;
            descriptors.push_back(desc);
        }
    }
    // Close the file
    file.close();
    
    for (size_t i = 0; i < negative_start_index; i++)
    {
        positive_desc.push_back(descriptors[i]);
    }  
    for (size_t i = negative_start_index; i < descriptors.size(); i++)
    {
        negative_desc.push_back(descriptors[i]);
    }
    return;
}

static void  display_descriptor(TrilateralMesh* m )
{
    TrilateralDescriptor desc; 
    if (descriptor_no < positive_desc.size())
    {
        desc = positive_desc[descriptor_no];
    }
    else
    {
        desc = positive_desc[descriptor_no - positive_desc.size()];
    }
    for (size_t i = 0; i < desc.visited_indices.size() ; i++)
    {
        int index = desc.visited_indices[i];

        m->raylib_mesh.colors[index * 4] = 255;
        m->raylib_mesh.colors[index * 4 + 1] = 0;
        m->raylib_mesh.colors[index * 4 + 2] = 0;
        m->raylib_mesh.colors[index * 4 + 3] = 255;
    }
    m->update_raylib_mesh();
}