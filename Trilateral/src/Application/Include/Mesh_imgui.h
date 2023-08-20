#pragma once
#include "../Include/MeshFactory.h"
#include "../Include/TrilateralMap.h"
#include "../Include/Sampling.h"
#include "../Include/Laplace-Beltrami.h"
#include <src/Application/Include/DominantSymmetry.h>
#include "../Include/CoreTypeDefs.h"
bool if_bilateral_map = true; 
bool if_isocurve_selected = false;
bool if_bilateral_map_selected = true; 
int point_1_index = 0 ;
int point_2_index = 100;
int point_3_index = 200;
float tau = 50.0f;
int partition_no = 10; 
std::vector<int> iso_curve_distances;
std::vector<float> iso_curve_histogram;
float scale = 1.0f;


bool is_mesh_quadified = false;
bool is_catmull_clark = false;
bool is_sqrt_3 = false; 
bool is_trilateral = false; 

int selected_mesh = 0; 

int* is_visited; // checks if vertex is visited 

bool is_visited_interior = false; 

// lines
int selected_mesh_for_points1 = 0;
int selected_mesh_for_points2 = 1;
bool is_draw_lines_activated = false;
bool is_draw_lines_activated_once = false; // indicator for one time buffering of points
bool is_polygon_filled = true ; 
std::vector<float> line_points; //line points will be global 

bool activate_histogram = false; 

std::vector<float> histogram; 

std::vector<float> lines; 

//int no_of_sampling_fps = 10; 
//int no_of_agd_points = 10; 
int no_of_points = 10;

Plane plane; 
std::vector<std::vector<int>> symmetry_paired_points;
std::vector<unsigned int> selectedIndices;
std::vector<unsigned int> indices;
std::vector<TrilateralDescriptor> trilateralDescVector; 

static std::string curtrilateralItem;
bool is_trilateral_generated = false;

float trilateralCurvatureWeight = 1; 
float trilateralGeodesicWeight = 1;
float trilateralAreaWeight = 1;

//spectral embedding
std::vector<glm::vec3> embed_vertices; 
void imgui_mesh_window(int& selected_mesh, MeshFactory& m_factory )
{

    ImGui::Begin("Input window");
    

    ImGui::InputInt("mesh no  ", &selected_mesh);

    ImGui::InputInt("vertex 1  ", &point_1_index);
    ImGui::InputInt("vertex 2  ", &point_2_index);
    ImGui::InputInt("vertex 3  ", &point_3_index);
    ImGui::InputInt("histogram_partition  ", &partition_no);
    ImGui::InputFloat("scale  ", &scale);

    if (is_polygon_filled)
    {
        if (ImGui::Button("Enable polygon mode "))
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            is_polygon_filled = false; 
        }
    }
    else
    {
        if (ImGui::Button("Disable polygon mode "))
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            is_polygon_filled = true;
        }
    }
    if (ImGui::Button("Trilateral  "))
    {
        //trilateral_map(m_factory , selected_mesh, point_1_index, point_2_index, point_3_index);
        
        is_visited = trialteral_ROI(m_factory , selected_mesh, point_1_index, point_2_index, point_3_index, partition_no , is_visited_interior);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Histogram  "))
    {
        //trilateral_map(m_factory , selected_mesh, point_1_index, point_2_index, point_3_index);

        histogram = histogramROi(m_factory, selected_mesh, point_1_index, point_2_index, point_3_index, partition_no,  is_visited , is_visited_interior);
        m_factory.remove_all();
        m_factory.add_all();
        activate_histogram = true; 
    }
    ImGui::InputInt("no of samples  ", &no_of_points);
    ImGui::SameLine();
    if (ImGui::Button("Laplacian generation "))
    {
        generate_L(&m_factory.mesh_vec[selected_mesh]);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("dominant symmetry plane  "))
    {
        plane =  generate_dominant_symmetry_plane(selected_mesh , m_factory);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("point matching using dominant symmetry plane "))
    {
        point_matching_with_dominant_symmetry_plane(m_factory, selected_mesh, &plane, no_of_points);
        m_factory.remove_all();
        m_factory.add_all();
    }
    //ImGui::InputInt("agd point no : ", &no_of_agd_points);
    /*if (ImGui::Button("Average AGD function"))
    {
        AGDIndices = AverageGeodesicFunction(m_factory , selected_mesh , no_of_agd_points);
        m_factory.remove_all();
        m_factory.add_all();
    }
    ImGui::InputInt("mgd point no : ", &no_of_mgd_points);
    if (ImGui::Button("Minimum Geodesic function"))
    {
        MGDIndices = minimumGeodesicFunction  (m_factory, selected_mesh, no_of_mgd_points, AGDIndices);
        m_factory.remove_all();
        m_factory.add_all();
    }*/
    if (ImGui::BeginCombo("trilateral generation using", curtrilateralItem.c_str() )) // The second parameter is the label previewed before opening the combo.
    {
        bool isSelected = false; 
        if (ImGui::Selectable( (const char*)"AGD", &isSelected))
        {
            selectedIndices= AverageGeodesicFunction(m_factory, selected_mesh, no_of_points);
            //trilateralDescVector = match_points_using_trilateral_decriptor(m_factory, selected_mesh, AGDIndices);
            curtrilateralItem = "AGD";
            is_trilateral_generated = true; 
        }
        if (ImGui::Selectable((const char*)"MGD", &isSelected))
        {
            selectedIndices = minimumGeodesicFunction(m_factory, selected_mesh, no_of_points, selectedIndices);
            curtrilateralItem = "MGD";
            is_trilateral_generated = true;

        }
        if( ImGui::Selectable((const char*)"FPS", &isSelected))
        {
            selectedIndices = furthest_point_sampling(&m_factory.mesh_vec[selected_mesh], no_of_points);
            curtrilateralItem = "FPS";
            is_trilateral_generated = true;
        }
        m_factory.remove_all();
        m_factory.add_all();
        ImGui::EndCombo();
    }
    if (ImGui::Button("point matching using AGD"))
    {
        trilateralDescVector = get_trilateral_points_using_min_distance(m_factory, selected_mesh, selectedIndices);
    }
    if (ImGui::Button("Reset Points"))
    {
        reset_points(m_factory, selected_mesh);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if(ImGui::Button("Generate spectral embedding"))
    {
        embed_vertices = generate_spectral_embedding(m_factory, selected_mesh , selectedIndices);
        m_factory.remove_all();
        m_factory.add_all();
    }
    ImGui::End();
        
}

//display the attributes of selected mesh 
void imgui_selected_mesh_properties_window(const int& selected_mesh , MeshFactory& m_factory)
{
    ImGui::Begin("Mesh properties ");
    
    ImGui::InputFloat(" X ", &m_factory.mesh_vec[selected_mesh].model_mat[3][0]);
    ImGui::InputFloat(" Y ", &m_factory.mesh_vec[selected_mesh].model_mat[3][1]);
    ImGui::InputFloat(" Z ", &m_factory.mesh_vec[selected_mesh].model_mat[3][2]);

    ImGui::End();
}

void imgui_trilateralConfiguration(const int& selected_mesh, MeshFactory& m_factory)
{
    ImGui::Begin("Trialteral Conf");
    ImGui::InputFloat("CurvatureWeight :", &trilateralCurvatureWeight);
    ImGui::InputFloat("GeodesicWeight :", &trilateralGeodesicWeight);
    ImGui::InputFloat("AreaWeight :", &trilateralAreaWeight);
    
    if (ImGui::Button("Generate Trilateral Point Pairs Using Minimum Distance"))
    {
        trilateralDescVector = get_trilateral_points_using_min_distance(m_factory, selected_mesh, selectedIndices);
    }
    if (ImGui::Button("Generate Trilateral Features Using Points "))
    {
        for (size_t i = 0; i < trilateralDescVector.size() ; i++)
        {
            trilateralDescVector[i] = generate_trilateral_descriptor(m_factory, (int&)selected_mesh, trilateralDescVector[i].p1, trilateralDescVector[i].p2, trilateralDescVector[i].p3, false);
        }
    }
    if (ImGui::Button("Point Matching Using trilateral Weights "))
    {
        point_match_trilateral_weights(m_factory, (int&)selected_mesh , trilateralDescVector , trilateralCurvatureWeight ,trilateralGeodesicWeight , trilateralAreaWeight , quadratic_dif);
    }
    ImGui::End();
}
