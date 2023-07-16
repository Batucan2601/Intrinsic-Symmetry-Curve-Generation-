#pragma once
#include "../Include/MeshFactory.h"
#include "../Include/TrilateralMap.h"
#include "../Include/Sampling.h"
#include "../Include/Laplace-Beltrami.h"
#include <src/Application/Include/DominantSymmetry.h>

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

int no_of_sampling_fps = 10; 
Plane plane; 
std::vector<std::vector<int>> symmetry_paired_points;
void imgui_mesh_window(int& selected_mesh, MeshFactory& m_factory )
{

    ImGui::Begin("Input window");
    

    ImGui::InputInt("mesh no  ", &selected_mesh);

    ImGui::InputInt("vertex 1  ", &point_1_index);
    ImGui::InputInt("vertex 2  ", &point_2_index);
    ImGui::InputInt("vertex 3  ", &point_3_index);
    ImGui::InputInt("histogram_partition  ", &partition_no);
    //ImGui::InputFloat("tau  ", &tau);
    ImGui::InputFloat("scale  ", &scale);
    /* if (ImGui::Button("bilateral map "))
    {
        if_isocurve_selected = false;
        compute_bilateral_map(m, point_1_index, point_2_index, tau, histogram);
    }
    ImGui::InputInt("vertex 1  ", &point_1_index);
    if (ImGui::Button("bilateral map according to point above "))
    {
        if_isocurve_selected = false;
        compute_bilateral_map_according_to_point(m, point_1_index, tau, histogram);
    }
    
    if (ImGui::Button("Isocurves with point above "))
    {
        iso_curve_distances = compute_iso_curves(m, point_1_index);
        if_isocurve_selected = true; 
    }
    if (ImGui::Button("compute isocurve histogram "))
    {
        iso_curve_histogram = compute_curve_distances(m , iso_curve_distances, point_1_index );
    }
    if (ImGui::Button("Compute all geodesic styles "))
    {
        compute_all(m);
    }
    if (ImGui::Button("Compute from histogram "))
    {
        find_point_from_histogram(m, iso_curve_histogram  );
    }
    if (ImGui::Button("Load model 2  "))
    {
        Mesh off((char*)"faust-registrations//tr_reg_000.off");
        m = off;
        buffer_mesh(m);
        buffer_vertices_with_triangles(m);
        if_isocurve_selected = false;
    }
    
    
    if (ImGui::Button("Quadify mesh "))
    {
        is_mesh_quadified = true; 
        quadify_mesh(m);
        buffer_quadified_mesh(m);
    }
    if (ImGui::Button("Catmull mesh "))
    {
        is_catmull_clark = true;
        catmull_clark_new(m);
        buffer_after_clark(m);
    }
    if (ImGui::Button("Color after Clark "))
    {
        color_weight_according_to_distance(m);
    }
    if (ImGui::Button("Calculate total area quad "))
    {
        catmull_calculate_triangle_area(m);
    }
    if (ImGui::Button("sqrt(3)"))
    {
        sqrt_3(m);
        buffer_sqrt31(m);
        is_sqrt_3 = true; 
    }
    if (ImGui::Button("color sqrt3 algorithm "))
    {
        color_sqrt3(m);
        sqrt3_calculate_triangle_area(m);
    }
    */
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

    ImGui::InputInt("mesh no for cross test 1  ", &selected_mesh_for_points1);
    ImGui::InputInt("mesh no for cross test 2  ", &selected_mesh_for_points2);
    if (ImGui::Button("Eslestir "))
    {
        match_points_from2_mesh(m_factory, selected_mesh_for_points1, selected_mesh_for_points2,  partition_no);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Mock eslestir"))
    {
        lines =  match_points_from2_mesh_mock(m_factory, selected_mesh_for_points1, selected_mesh_for_points2, partition_no);
        is_draw_lines_activated = true; 
        is_draw_lines_activated_once = true; 
    }
    if (ImGui::Button("Generate Trilateral Descriptor "))
    {
       TrilateralDescriptor trilateral_descriptor =  generate_trilateral_descriptor(m_factory, selected_mesh, point_1_index, point_2_index, point_3_index , false );
       std::cout << "area " << trilateral_descriptor.area << " " << trilateral_descriptor.lenght_1_2 << " " << trilateral_descriptor.lenght_1_3 << "  " << trilateral_descriptor.lenght_2_3 << " " << trilateral_descriptor.curvature_1_2 << "  " << trilateral_descriptor.curvature_1_3 << "  " << trilateral_descriptor.lenght_2_3 << std::endl;
    }
    if (ImGui::Button("Brute Force Trilateral Symmetry "))
    {
        brute_force_symmetry_extraction(m_factory, selected_mesh);
        m_factory.remove_all();
        m_factory.add_all();
    }
    ImGui::InputInt("no of samples  ", &no_of_sampling_fps);
    ImGui::SameLine();
    if (ImGui::Button("FPS sampling"))
    {
        furthest_point_sampling(&m_factory.mesh_vec[selected_mesh] , no_of_sampling_fps);
        m_factory.remove_all();
        m_factory.add_all();
    }
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
        point_matching_with_dominant_symmetry_plane(m_factory, selected_mesh, &plane, no_of_sampling_fps);
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