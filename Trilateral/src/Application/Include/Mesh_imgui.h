#pragma once
#include "../Include/MeshFactory.h"
#include "../Include/TrilateralMap.h"
#include "../Include/Sampling.h"
#include "../Include/Laplace-Beltrami.h"
#include <src/Application/Include/DominantSymmetry.h>
#include "../Include/CoreTypeDefs.h"
#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
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
std::vector<std::pair<unsigned int, unsigned int >> calculated_symmetry_pairs; 
Mesh m1, m2;
std::vector<int> m1_map_indices, m2_map_indices;
std::string KIDS_text_file_name;
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
        
        //is_visited = trialteral_ROI(m_factory , selected_mesh, point_1_index, point_2_index, point_3_index, partition_no , is_visited_interior);
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
    if (ImGui::Button("separate mesh with dominant symmetry plane "))
    {
        generate_two_separate_mesh_using_dominant_symmetry_plane(plane , &m_factory.mesh_vec[selected_mesh] , &m1 , &m2 , &m1_map_indices , &m2_map_indices ) ;
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("match two separated meshes with FPS  "))
    {
        //match_two_meshes_with_fps(&m_factory.mesh_vec[selected_mesh] ,&m1, &m2, &m1_map_indices, &m2_map_indices, no_of_points);
        trilateralDescVector = match_two_meshes_with_fps(&m_factory.mesh_vec[selected_mesh], &plane, no_of_points);
        is_trilateral_generated = true;
    }
    /*if (ImGui::Button("point matching using dominant symmetry plane "))
    {
        point_matching_with_dominant_symmetry_plane(m_factory, selected_mesh, &plane, no_of_points);
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
        if (ImGui::Selectable((const char*)"Random Symmetry Pairs", &isSelected))
        {
            selectedIndices = random_symmetry_indices_sampling(&m_factory.mesh_vec[selected_mesh], no_of_points);
            curtrilateralItem = "Random Symmetry Pairs";
            is_trilateral_generated = true;
        }
        m_factory.remove_all();
        m_factory.add_all();
        ImGui::EndCombo();
    }
    if (ImGui::Button("trilateral drawing "))
    {
        trilateral_map_drawing_using_three_points(m_factory, selected_mesh, point_1_index, point_2_index, point_3_index);
    }
    if (ImGui::Button("point matching using AGD"))
    {
        trilateralDescVector = get_trilateral_points_using_closest_pairs(m_factory, selected_mesh, selectedIndices);
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
    if (ImGui::Button("Read symmetry values"))
    {
        read_symmetry_format((char*)"../../Trilateral/Mesh/off/sym.txt", &m_factory.mesh_vec[selected_mesh]);
    }
    if (ImGui::Button("Symmetry Plane using Isomap"))
    {
        plane = generate_isomap_embedding(&m_factory.mesh_vec[selected_mesh] ,false , 1);
        Mesh plane_mesh = generate_mesh_from_plane( &plane , &plane.point);
        m_factory.add_mesh(plane_mesh);
        
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("generate symmetry plane with classical MDS "))
    {
        generate_symmetry_plane_dividing_classical_MDS(&m_factory.mesh_vec[selected_mesh]);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("generate symmetry plane with landmark MDS "))
    {
        Mesh  landmark_mesh = compute_landmark_MDS(&m_factory.mesh_vec[selected_mesh] , 3 );
        m_factory.mesh_vec.clear();
        m_factory.add_mesh(landmark_mesh);

        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("generate  trilateral descriptors from symmetry plane with landmark MDS "))
    {
        //Plane plane = trilateral_symmetry_with_landmark_MDS_with_plane(&m_factory.mesh_vec[selected_mesh], 3);
        float error_percentage;
        Plane plane = trilateral_symmetry_with_landmark_MDS_with_plane(&m_factory.mesh_vec[selected_mesh], 3 , 100 , 100 , error_percentage);
        Mesh plane_mesh = generate_mesh_from_plane(&plane, &plane.point);
        m_factory.add_mesh(plane_mesh);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::InputText("textfile_name_for_kids_dataset" , &KIDS_text_file_name ))
    {

    }
    if (ImGui::Button("generate  trilateral descriptors w sym plane w KIDS dataset "))
    {
        // total of 15 + 15 meshes 
        std::vector<Mesh> mesh_vector; 
        for (size_t i = 0; i < 15 +15; i++)
        {
            std::string path("../../Trilateral/Mesh/off/");
            std::string isometry_batch_no( "000" + std::to_string((i / 15) + 1));
            std::string isometry_no(std::to_string(i % 15 + 1) );
            // read the meshes.
            path = path + isometry_batch_no + ".isometry." + isometry_no + ".off";
            Mesh m((char*)path.c_str() );
            //read the symmetry format
            read_symmetry_format((char*)"../../Trilateral/Mesh/off/sym.txt", &m);
            mesh_vector.push_back(m);
        }
        create_trilateral_sym_w_landmarl_with_planes(mesh_vector, 3, 100, 100, "../../Results/" + KIDS_text_file_name);
    }
    if (ImGui::Button("get n ring "))
    {
        get_N_ring_area(&m_factory.mesh_vec[selected_mesh], 100, 1);
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
    
    Mesh m = m_factory.mesh_vec[selected_mesh];

    // get each point  p_i and 2 others with minimal geoesic distance for p_i
    if (ImGui::Button("Generate Trilateral Point Pairs Using Minimum Distance"))
    {
        trilateralDescVector = get_trilateral_points_using_closest_pairs(m_factory, selected_mesh, selectedIndices);
    }
    if (ImGui::Button("Point Matching Using trilateral Weights "))
    {
        calculated_symmetry_pairs = point_match_trilateral_weights(&m, trilateralDescVector , trilateralCurvatureWeight ,trilateralGeodesicWeight , trilateralAreaWeight );
    }
    if (ImGui::Button("Display Accuracy"))
    {
        display_accuracy(&m, calculated_symmetry_pairs);
    }
    ImGui::End();
}


float N; 
#define BOOL_CONTROL_SIZE 2 
bool n_1_bool_control[BOOL_CONTROL_SIZE];
std::string method_name;
//parameters
#define NUMBER_OF_PARAMETERS 9 
bool  parameter_checkbox[NUMBER_OF_PARAMETERS];
float parameter_weights[NUMBER_OF_PARAMETERS];
std::string parameter_names[NUMBER_OF_PARAMETERS];
void imgui_N_Lateral_Parameters(const int& selected_mesh, MeshFactory& m_factory)
{
    float N;

    //init arrays
    for (size_t i = 0; i < BOOL_CONTROL_SIZE; i++)
    {
        n_1_bool_control[i] = false; 
    }
    for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
    {
        parameter_checkbox[i] = false;
    }
    for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
    {
        parameter_weights[i] = 1.0;
    }
    ImGui::Begin("N lateral Params ");
    ImGui::InputFloat("N parameter in N lateral:", &N);

// point selection algorithm
    ImGui::Text(" Select which way to fetch other n-1 lateral points for each.");
    if (ImGui::Checkbox("closest points", &n_1_bool_control[0]))
    {
        n_1_bool_control[0] != n_1_bool_control[0];
        method_name = "closest points";
    }
    if (ImGui::Checkbox("furthest points", &n_1_bool_control[1]))
    {
        n_1_bool_control[1] != n_1_bool_control[1];
        method_name = "furthest points";
    }
// parameters
    parameter_names[0] = "area";
    parameter_names[1] = "euclidian distance";
    parameter_names[2] = "geodesic distance";
    parameter_names[3] = "curvature";
    parameter_names[4] = "Heat Kernel Signature";
    parameter_names[5] = "X";
    parameter_names[6] = "X";
    parameter_names[7] = "X";
    parameter_names[8] = "X";

    ImGui::Text("Select required paramters and give them weights ");
    
    for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
    {
        if (ImGui::Checkbox(parameter_names[i].c_str(), &parameter_checkbox[i])) {}
        ImGui::SameLine();
        if (ImGui::InputFloat("weight = ", &parameter_weights[i])) {}
    }
    
    Mesh* m = &m_factory.mesh_vec[selected_mesh];
    if (ImGui::Button("Start algorithm for current mesh"))
    {
        //start_n_lateral_algorithm(m, N, method_name, parameter_checkbox, parameter_weights, parameter_names);
    }
    if (ImGui::Button("Start algorithm for all of the dataset"))
    {

    }
    ImGui::End();

}