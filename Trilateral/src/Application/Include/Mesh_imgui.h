#pragma once
#include "../Include/MeshFactory.h"
#include "../Include/TrilateralMap.h"
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
#include "../Include/SkeletalNLateral.h"
#include "../Include/ROI.h"

bool if_bilateral_map = true;
bool if_isocurve_selected = false;
bool if_bilateral_map_selected = true;
int point_1_index = 0;
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

std::vector<int> is_visited; // checks if vertex is visited 

bool is_visited_interior = false;

// lines
int selected_mesh_for_points1 = 0;
int selected_mesh_for_points2 = 1;
bool is_draw_lines_activated = false;
bool is_draw_lines_activated_once = false; // indicator for one time buffering of points
bool is_polygon_filled = true;
bool is_normals_shown = false;
bool is_skeletalNLateral_created = false;

std::vector<float> line_points; //line points will be global 

bool activate_histogram = false;

std::vector<float> histogram;

std::vector<float> lines;

//int no_of_sampling_fps = 10; 
//int no_of_agd_points = 10; 
int no_of_points = 10;
int no_of_hist_division = 10;
int no_of_sym_plane_iterations = 2; 
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
TrilateralMesh m1, m2;
std::vector<int> m1_map_indices, m2_map_indices;
std::string KIDS_text_file_name;

//fuzziness
float fuzziness = 1;




static std::vector<float> bounding_box;
static std::map<std::string, glm::vec3 > key_points;
static std::vector<float> shape_diameter;
static std::vector<float> skeleton_lines;
static Skeleton skeleton;
static BackBone best_backbone;
static std::vector<std::pair<unsigned int, unsigned int >> skeleton_best_end_point_pairs;
static std::vector<std::pair<int, int>>  skeletalNLateral_end_point_results;

static std::vector<TrilateralDescriptor> desc_r;
static std::vector<TrilateralDescriptor> desc_l;
void imgui_KIDS_skeleton(const int& selected_mesh, MeshFactory& m_factory)
{
    if (ImGui::Button("Generate Bounding Box For mesh"))
    {
        bounding_box = generate_bounding_box("0001.isometry.12.txt");
    }
    if (ImGui::Button("Generate Skeleton Keypoints for mesh"))
    {
        //key_points = generate_skeleton_keypoints("0001.isometry.1.txt");
        skeleton_lines = generate_skeleton_lines("0001.isometry.12.txt");
    }
    if (ImGui::Button("Match skeleton"))
    {
        //match_skeleton_keypoints(m_factory, &m_factory.mesh_vec[selected_mesh], bounding_box, key_points);
        match_skeleton_lines(m_factory, &m_factory.mesh_vec[selected_mesh], bounding_box, skeleton_lines);
        m_factory.remove_all();
        m_factory.add_all();
    }
    ImGui::LabelText("Skeleton generation with Cohen-Or's method", "Value");
    if (ImGui::Button("Generate Skeleton"))
    {
       // skeleton = skeleton_read_swc_file(m_factory, "0001.isometry.8.swc");
    }
    if (ImGui::Button("Generate Backbone"))
    {
        std::vector<unsigned int> right_mesh_end_points;
        std::vector<unsigned int> left_mesh_end_points;
        skeleton_generate_backbone(m_factory, skeleton, selected_mesh, best_backbone, right_mesh_end_points, left_mesh_end_points);
        //skeleton = skeleton_read_swc_file(m_factory, "0001.isometry.8.swc");
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("embed mesh to 2D "))
    {
        try {
            embed_mesh_endpoints_to_2d(m_factory.mesh_vec[selected_mesh], skeleton, N_LATERAL_PARAMETERS);
        }
        catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << std::endl;
        }
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Start N-Lateral algorithm for skeleton"))
    {
        std::vector<unsigned int> right_mesh_end_points;
        std::vector<unsigned int> left_mesh_end_points;
        skeleton_generate_backbone(m_factory, skeleton, selected_mesh, best_backbone, right_mesh_end_points, left_mesh_end_points);
        skeleton_get_N_Lateral_points(m_factory, skeleton, selected_mesh, best_backbone, skeleton_best_end_point_pairs
            , right_mesh_end_points, left_mesh_end_points);

        start_n_lateral_algorithm_with_skeleton_end_points(&m_factory.mesh_vec[selected_mesh], N_LATERAL_PARAMETERS,
            left_mesh_end_points, right_mesh_end_points);

        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Match End Points With N-Lateral"))
    {
        skeletalNLateral_end_point_results = SkeletalNLateral_compare_endpoints_with_SkeletalNlateral(skeleton,
            &m_factory.mesh_vec[selected_mesh], N_LATERAL_PARAMETERS.N, shape_diameter);
        SkeletalNLateral_generate_buffer(skeleton, skeletalNLateral_end_point_results);
        SkeletalNLateral_buffer();
        is_skeletalNLateral_created = true;
    }
    if (ImGui::Button("Shape diameter function "))
    {
        std::vector<unsigned int> indices;
        for (size_t i = 0; i < skeleton.endPoints.size(); i++)
        {
            int mesh_index1 = mesh_get_closest_index(&m_factory.mesh_vec[selected_mesh]
                , skeleton.skeletonFormat[i].point);
            indices.push_back(mesh_index1);

        }
        ShapeDiameter_calculate(&m_factory.mesh_vec[selected_mesh], indices, shape_diameter);

        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Best backbone generation using mid point"))
    {
        std::vector<unsigned int> right_mesh_end_points;
        std::vector<unsigned int> left_mesh_end_points;
        skeleton_generate_backbone_w_midpoint(m_factory, skeleton, selected_mesh, best_backbone, right_mesh_end_points, left_mesh_end_points);
        //skeleton = skeleton_read_swc_file(m_factory, "0001.isometry.8.swc");
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Generate Spectral Skeleton"))
    {
       // skeleton = skeleton_read_swc_file(m_factory, "skeletonSpectral.swc");
    }
}


void imgui_mesh_window(int& selected_mesh, MeshFactory& m_factory)
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
    if (ImGui::Button("Show normals"))
    {
        is_normals_shown = !is_normals_shown;
    }
    if (ImGui::Button("Trilateral  "))
    {
        //trilateral_map(m_factory , selected_mesh, point_1_index, point_2_index, point_3_index);

        is_visited = ROI_trilateral(&m_factory.mesh_vec[selected_mesh], point_1_index, point_2_index, point_3_index, partition_no, true);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Histogram  "))
    {
        //trilateral_map(m_factory , selected_mesh, point_1_index, point_2_index, point_3_index);

        //histogram = histogramROi(m_factory, selected_mesh, point_1_index, point_2_index, point_3_index, partition_no,  is_visited );
        m_factory.remove_all();
        m_factory.add_all();
        activate_histogram = true;
    }
    ImGui::InputInt("no of samples  ", &no_of_points);

    ImGui::InputInt("no of sym plane iterations", &no_of_sym_plane_iterations);
    ImGui::SameLine();
    if (ImGui::Button("dominant symmetry plane  "))
    {

    }
    if (ImGui::Button("separate mesh with dominant symmetry plane "))
    {
        dom_sym_generate_two_separate_mesh_using_dominant_symmetry_plane(plane, &m_factory.mesh_vec[selected_mesh], &m1, &m2, &m1_map_indices, &m2_map_indices);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("match two separated meshes with FPS  "))
    {
        //match_two_meshes_with_fps(&m_factory.mesh_vec[selected_mesh] ,&m1, &m2, &m1_map_indices, &m2_map_indices, no_of_points);
        trilateralDescVector = match_two_meshes_with_fps(&m_factory.mesh_vec[selected_mesh], &plane, no_of_points);
        is_trilateral_generated = true;
    }
    if (ImGui::BeginCombo("trilateral generation using", curtrilateralItem.c_str())) // The second parameter is the label previewed before opening the combo.
    {
        bool isSelected = false;
        if (ImGui::Selectable((const char*)"AGD", &isSelected))
        {
            selectedIndices = AverageGeodesicFunction(m_factory, selected_mesh, no_of_points);
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
        if (ImGui::Selectable((const char*)"FPS", &isSelected))
        {
            selectedIndices = furthest_point_sampling(&m_factory.mesh_vec[selected_mesh], no_of_points, true);
            curtrilateralItem = "FPS";
            is_trilateral_generated = true;
            m_factory.remove_all();
            m_factory.add_all();
        }
        if (ImGui::Selectable((const char*)"Random Symmetry Pairs", &isSelected))
        {
            selectedIndices = random_symmetry_indices_sampling(&m_factory.mesh_vec[selected_mesh], no_of_points);
            curtrilateralItem = "Random Symmetry Pairs";
            is_trilateral_generated = true;
        }

        ImGui::EndCombo();
    }
    if (ImGui::Button("trilateral drawing "))
    {
        trilateral_map_drawing_using_three_points(m_factory, selected_mesh, point_1_index, point_2_index, point_3_index);
        m_factory.remove_all();
        m_factory.add_all();
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
    if (ImGui::Button("Generate spectral embedding"))
    {
        embed_vertices = generate_spectral_embedding(m_factory, selected_mesh, selectedIndices);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Read symmetry values"))
    {
        read_symmetry_format((char*)"../../Trilateral/TrilateralMesh/off/sym.txt", &m_factory.mesh_vec[selected_mesh]);
    }
    if (ImGui::Button("Symmetry Plane using Isomap"))
    {
        plane = generate_isomap_embedding(&m_factory.mesh_vec[selected_mesh], false, 1);
        TrilateralMesh plane_mesh = generate_mesh_from_plane(&plane, &plane.point);
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
        TrilateralMesh  landmark_mesh = compute_landmark_MDS(&m_factory.mesh_vec[selected_mesh], 3);
        m_factory.mesh_vec.clear();
        m_factory.add_mesh(landmark_mesh);

        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("generate  trilateral descriptors from symmetry plane with landmark MDS "))
    {
        //Plane plane = trilateral_symmetry_with_landmark_MDS_with_plane(&m_factory.mesh_vec[selected_mesh], 3);
        float error_percentage;
        Plane plane = trilateral_symmetry_with_landmark_MDS_with_plane(&m_factory.mesh_vec[selected_mesh], 3, 100, 100, error_percentage);
        TrilateralMesh plane_mesh = generate_mesh_from_plane(&plane, &plane.point);
        m_factory.add_mesh(plane_mesh);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::InputText("textfile_name_for_kids_dataset", &KIDS_text_file_name))
    {

    }
    if (ImGui::Button("generate  trilateral descriptors w sym plane w KIDS dataset "))
    {
        // total of 15 + 15 meshes 
        std::vector<TrilateralMesh> mesh_vector;
        for (size_t i = 0; i < 15 + 15; i++)
        {
            std::string path("../../Trilateral/TrilateralMesh/off/");
            std::string isometry_batch_no("000" + std::to_string((i / 15) + 1));
            std::string isometry_no(std::to_string(i % 15 + 1));
            // read the meshes.
            path = path + isometry_batch_no + ".isometry." + isometry_no + ".off";
            TrilateralMesh m((char*)path.c_str());
            //read the symmetry format
            read_symmetry_format((char*)"../../Trilateral/TrilateralMesh/off/sym.txt", &m);
            mesh_vector.push_back(m);
        }
        create_trilateral_sym_w_landmarl_with_planes(mesh_vector, 3, 100, 100, "../../Results/" + KIDS_text_file_name);
    }
    if (ImGui::Button("get n ring "))
    {
        get_N_ring_area(&m_factory.mesh_vec[selected_mesh], 100, 1);
    }
    ImGui::InputInt(" histogram sampling ", &no_of_hist_division);
    if (ImGui::Button("FPS and histogram matching "))
    {
        trilateral_FPS_histogram_matching(m_factory, selected_mesh, no_of_points, no_of_hist_division, true);
        m_factory.remove_all();
        m_factory.add_all();
        //void trilateral_FPS_histogram_matching(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no)

    }
    if (ImGui::Button("FPS and histogram matching w/ spin images "))
    {
        trilateral_FPS_histogram_matching_w_spin_image(m_factory, selected_mesh, no_of_points, no_of_hist_division);
        m_factory.remove_all();
        m_factory.add_all();
        //void trilateral_FPS_histogram_matching(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no)

    }
    if (ImGui::Button("FPS and histogram matching w/ spin images w/ Multidimensional scaling "))
    {
        trilateral_FPS_histogram_matching_w_spin_image_MDS(m_factory, selected_mesh, no_of_points, no_of_hist_division);
        m_factory.remove_all();
        m_factory.add_all();
        //void trilateral_FPS_histogram_matching(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no)

    }
    if (ImGui::Button("mesh and histogram matching w/ principal components "))
    {
        trilateral_FPS_histogram_matching_w_principal_comp(m_factory, selected_mesh, no_of_points, no_of_hist_division);
        m_factory.remove_all();
        m_factory.add_all();
    }
    ImGui::InputFloat("fuzziness", &fuzziness);
    if (ImGui::Button("Calculate fuzzy geodesic areas with trilateral points "))
    {
        trilateral_fuzzyGeodesic(m_factory, selected_mesh, point_1_index, point_2_index, point_3_index, fuzziness);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("FPS matching w/ fuzzy geodesic "))
    {
        trilateral_FPS_matching_w_fuzzy_geodesic(m_factory, selected_mesh, no_of_points, fuzziness, true);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("skeleton endpoint matching w/ fuzzy geodesic "))
    {
        trilateral_w_skeleton_endpoints(m_factory, selected_mesh, fuzziness, skeleton, false);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("FPS matching using dominant sym plane and LMDS "))
    {
        std::string empty = " ";
        trilateral_self_matching_with_dominant_sym(m_factory, selected_mesh, no_of_points, empty , no_of_sym_plane_iterations, false  );
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("point matching with dominant plane and skeleton"))
    {
       // trilateral_point_matching_with_skeleton_endpoints(m_factory, selected_mesh, skeleton );
        m_factory.remove_all();
        m_factory.add_all();

    }
    if (ImGui::Button("point matching with dominant plane and skeleton and HKS "))
    {
     //   trilateral_point_matching_with_skeleton_endpoints_w_HKS(m_factory, selected_mesh, skeleton , desc_l , desc_r , plane);
        m_factory.remove_all();
        m_factory.add_all();

    }
    if (ImGui::Button(""))
    {
     //   trilateral_point_matching_with_skeleton_endpoints_w_HKS(m_factory, selected_mesh, skeleton, desc_l , desc_r , plane);
        m_factory.remove_all();
        m_factory.add_all();

    }
    if (ImGui::Button("Heat Kernal Signature"))
    {
        //HKS_extract_kernel_signature(&m_factory.mesh_vec[selected_mesh]);
        HKS_read_kernel_signature(&m_factory.mesh_vec[selected_mesh]);
    }
    ImGui::End();

}

//display the attributes of selected mesh 
void imgui_selected_mesh_properties_window(const int& selected_mesh, MeshFactory& m_factory)
{
    ImGui::Begin("TrilateralMesh properties ");

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

    TrilateralMesh m = m_factory.mesh_vec[selected_mesh];

    // get each point  p_i and 2 others with minimal geoesic distance for p_i
    if (ImGui::Button("Generate Trilateral Point Pairs Using Minimum Distance"))
    {
        trilateralDescVector = get_trilateral_points_using_closest_pairs(m_factory, selected_mesh, selectedIndices);
    }
    if (ImGui::Button("Point Matching Using trilateral Weights "))
    {
        calculated_symmetry_pairs = point_match_trilateral_weights(&m, trilateralDescVector, trilateralCurvatureWeight, trilateralGeodesicWeight, trilateralAreaWeight);
    }
    if (ImGui::Button("Display Accuracy"))
    {
        display_accuracy(&m, calculated_symmetry_pairs);
    }
    ImGui::End();
}


void imgui_N_Lateral_Parameters(const int& selected_mesh, MeshFactory& m_factory)
{

    ImGui::Begin("N lateral Params ");
    ImGui::InputInt("N parameter in N lateral:", &N_LATERAL_PARAMETERS.N);
    ImGui::InputInt("no of n latearl points ", &N_LATERAL_PARAMETERS.no_of_N_lateral_pairs);

    // point selection algorithm
    ImGui::Text(" Select which way to fetch other n-1 lateral points for each.");
    static const char* current_point_Selection_method = N_LATERAL_PARAMETERS.n_lateral_construction_methods[0].c_str();
    if (ImGui::BeginCombo("combobox", current_point_Selection_method))
    {
        for (int n = 0; n < N_LATERAL_PARAMETERS.N_LATERAL_CONSTRUCTION_METHOD_NO; n++)
        {
            bool is_selected = (current_point_Selection_method == N_LATERAL_PARAMETERS.n_lateral_construction_methods[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(N_LATERAL_PARAMETERS.n_lateral_construction_methods[n].c_str(), is_selected))
                current_point_Selection_method = N_LATERAL_PARAMETERS.n_lateral_construction_methods[n].c_str();
            if (is_selected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

    N_LATERAL_PARAMETERS.current_n_lateral_construction_method = current_point_Selection_method;

    ImGui::Text("Select required paramters and give them weights ");

    for (size_t i = 0; i < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; i++)
    {
        bool temp = N_LATERAL_PARAMETERS.parameter_checkbox[i];
        if (ImGui::Checkbox(N_LATERAL_PARAMETERS.parameter_names[i].c_str(), &temp)) {}
        ImGui::SameLine();
        if (ImGui::InputFloat("weight = ", &N_LATERAL_PARAMETERS.parameter_weights[i])) {}

        if (i == N_LATERAL_PARAMETERS.K_RING_POS)
        {
            ImGui::SameLine();
            int k_no;
            if (ImGui::InputInt("K in K ring = ", &k_no))
            {
                N_LATERAL_PARAMETERS.parameter_names[N_LATERAL_PARAMETERS.K_RING_POS] += std::to_string(k_no);
            }

        }

        N_LATERAL_PARAMETERS.parameter_checkbox[i] = temp;
    }
    if (ImGui::Button("Read symmetry values"))
    {
        read_symmetry_format((char*)"../../Trilateral/Mesh/off/sym.txt", &m_factory.mesh_vec[selected_mesh]);
    }
    if (ImGui::Button("Start algorithm for current mesh"))
    {
        TrilateralMesh* m = &m_factory.mesh_vec[selected_mesh];


        start_n_lateral_algorithm(m_factory , selected_mesh, N_LATERAL_PARAMETERS);
        m_factory.remove_all();
        m_factory.add_all();
    }
    if (ImGui::Button("Start algorithm for all of the dataset"))
    {

        // total of 15 + 15 meshes 
        std::vector<TrilateralMesh> mesh_vector;
        for (size_t i = 0; i < 15 + 15; i++)
        {
            std::string path("../../Trilateral/TrilateralMesh/off/");
            std::string isometry_batch_no("000" + std::to_string((i / 15) + 1));
            std::string isometry_no(std::to_string(i % 15 + 1));
            // read the meshes.
            path = path + isometry_batch_no + ".isometry." + isometry_no + ".off";
            TrilateralMesh m((char*)path.c_str());
            //read the symmetry format
            read_symmetry_format((char*)"../../Trilateral/TrilateralMesh/off/sym.txt", &m);
            mesh_vector.push_back(m);
            //start_n_lateral_algorithm(&m, N_LATERAL_PARAMETERS);
        }

    }

    ImGui::End();

}
static float camera_pos_x = 0.0f;
static float camera_pos_y = 0.0f;
static float camera_pos_z = 0.0f;

//testing
static int dijkstra_index;
static std::vector<int> vertex_list;
static std::vector<float> dijkstra_distances;
void imgui_debug_layer(int& selected_mesh, MeshFactory& mesh_fac, glm::vec3& cameraPos, glm::vec3 cameraDir, glm::vec3 up)
{


    ImGui::Begin("DEBUG layer");
    ImGui::Text(" camerapos is %f %f %f ", cameraPos.x, cameraPos.y, cameraPos.z);
    ImGui::Text(" cameraDir is %f %f %f ", cameraDir.x, cameraDir.y, cameraDir.z);
    ImGui::Text(" cameraUp is %f %f %f ", up.x, up.y, up.z);
    ImGui::InputFloat("camera X ", &camera_pos_x);
    ImGui::InputFloat("camera Y ", &camera_pos_y);
    ImGui::InputFloat("camera Z ", &camera_pos_z);
    if (ImGui::Button("Teleport to coordinates above "))
    {
        cameraPos = glm::vec3(camera_pos_x, camera_pos_y, camera_pos_z);
    }
    ImGui::InputInt("End point Index", &dijkstra_index);
    if (ImGui::Button("Calculate skeletal dijkstra"))
    {
        //get the n'th end point
        SkeletonEndPoint end_point_index = skeleton.endPoints[dijkstra_index];
        skeleton_calculate_dijkstra(skeleton, end_point_index.index,
            vertex_list, dijkstra_distances);
        //get smallest and color it
        int smallest_index = -1;
        float smallest_dist = INFINITY;
        for (size_t i = 0; i < skeleton.endPoints.size(); i++)
        {
            SkeletonEndPoint end_point_index_temp = skeleton.endPoints[i];
            if (dijkstra_distances[end_point_index_temp.index] < smallest_dist && end_point_index_temp.index != end_point_index.index)
            {
                smallest_index = end_point_index_temp.index;
                smallest_dist = dijkstra_distances[end_point_index_temp.index];
            }
        }

        mesh_fac.mesh_skeleton_vec.skeleton_points[end_point_index.index * 6 + 3] = 0.0f;
        mesh_fac.mesh_skeleton_vec.skeleton_points[end_point_index.index * 6 + 3 + 1] = 255.0f;
        mesh_fac.mesh_skeleton_vec.skeleton_points[end_point_index.index * 6 + 3 + 2] = 0.0f;

        mesh_fac.mesh_skeleton_vec.skeleton_points[smallest_index * 6 + 3] = 0.0f;
        mesh_fac.mesh_skeleton_vec.skeleton_points[smallest_index * 6 + 3 + 1] = 255.0f;
        mesh_fac.mesh_skeleton_vec.skeleton_points[smallest_index * 6 + 3 + 2] = 0.0f;

        skeleton_buffer(mesh_fac);

    }
    ImGui::InputInt("Save TrilateralMesh", &dijkstra_index);
    ImGui::End();



}


