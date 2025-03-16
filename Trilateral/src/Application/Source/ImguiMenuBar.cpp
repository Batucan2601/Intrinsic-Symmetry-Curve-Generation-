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
#include "../Include/SCBDataset.h"
#include "ImGuiFileDialog.h"
#include "raymath.h"
#include <fstream>  // Include for file stream handling
#include "../Include/NLateralMapping.h"
#include "../Include/Geodesic.h"
#include "../Include/CurvatureGeneration.h"
#include "../Include/Voronoi.h"


static void imgui_menubar_save_mesh(TrilateralMesh* m );
static void dominant_symmetry_plane(TrilateralMesh* m);
static void skeleton_generation(TrilateralMesh* m);
static void trilateral_functions(TrilateralMesh* m);
static void dvorak_functions(TrilateralMesh* m);
static void distribution_functions(TrilateralMesh* m);
static void dataset(TrilateralMesh* m);
static void shape_diameter(TrilateralMesh* m);
static void SCAPE_dataset(TrilateralMesh* m);
static void KIDS_dataset(TrilateralMesh* m);
static void TOSCA_dataset(TrilateralMesh* m);
static void geodesic(TrilateralMesh* m);
static void curvature_creation(TrilateralMesh* m);

static void laplace_beltrami_operations(TrilateralMesh* m);
static void Nlateral_functions(TrilateralMesh* m);
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
static bool is_writing_ndsc = false;
static bool is_reading_ndsc = false;
static bool is_mesh_wires = false;
static bool is_draw_plane = false;
static bool is_draw_skeleton = false;
static bool is_draw_curvature = false;
static bool is_draw_curvature_index = false;
static bool is_draw_resemblance_pairs = true;
static bool is_draw_mesh = true;
static bool is_draw_normals = false; 
static bool is_draw_agd = false; 
static bool is_draw_desc = false; 
static bool is_draw_hist = false; 
static bool is_draw_curves = false; 
static float agd_sphere_radius = 1; 
static float line_thickness_3d = 0.01; 
static int descriptor_no = 0;
static int descriptor_no_voronoi1 = 0;
static int descriptor_no_voronoi2 = 0;
static Plane plane;
static Skeleton skeleton;
static BackBone backbone;
static Curvature curvature; 
static std::pair<Curvature, Curvature> curvature_front_and_back;
static std::vector<Curve> curves;
static int dvorak_no_of_significant_points = 0;
static int no_of_dist_points = 0;
static float dvorak_geodesic_dist_param = 0;
static float max_geodesic_param = 0;
static float longest_distance = 0;
static int mesh_index = 0;
std::vector<TrilateralDescriptor> positive_desc; 
std::vector<TrilateralDescriptor> negative_desc; 
std::vector<unsigned int> distributed_indices;
static Eigen::SparseMatrix<double> L; //laplacian 
static int N = 0;
std::vector<NLateralDescriptor> nlateral_descriptors;
std::pair<std::vector<NLateralDescriptor>, std::vector<NLateralDescriptor> > nlateral_descriptors_pos_neg;
std::vector<NodeAffinityParams> skeleton_params;
float hks_dif_param = 0.1;
float curv_param = 0.5;
float closeness_param = 0.2;
float skel_dist_param = 0.2;
float skel_depth_param = 5;
float ratio_dif_param= 0.2;
float proximity_param = 1;
float area_dif_param = 0.2;
float skel_point_dist_param = 0.2; 
float fuzzy_param = 0.1;
float min_geo_tau = 0.7;
int avg_geo_N_ring = 2;
float nlateral_tri_hist_param = 0.2;
float distance_to_mid_param = 0.7; 
float voronoi_param =  1 * 1e-2; 
float sdf_param = 0.2;
int hist_param = 5;
int min_agd_param  = 3;
std::vector<unsigned int> original_agd_vertices; 
std::vector<unsigned int> avg_dijk_indices;
int geodesic_point_no = 0;
std::vector<unsigned int> min_dijk_indices;
std::vector<unsigned int> voronoi_set;
int no_of_points_voronoi = 50;
int voronoi_no = 0;
Voronoi voronoi; 
//curvature
float quality_param = 0.7; 
int curv_point_index = 0;
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
            if (ImGui::BeginMenu("Laplace beltrami"))
            {
                laplace_beltrami_operations(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Mesh Drawing"))
            {
                mesh_drawing();
                ImGui::EndMenu();
            }
            if(ImGui::BeginMenu("NLateral"))
            {
                Nlateral_functions(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("dataset"))
            {
                dataset(m);
                SCAPE_dataset(m);
                TOSCA_dataset(m);
                KIDS_dataset(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Geodesic distances"))
            {
                geodesic(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Curvature creation "))
            {
                curvature_creation(m);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Shape diameter"))
            {
                shape_diameter(m);
                ImGui::EndMenu();
            }
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
    display_file_dialogs(m);
    if (is_draw_hist)
    {
        Nlateral_display_histogram(m, nlateral_descriptors, descriptor_no);
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
        //HKS_read_kernel_signature(m );
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
    if (ImGui::MenuItem("End point matching with  skeelton endpoint poins Optimal transform"))
    {
        trilateral_point_matching_with_skeleton_endpoints_and_OT(m,skeleton, positive_desc, negative_desc, plane, dvorak_no_of_significant_points, convergence_ratio);
        is_draw_plane = true;
    }
    if (ImGui::MenuItem("End point matching with Dvorak significant poins Optimal transform"))
    {
        trilateral_point_matching_with_gaussian_endpoints_and_OT(m, skeleton,positive_desc,negative_desc, plane, dvorak_no_of_significant_points, convergence_ratio);
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
static void Nlateral_functions(TrilateralMesh* m)
{
    ImGui::InputInt("N == ", &N);
    if (ImGui::BeginMenu("Nlateral params "))
    {
        ImGui::InputInt("gaussian points ", &dvorak_no_of_significant_points);
        ImGui::InputFloat("Sweep distance", &dvorak_geodesic_dist_param);
        ImGui::InputFloat("hks difference ", &hks_dif_param);
        ImGui::InputFloat("closeness param ", &closeness_param);
        ImGui::InputFloat("area param ", &area_dif_param);
        ImGui::InputFloat("min geo tau param", &min_geo_tau);
        ImGui::InputInt("N ring for geodesic", &avg_geo_N_ring);
        ImGui::InputFloat("distance to mid point", &distance_to_mid_param);
        ImGui::InputFloat("sdf param ", &sdf_param);
        ImGui::InputFloat("fuzzy param ", &fuzzy_param);
        ImGui::InputInt("hist no ", &hist_param);
        ImGui::InputInt("min agd no ", &min_agd_param);

        ImGui::EndMenu();
    }
    if (ImGui::MenuItem("generate descriptors with midpoint "))
    {
        nlateral_descriptors = NLateralMapping_generate_via_voronoi_midpoints(m, avg_dijk_indices, dvorak_geodesic_dist_param, min_geo_tau,fuzzy_param
        , distance_to_mid_param, hks_dif_param, closeness_param,sdf_param, hist_param , min_agd_param, longest_distance, original_agd_vertices,
            voronoi_param);
    }
    if (ImGui::MenuItem("get extra matches with curvature "))
    {
       // NLateralMapping_get_new_matchings(m, curvature, nlateral_descriptors , distance_to_mid_param , sdf_param , hks_dif_param , hist_param);
       NLateralMapping_get_new_matchings(m, curvature, nlateral_descriptors , distance_to_mid_param , hks_dif_param );
    }
    if (ImGui::MenuItem("Color distance from unmathced "))
    {
        NLateralMapping_get_unmathced_areas(m, curvature, 0.1 , true);
    }
    if (ImGui::MenuItem("Color midpoints of resemblances  "))
    {
        m->color_midpoints(ORANGE);
    }
    if (ImGui::BeginMenu("NLateral Descriptor"))
    {
        if (ImGui::MenuItem("Save nlateral descriptors"))
        {
            is_writing_ndsc = true;
        }
        if (ImGui::MenuItem("Load nlateral descriptors"))
        {
            is_reading_ndsc = true;
        }
        if (ImGui::InputInt("descriptor no ", &descriptor_no))
        {
            Nlateral_display_desc(m, nlateral_descriptors , descriptor_no);
            is_draw_desc = true;
        }
        if (ImGui::MenuItem("Display histogram "))
        {
            is_draw_hist = true;
        }
        if (ImGui::MenuItem("Display descriptor all"))
        {
            display_descriptor_all(m);
        }
        if (ImGui::MenuItem("Descriptor to backbone generation"))
        {
          skeleton_color_path_to_backbone(m, skeleton, backbone, nlateral_descriptors[descriptor_no].skeleton_index);
        }
        if (ImGui::MenuItem("Write matching points"))
        {
            Nlateral_write_matching_points(m , nlateral_descriptors);
        }
        if (ImGui::MenuItem("Read matching points"))
        {
            Nlateral_read_matching_points(m, nlateral_descriptors);
        }
        ImGui::EndMenu();

    }
    if (ImGui::BeginMenu("Voronoi test"))
    {
        if (ImGui::InputInt("descriptor 1 no ", &descriptor_no_voronoi1))
        {
        }
        if (ImGui::InputInt("descriptor 2 no ", &descriptor_no_voronoi2))
        {
        }
        if (ImGui::InputFloat("Voronoi", &voronoi_param))
        {

        }
        if (ImGui::MenuItem("show voronoi w/agd index"))
        {
            voronoi = Voronoi(m, avg_dijk_indices[descriptor_no_voronoi1],
                avg_dijk_indices[descriptor_no_voronoi2], voronoi_param);
            //voronoi_set = NLateral_show_voronoi(m, nlateral_descriptors[descriptor_no_voronoi1].indices[0], nlateral_descriptors[descriptor_no_voronoi2].indices[0]);
            voronoi.generate_voronoi_parts();
            voronoi.color();
        }
        if (ImGui::MenuItem("show voronoi w/index"))
        {
            voronoi = Voronoi(m, descriptor_no_voronoi1,
            descriptor_no_voronoi2, voronoi_param);
            //voronoi_set = NLateral_show_voronoi(m, nlateral_descriptors[descriptor_no_voronoi1].indices[0], nlateral_descriptors[descriptor_no_voronoi2].indices[0]);
            voronoi.generate_voronoi_parts();
            voronoi.color();
        }
        if (ImGui::MenuItem("Voronoi of midpoints"))
        {
            voronoi_set = NLateral_show_voronoi_midpoints(m);
        }
        if (ImGui::MenuItem("Get closest curve  with matches "))
        {
            voronoi = Voronoi_get_closest_voronoi(m, voronoi_param);
        }
        if (ImGui::MenuItem("Prune voronoi "))
        {
           Voronoi_prune_voronoi(m, voronoi , voronoi_param);
        }
        if (ImGui::MenuItem("Recalculate with better matches"))
        {
            voronoi = Voronoi_destroy_wrong_matches_and_recalculate(m,voronoi_param , voronoi );
        }
        if (ImGui::MenuItem("recalculate with distances "))
        {
            voronoi = Voronoi_check_pair_closeness_and_recalculate(m, voronoi_param,  distance_to_mid_param, voronoi);
        }
        if (ImGui::MenuItem("show connected voronoi set "))
        {
            NLateral_generate_voronoi_curve(m, voronoi_set, true);
        }
        ImGui::InputInt("no of points added to voronoi ", &no_of_points_voronoi);
        if (ImGui::MenuItem("Execute additive voronoi "))
        {
            Voronoi_algorithm_in_action(m, voronoi_param, hks_dif_param, sdf_param, distance_to_mid_param, hist_param, fuzzy_param, no_of_points_voronoi , avg_dijk_indices);
        }
        if (ImGui::MenuItem("Color every pairs of voronoi  "))
        {
            Voronoi_color_every_pairs_voronoi(m, voronoi_param);
        }
        if (ImGui::InputInt("pair no  ", &voronoi_no))
        {
            Voronoi_show_voronoi(m, voronoi_no, voronoi_param);

        }
        if (ImGui::MenuItem("Voronoi connect boundary "))
        {
            voronoi.connect_boundary();
        }
        ImGui::EndMenu();
    }
  
}
static void dataset(TrilateralMesh* m)
{
    if (ImGui::BeginMenu("select Dataset"))
    {
        if (ImGui::MenuItem("SCAPE"))
        {
            SCB_select_dataset(SCAPE);
        }
        if (ImGui::MenuItem("KIDS"))
        {
            SCB_select_dataset(KIDS);
        }
        if (ImGui::MenuItem("TOSCA"))
        {
            SCB_select_dataset(TOSCA);
        }
        ImGui::EndMenu();
    }
    ImGui::InputInt("Mesh Index", &mesh_index);
    if (ImGui::MenuItem("display mesh "))
    {
        SCB_select_mesh(*m, mesh_index, skeleton);
    }

}
static void shape_diameter(TrilateralMesh* m)
{
    if (ImGui::MenuItem("color shape diameter"))
    {
        ShapeDiameter_color(m);
    }
    if (ImGui::MenuItem("Draw normals"))
    {
        is_draw_normals = true;
    }
}
static void KIDS_dataset(TrilateralMesh* m)
{
    if (ImGui::MenuItem("read KIDS"))
    {
        SCB_read_KIDS();
    }
}
static void SCAPE_dataset(TrilateralMesh* m)
{
    if (ImGui::MenuItem("read SCAPE"))
    {
        SCB_read_SCAPE();
    }
    
}
static void TOSCA_dataset(TrilateralMesh* m)
{
    if (ImGui::MenuItem("read TOSCA"))
    {
        SCB_read_TOSCA();
    }
}

static void geodesic(TrilateralMesh* m)
{
    ImGui::InputInt("No of points" , &no_of_dist_points);
    ImGui::InputFloat("Sweep distance" , &dvorak_geodesic_dist_param);
    ImGui::InputFloat("max geodesic param" , &max_geodesic_param);
    if (ImGui::MenuItem("average geodesic Function"))
    {
        avg_dijk_indices = Geodesic_avg_dijkstra(m, no_of_dist_points, dvorak_geodesic_dist_param, avg_geo_N_ring, true);
    }
    if (ImGui::MenuItem("average geodesic Function modified"))
    {
        float biggest; 
        avg_dijk_indices = Geodesic_avg_dijkstra_modified(m, dvorak_geodesic_dist_param, avg_geo_N_ring, true, biggest);
    }
    if (ImGui::MenuItem("average geodesic Function modified with points "))
    {
        avg_dijk_indices = Geodesic_avg_dijkstra_modified_to_points(m, avg_dijk_indices,
        no_of_dist_points, dvorak_geodesic_dist_param, avg_geo_N_ring, true);
    }
    if (ImGui::MenuItem("minimum geodesic Function"))
    {
        avg_dijk_indices = Geodesic_min_dijkstra(m, avg_dijk_indices, dvorak_geodesic_dist_param, min_geo_tau, true);
    }
    if (ImGui::MenuItem("Biggest Geodesics "))
    {
        avg_dijk_indices = Geodesic_find_biggest_AGD(m, dvorak_geodesic_dist_param, max_geodesic_param);
    }
    if (ImGui::MenuItem("Write Sampled points "))
    {
        Geodesic_write_sampled_points(m, avg_dijk_indices);
    }
    if (ImGui::MenuItem("Read Sampled points "))
    {
        Geodesic_read_sampled_points(m, avg_dijk_indices);
    }
    if (ImGui::InputInt("Show the sampled point " , &geodesic_point_no))
    {

    }
    if (ImGui::MenuItem("FPS with midpoint sampling"))
    {
        float biggest; 
        unsigned int mid1, mid2; 
        avg_dijk_indices = midpoint_sampling(m,0.01 ,biggest,mid1 , mid2);
    }
    if(ImGui::MenuItem("Geodesic crete curve with midpoints "))
    {
        unsigned int m1, m2;
        Geodesic_generate_secondary_curve_w_midpoints(m, m1, m2);
    }
    if (ImGui::MenuItem("Generate multiple curves "))
    {
        unsigned int m1, m2;
        Geodesic_generate_multiple_secondary_curve(m, m1, m2);
    }
    static int path_p1;
    static int path_p2;
    ImGui::InputInt("point 1 ", &path_p1);
    ImGui::InputInt("point 2 ", &path_p2);
    if (ImGui::MenuItem("Show path "))
    {
        Geodesic_color_path(m , path_p1, path_p2);
    }
    if (ImGui::MenuItem("Color according to midpoints"))
    {
        Geodesic_color_according_to_midpoints(m);
    }
}

static void curvature_creation(TrilateralMesh* m )
{
    ImGui::InputFloat(" normal degree dif", &quality_param);
    if (ImGui::MenuItem(" Create curature from sym points "))
    {
        curvature = CurvatureGeneration_generate(m ,avg_dijk_indices,hks_dif_param, quality_param);
        is_draw_curvature = true;

    }
    if (ImGui::MenuItem(" Create curvature full with updates "))
    {
        curvature_front_and_back =  CurvatureGeneration_generate_full_curv(m, avg_dijk_indices, hks_dif_param, quality_param, distance_to_mid_param
        ,closeness_param,fuzzy_param ,longest_distance, original_agd_vertices);
        is_draw_curvature = true; 
    }
    if (ImGui::MenuItem(" Color each side "))
    {
        Curvature_color_sides(m, curvature);
    }
    if (ImGui::MenuItem(" curvature update once  "))
    {
        CurvatureGeneration_update(m, curvature_front_and_back.first, avg_dijk_indices, hks_dif_param, quality_param);
    }
    if (ImGui::MenuItem(" get new curve point with PCA "))
    {
        curvature = CurvatureGeneration_move_PCA_and_connect(m, curvature_front_and_back.first,
        curvature_front_and_back.second);
    }
    if (ImGui::MenuItem(" curvature update multiple times  "))
    {
        CurvatureGeneration_update_w_quality(m, curvature_front_and_back.first, avg_dijk_indices, hks_dif_param, quality_param);
    }
    ImGui::InputInt(" point index", &curv_point_index );
    if (ImGui::MenuItem(" Show point with index "))
    {
        is_draw_curvature_index = true; 
    }
    if (ImGui::MenuItem(" Generate minimum geodesic point "))
    {
        unsigned int p1; 
        unsigned int p2; 
        float biggest;
        Geodesic_mid_point_w_AGD(m, p1, p2,biggest);
    }
    if (ImGui::MenuItem("Laplacian smoothing"))
    {
        CurvatureGeneration_laplacian_smoothing(m, curvature, quality_param);
    }
    if (ImGui::MenuItem("Curve smoothing single "))
    {
        CurvatureGeneration_curve_smoothing(m, curvature, quality_param);
    }
    if (ImGui::MenuItem("Curve smooth all "))
    {
        while(CurvatureGeneration_curve_smoothing(m, curvature, quality_param));
    }
    if (ImGui::MenuItem("Add new matching"))
    {
        CurvatureGeneration_add_new_matching( m, curvature,avg_dijk_indices,curvature.midpoint_index,curvature.midpoint_inv_index, 
        quality_param ,hks_dif_param , distance_to_mid_param,sdf_param,fuzzy_param,longest_distance ,original_agd_vertices);
    }
    if (ImGui::MenuItem("Display matching points midpoints"))
    {
        curves = CurvatureGeneration_generate_curve_paths(m);
        is_draw_curves = true; 
    }
    if (ImGui::MenuItem("Curve prune front"))
    {
        curve_prune(m, curvature, true, false);
    }
    if (ImGui::MenuItem("Curve prune back"))
    {
        curve_prune(m, curvature, false, true);
    }
    if (ImGui::MenuItem("Generate curve with res pairs"))
    {
        curvature = CurvatureGeneration_generate_curve_w_sym_pairs(m);
        is_draw_curvature = true; 

    }
}

static std::pair<Eigen::VectorXd, Eigen::MatrixXd>  eigen_pairs;
int time_step = 0;
static std::vector<std::pair<int, float>> hks_pair; 
static void laplace_beltrami_operations(TrilateralMesh* m)
{
    if (ImGui::MenuItem(" Generate Laplace beltrami "))
    {
        L = laplace_beltrami(m);
    }
    if (ImGui::MenuItem(" EigenDecompose "))
    {
        eigen_pairs = laplace_beltrami_eigendecompose(L,30);
    }
    ImGui::InputInt("time step == ", &time_step);
    if (ImGui::MenuItem(" Heat kernel signature "))
    {
        std::vector<double> time_cons = { 0.1 ,0.2 , 0.3 , 0.4 , 1.0 , 2.0 , 3.0 , 4.0, 5.0 , 6.0,7.0,8.0, 9.0,10.0,15.0,20.0 ,200.0,2000.0};
        HKS_compute_kernel(m, eigen_pairs , time_cons , time_step);
    }
    if (ImGui::MenuItem(" Heat kernel signature on a descriptor "))
    {
        HKS_hks_on_descriptor(m, positive_desc[0]);

    }
    if (ImGui::MenuItem(" Heat kernel signature read HKS "))
    {
        //HKS_read_kernel_signature(m );
    }
    if (ImGui::MenuItem(" HKS significant points"))
    {
        hks_pair = HKS_extraction_significant_points(m , dvorak_no_of_significant_points);
    }
    if (ImGui::MenuItem(" HKS sweep distances"))
    {
        HKS_sweep_distance(m, hks_pair, dvorak_geodesic_dist_param);
    }
}



#pragma region draw section
static void draw_dom_sym();
static void draw_skeleton();
static void draw_resemblance_pairs(TrilateralMesh* m );
static void draw_mesh(TrilateralMesh* m, Shader& shader);
static void draw_curvature(TrilateralMesh* m );
static void draw_curves(TrilateralMesh* m);
static void draw_normals(TrilateralMesh* m);
static void draw_spheres(TrilateralMesh* m, float radius);
static void draw_descriptor(TrilateralMesh* m);


void draw_all_shader(TrilateralMesh* m , Shader& shader )
{
    draw_mesh(m , shader);
}
void draw_all(TrilateralMesh* m)
{
    draw_dom_sym();
    draw_skeleton();
    draw_resemblance_pairs(m);
    draw_curvature(m);
    draw_normals(m);
    draw_spheres(m, agd_sphere_radius);
    draw_curves(m);
    draw_descriptor(m);
}
static void draw_curves(TrilateralMesh* m)
{
    if (is_draw_curves)
    {
        for (size_t i = 0; i < curves.size(); i++)
        {
            
            for (size_t j = 0; j < curves[i].curve_path.size() - 1; j++)
            {
                int index_j = curves[i].curve_path[j];
                int index_j_1 = curves[i].curve_path[j + 1];
                /*DrawLine3D(CoreType_conv_glm_raylib_vec3(m->vertices[index_j]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[index_j_1]), YELLOW);*/
                DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[index_j]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[index_j_1]), line_thickness_3d, line_thickness_3d, 8, GREEN);
            }
            //draw midpoints
            unsigned int index = curves[i].curve_path[curves[i].curve_path.size()/2];
            DrawSphere(CoreType_conv_glm_raylib_vec3(m->vertices[index]), 0.01, RED);
        }
    }
    
}
static void draw_spheres(TrilateralMesh* m, float radius)
{
    if (is_draw_agd)
    { 
        for (size_t i = 0; i < avg_dijk_indices.size(); i++)
        {
            int index = avg_dijk_indices[i];
            if (i == geodesic_point_no)
            {
                DrawSphere(CoreType_conv_glm_raylib_vec3(m->vertices[index]), radius, BLUE);
            }
            else
            {
                DrawSphere(CoreType_conv_glm_raylib_vec3(m->vertices[index]), radius, RED);
            }
        }
        if (!(avg_dijk_indices.size() == 0))
        {
            return; 
        }
        //draw the results from nlateraldess
        for (size_t i = 0; i < nlateral_descriptors.size(); i++)
        {
            int index = nlateral_descriptors[i].indices[0];
            DrawSphere(CoreType_conv_glm_raylib_vec3(m->vertices[index]), radius, RED);
        }
    }
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
static void draw_normals(TrilateralMesh* m)
{
    if (is_draw_normals)
    {
        for (size_t i = 0; i < m->normals_display.size(); i+= 12 )
        {
            glm::vec3 p1(m->normals_display[i], m->normals_display[i +1], m->normals_display[i + 2]);
            glm::vec3 p2(m->normals_display[i + 6 ], m->normals_display[i + 7], m->normals_display[i + 8]);
            glm::vec3 dir = p2 - p1; 
            dir = dir  * 1e-2f;
            p2 = p1 + dir; 
            glm::vec3 mid = (p2 + p1) / 2.0f;
            DrawLine3D(CoreType_conv_glm_raylib_vec3(p1), CoreType_conv_glm_raylib_vec3(mid) 
                , RED);
            DrawLine3D(CoreType_conv_glm_raylib_vec3(mid), CoreType_conv_glm_raylib_vec3(p2)
                , GREEN);

        }
    }
}
static void draw_curvature(TrilateralMesh* m )
{
    if (is_draw_curvature)
    {
       


        for (size_t i = 0; i < curvature_front_and_back.first.paths.size(); i++)
        {
            curvature_front_and_back.first;
            for (size_t j = 0; j < curvature_front_and_back.first.paths[i].size() - 1; j++)
            {
                unsigned int p1 = curvature_front_and_back.first.paths[i][j];
                unsigned int p2 = curvature_front_and_back.first.paths[i][j + 1];
                DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[p1]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[p2]), line_thickness_3d, line_thickness_3d, 8, YELLOW);
            }
        }
        for (size_t i = 0; i < curvature_front_and_back.second.paths.size(); i++)
        {
            for (size_t j = 0; j < curvature_front_and_back.second.paths[i].size() - 1; j++)
            {
                unsigned int p1 = curvature_front_and_back.second.paths[i][j];
                unsigned int p2 = curvature_front_and_back.second.paths[i][j + 1];
                DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[p1]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[p2]), line_thickness_3d, line_thickness_3d, 8, YELLOW);
            }
        }

        for (size_t i = 0; i < curvature.paths.size(); i++)
        {
            for (size_t j = 0; j < curvature.paths[i].size() - 1; j++)
            {
                unsigned int p1 = curvature.paths[i][j];
                unsigned int p2 = curvature.paths[i][j + 1];
                DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[p1]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[p2]), line_thickness_3d, line_thickness_3d, 8, RED);
            }
        }
        for (size_t i = 0; i < curvature.curve_points.size(); i++)
        {
            int index = curvature.curve_points[i].mid_point;
            Color c = RED;
            c.r = c.r - (curvature.curve_points.size() * 255.0 / i);
            c.b = c.b + (curvature.curve_points.size() * 255.0 / i);
            DrawSphere(CoreType_conv_glm_raylib_vec3(m->vertices[index]), 0.01, c);
        }
 
    }
    if (is_draw_curvature_index)
    {
        int index = curvature.curve_points[curv_point_index].mid_point;
        int first = curvature.curve_points[curv_point_index].correspondence.first;
        int second = curvature.curve_points[curv_point_index].correspondence.second;
        DrawSphere(CoreType_conv_glm_raylib_vec3(m->vertices[index]), 0.01, WHITE);
        DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[first]), CoreType_conv_glm_raylib_vec3(m->vertices[second]),
            line_thickness_3d, line_thickness_3d, 8,WHITE );
    }

}
static void draw_resemblance_pairs(TrilateralMesh* m )
{
    if (is_draw_resemblance_pairs)
    {
        for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
        {
            int index1 = m->calculated_symmetry_pairs[i].first;
            int index2 = m->calculated_symmetry_pairs[i].second;
            if (descriptor_no == i)
            {
                DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[index1]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[index2]), line_thickness_3d, line_thickness_3d, 8, RED);
            }
            else
            {
                DrawCylinderEx(CoreType_conv_glm_raylib_vec3(m->vertices[index1]),
                    CoreType_conv_glm_raylib_vec3(m->vertices[index2]), line_thickness_3d, line_thickness_3d, 8, BLUE);
            }
        }
    }
   
}
static void draw_mesh(TrilateralMesh* m , Shader& shader )
{
    Material mat = LoadMaterialDefault();
    mat.shader = shader;

    if (is_draw_mesh)
    {
        if (!is_mesh_wires)
        {
            DrawMesh(m->raylib_mesh, mat, MatrixIdentity());
        }
        else
        {
            DrawMeshWires(m->raylib_mesh, Vector3{ 0,0,0 }, 1, BLACK);
        }
    }
    free(mat.maps); 
    
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
    drawFileDialog(file_path, file_path_name, ".dsc", is_writing_ndsc);
    if (file_path_name != "")
    {
        NLateralDescriptor_write(file_path_name, m, nlateral_descriptors);
        file_path_name = "";
    }
    drawFileDialog(file_path, file_path_name, ".dsc", is_reading_ndsc);
    if (file_path_name != "")
    {
        NLateralDescriptor_read(file_path_name, m , nlateral_descriptors);
        file_path_name = "";
    }
}
static void draw_descriptor(TrilateralMesh* m)
{
    if (is_draw_desc)
    {
        NLateralDescriptor* desc = &nlateral_descriptors[descriptor_no];
        
        Vector3 p1 = CoreType_conv_glm_raylib_vec3( m->vertices[desc->indices[0]]);
        Vector3 p2 = CoreType_conv_glm_raylib_vec3(m->vertices[desc->indices[1]]);
        Vector3 p3 = CoreType_conv_glm_raylib_vec3(m->vertices[desc->indices[2]]);
        DrawSphere(p1,agd_sphere_radius , WHITE );
        DrawSphere(p2 , agd_sphere_radius,BLUE);
        DrawSphere(p3 , agd_sphere_radius, RED);
    }
}
static void mesh_drawing()
{
    ImGui::InputFloat("line thickness"  , &line_thickness_3d );
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
    ImGui::InputFloat(" Radius for AGD Spheres "  , &agd_sphere_radius);
    if (is_draw_agd)
    {
        if (ImGui::MenuItem("Disable AGD Drawing"))
        {
            is_draw_agd = !is_draw_agd;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable AGD Drawing "))
        {
            is_draw_agd = !is_draw_agd;
        }
    }
    if (is_draw_resemblance_pairs)
    {
        if (ImGui::MenuItem("Disable resemblance Drawing"))
        {
            is_draw_resemblance_pairs = !is_draw_resemblance_pairs;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable resemblance Drawing "))
        {
            is_draw_resemblance_pairs = !is_draw_resemblance_pairs;
        }
    }
    if (is_draw_curves)
    {
        if (ImGui::MenuItem("Disable curve Drawing"))
        {
            is_draw_curves = !is_draw_curves;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable curve Drawing "))
        {
            is_draw_curves = !is_draw_curves;
        }
    }
    if (is_draw_curvature)
    {
        if (ImGui::MenuItem("Disable curvature Drawing"))
        {
            is_draw_curvature = !is_draw_curvature;
        }
    }
    else
    {
        if (ImGui::MenuItem("Enable curvature Drawing "))
        {
            is_draw_curvature = !is_draw_curvature;
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
        if (ImGui::MenuItem("Uniform Point Sampling "))
        {
            distributed_indices = uniform_point_sampling(m, no_of_dist_points, true);
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

    m->raylib_mesh.colors[desc.p1 * 4] = 255;
    m->raylib_mesh.colors[desc.p1 * 4 + 1] = 255;
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

