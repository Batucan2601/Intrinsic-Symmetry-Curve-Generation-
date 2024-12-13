#pragma once
#include "glm/glm.hpp"
#include <map>
#include "TrilateralMesh.h"
#include "Histogram.h"
#include "../Include/SkeletonTypes.h"
#include "glm/gtc/matrix_transform.hpp"
#include <glm/gtc/type_ptr.hpp>
#include "eigen/Eigen/dense"
#include "eigen/Eigen/sparse"
struct NLateralDescriptor
{
	std::vector<unsigned int> point_indices;
	std::vector<std::vector<double>> euclidian_distances;
	std::vector<std::vector<double>>  geodesic_distances;
	std::vector<std::vector<double>>  curvatures;
	// extras
	std::vector<double> k_ring_areas;
	TrilateralMesh mesh;
	NLateralDescriptor(TrilateralMesh& mesh, const std::vector<unsigned int>& point_indices, int N);
	NLateralDescriptor();

	void get_euclidian_distances();
	void get_geodesic_distances();
	void get_curvatures();
	void get_k_ring_areas(int K);
	void get_ROI();

	// new implementation
	void create_histogram_area(TrilateralMesh* m , int N);
	void create_histogram_HKS(TrilateralMesh* m , int N);
	void create_histogram_SDF(TrilateralMesh* m, int N);
	int N; //this was present
	std::vector<unsigned int> indices; // fist one is the origin point 
	std::vector<std::vector<std::vector<int>>> paths;  // from indices[0] to others
	std::vector<float> distances;  // distances to others
	std::vector<unsigned int> vertices_inside;
	std::vector<unsigned int> triangles_inside;
	Eigen::VectorXd weight;
	std::vector<float> skel_dist_mid;
	unsigned int skeleton_index;
	std::vector<unsigned int> depth; 
	float n_ring_area;
	double area;
	double skel_point_dist;
	double paths_ratio; 
	Histogram histogram;
	float max_distance; 
};

struct NLateralDescriptorRestrictions
{
	int no_of_points;
	
	float sweep_distance;
	bool is_sweep_distance;
	float hks_dif_param;
	bool is_hks_dif; 
	float skel_dist_param;
	bool is_skel_dist;
	float n_ring_param;
	bool is_n_ring; 
	float area_dif_param;
	bool is_area_dif; 
	float skel_point_dist_param; 
	bool is_skel_point_dist;
	float paths_dif_param;
	bool is_paths_dif; 
	float min_geo_tau;
	bool is_min_geo_tau; 
	
};
struct NLateralParameters
{
	//constants
	static const int NO_OF_PARAMETERS  =  9;
	static const int K_RING_POS =  5;
	static const int N_LATERAL_CONSTRUCTION_METHOD_NO =  2;

	//variables that can be changed
	int no_of_N_lateral_pairs =  100;
	std::string current_n_lateral_construction_method;
	int N;
	static const int N_RING_NO = 1 ; // how many layer of breadth first search.

	std::vector<std::string> parameter_names;
	std::vector<float> parameter_weights;
	std::map<std::string, float> parameter_maximums; //holds maximums of each descriptor to norrmalzie it
	std::vector<std::string> n_lateral_construction_methods;

	// logic of imgui
	std::vector<bool> parameter_checkbox;

	NLateralParameters();
};

static NLateralParameters N_LATERAL_PARAMETERS; 



NLateralDescriptor generate_NLateralDescriptor(TrilateralMesh* m, const std::vector<unsigned int>& mesh_indices, const std::vector<bool>& parameter_checkbox
	, const std::vector<bool>& parameter_weights, const std::vector<std::string>& parameter_names);

std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_furthest_pairs(TrilateralMesh* m, std::vector<unsigned int>& indices,  NLateralParameters N_LATERAL_PARAMETERS);
std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_closest_pairs(TrilateralMesh* m, std::vector<unsigned int>& indices, NLateralParameters N_LATERAL_PARAMETERS);

std::vector <std::pair<unsigned int, unsigned int>> point_match_n_lateral_descriptors( TrilateralMesh* m ,const std::vector<NLateralDescriptor>& nlateral_vec_left, const std::vector<NLateralDescriptor>& n_lateral_vec_right
, NLateralParameters N_LATERAL_PARAMETERS);

//void start_n_lateral_algorithm(MeshFactory& mesh_fac, int selected_mesh, NLateralParameters N_LATERAL_PARAMETERS);
//void start_n_lateral_algorithm_for_mesh(std::vector<SkeletonFormat>& mesh, NLateralParameters N_LATERAL_PARAMETERS);

void NLateral_parameters_calculate_maximums(TrilateralMesh* m, NLateralParameters& N_LATERAL_PARAMETERS , std::vector<unsigned int>&left , std::vector<unsigned int>&  right);

void start_n_lateral_algorithm_with_skeleton_end_points(TrilateralMesh* m, NLateralParameters& N_LATERAL_PARAMETERS,
	std::vector<unsigned int>& mesh_left_endpoints, std::vector<unsigned int>& mesh_right_endpoints);
NLateralDescriptor NLateral_generate_descriptor(TrilateralMesh* m, const std::vector<unsigned int>& mesh_indices);

std::vector<NLateralDescriptor> NLateral_generate_closest_points(TrilateralMesh* m,  std::vector<unsigned int>& indices, 
int N, int histogram_size );

std::vector<unsigned int> Nlateral_check_vertices_visited(TrilateralMesh* m, NLateralDescriptor& desc);
std::vector<unsigned int> Nlateral_check_triangles_visited(TrilateralMesh* m, NLateralDescriptor& desc);
float Nlateral_get_maximum_dist(TrilateralMesh* m, NLateralDescriptor& desc);
void NLateral_compute_skel_point_dist(TrilateralMesh* m, Skeleton& skel, NLateralDescriptor& desc);
void Nlateral_display_desc(TrilateralMesh* m, std::vector<NLateralDescriptor>& descs, int index);
void Nlateral_display_desc(TrilateralMesh* m, std::pair<std::vector<NLateralDescriptor>, std::vector<NLateralDescriptor>>& descs, int index);

void NLateralDescriptor_write(std::string filename, TrilateralMesh* m, std::vector<NLateralDescriptor>& desc);
void NLateralDescriptor_read(std::string filename, TrilateralMesh* m, std::vector<NLateralDescriptor>& desc);

std::vector<unsigned int> NLateral_sweepDistance(TrilateralMesh* m , std::vector<unsigned int> indices, float sweep_distance );
bool NLateral_check_path_lengths(NLateralDescriptor& desc1, NLateralDescriptor& desc2, float similarity);

double NLateral_get_paths_ratio(TrilateralMesh* m,NLateralDescriptor& desc);
unsigned int NLateral_get_closest_index_to_midpoint(TrilateralMesh* m, std::vector<unsigned int>& points);

bool NLateral_compare_HKS(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float hks_perc , std::ofstream& file);
bool NLateral_compare_skeldist_mid(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float skel_percentage, float maximum_skel, std::ofstream& file);
bool NLateral_compare_Nring(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float maximum_n_ring, std::ofstream& file);
bool NLateral_compare_area(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float maximum_area, float area_percentage, std::ofstream& file);
bool NLatera_compare_depth(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, int depth_dif_param, std::ofstream& file);
bool NLateral_compare_trilateral_with_midpoint(TrilateralMesh* m, unsigned int p1, unsigned int p2, unsigned int p_middle, float dissimilarity
,std::ofstream& file );
bool NLateral_compare_distance_to_midpoint(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, unsigned int midpoint_index
,float distance_to_mid_param, std::ofstream& file);
bool NLateral_compare_SDF(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float maximum_sdf ,
float sdf_param , std::ofstream& file );

bool Nlateral_check_endpoint(TrilateralMesh* m, Skeleton& skel, NLateralDescriptor& desc1, NLateralDescriptor& desc2);
bool NLateral_compare_position_to_midpoint(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, unsigned int midpoint_index,
	float distances_from_mid, float distances_between_desc , std::ofstream& file);

