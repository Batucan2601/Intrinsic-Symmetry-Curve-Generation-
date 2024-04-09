#include "glm/glm.hpp"
#include <map>
#include "Mesh.h"
#include "../Include/Skeleton.h"
struct NLateralDescriptor
{
	int N;
	std::vector<unsigned int> point_indices;
	std::vector<std::vector<double>> euclidian_distances;
	std::vector<std::vector<double>>  geodesic_distances;
	std::vector<std::vector<double>>  curvatures;
	double area;
	// extras
	std::vector<double> k_ring_areas;
	Mesh mesh;
	NLateralDescriptor(Mesh& mesh, const std::vector<unsigned int>& point_indices, int N);

	void get_euclidian_distances();
	void get_geodesic_distances();
	void get_curvatures();
	void get_k_ring_areas(int K);
	void get_ROI();

	
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

NLateralDescriptor generate_NLateralDescriptor(Mesh* m, const std::vector<unsigned int>& mesh_indices, const std::vector<bool>& parameter_checkbox
	, const std::vector<bool>& parameter_weights, const std::vector<std::string>& parameter_names);

std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_furthest_pairs(Mesh* m, std::vector<unsigned int>& indices,  NLateralParameters N_LATERAL_PARAMETERS);
std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_closest_pairs(Mesh* m, std::vector<unsigned int>& indices, NLateralParameters N_LATERAL_PARAMETERS);

std::vector <std::pair<unsigned int, unsigned int>> point_match_n_lateral_descriptors( Mesh* m ,const std::vector<NLateralDescriptor>& nlateral_vec_left, const std::vector<NLateralDescriptor>& n_lateral_vec_right
, NLateralParameters N_LATERAL_PARAMETERS);

void start_n_lateral_algorithm(Mesh* m , NLateralParameters N_LATERAL_PARAMETERS);
void start_n_lateral_algorithm_for_mesh(std::vector<SkeletonFormat>& mesh, NLateralParameters N_LATERAL_PARAMETERS);

void NLateral_parameters_calculate_maximums(Mesh* m, NLateralParameters& N_LATERAL_PARAMETERS , std::vector<unsigned int>&left , std::vector<unsigned int>&  right);

void start_n_lateral_algorithm_with_skeleton_end_points(Mesh* m, NLateralParameters& N_LATERAL_PARAMETERS,
	std::vector<unsigned int>& mesh_left_endpoints, std::vector<unsigned int>& mesh_right_endpoints);
