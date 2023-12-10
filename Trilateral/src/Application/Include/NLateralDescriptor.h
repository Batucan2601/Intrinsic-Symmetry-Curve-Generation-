#include "glm/glm.hpp"
#include "Mesh.h"
struct NLateralDescriptor
{
	int N;
	std::vector<unsigned int> point_indices;
	std::vector<std::vector<double>> euclidian_distances;
	std::vector<std::vector<double>>  geodesic_distances;
	std::vector<std::vector<double>>  curvatures;
	// extras
	std::vector<double> k_ring_areas;
	Mesh* mesh;
	NLateralDescriptor(Mesh& mesh, const std::vector<unsigned int>& point_indices, int N);

	void get_euclidian_distances();
	void get_geodesic_distances();
	void get_curvatures();
	void get_k_ring_areas(int K);
};

NLateralDescriptor generate_NLateralDescriptor(Mesh* m, const std::vector<unsigned int>& mesh_indices, const std::vector<bool>& parameter_checkbox
	, const std::vector<bool>& parameter_weights, const std::vector<std::string>& parameter_names);