#pragma once 
#include "glm/glm.hpp"
#include "TrilateralMesh.h"
#include "Histogram.h"
#include "raylib.h"

#include <algorithm>
#include <vector>



typedef struct {
	glm::vec3 point;
	glm::vec3 normal;
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;
	glm::vec3 p4;
}Plane;

struct TrilateralDescriptor
{
	double area; // ROI
	double total_length;
	double geodesic_lenght_1_2;//  geodesic length between 1 - 2
	double geodesic_lenght_1_3;//  length between 1 - 3
	double geodesic_lenght_2_3;//  length between 2 - 3
	double euclidian_lenght_1_2;//  euclidian length between 1 - 2
	double euclidian_lenght_1_3;//  euclidian length between 1 - 2
	double euclidian_lenght_2_3;//  euclidian length between 1 - 2
	double curvature_1_2;
	double curvature_1_3;
	double curvature_2_3;
	unsigned int p1; //point indices
	unsigned int p2;
	unsigned int p3;
	// extras
	float n_ring_area_p1;
	float n_ring_area_p2;
	float n_ring_area_p3;
	// histogram
	Histogram histogram; 

	TrilateralDescriptor();

};

enum  ComparisonMethod
{
	absoulute_dif,
	quadratic_dif,
};
enum PointStatus
{
	INSIDE,
	OUTSIDE,
	EDGE,
	MIDPOINT,
};

//plane stuff
TrilateralMesh generate_mesh_from_plane( Plane* plane , glm::vec3 * m  );
float get_point_status_from_plane(Plane* plane, glm::vec3* point);
glm::vec3 project_point_to_plane(Plane* plane, glm::vec3* point);
glm::vec3 symmetry_point_from_plane(Plane* plane, glm::vec3* point);
Plane generate_plane_from_two_vectors(const glm::vec3& vec1 , const glm::vec3& vec2);
Plane generate_plane_from_three_points(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);
Plane generate_plane_from_formula(const float& A , const float& B, const float& C, const float& D);


bool is_triangles_intersect(const glm::vec3& p1_1, const glm::vec3& p2_1, const glm::vec3& p3_1, const glm::vec3& p1_2, const glm::vec3& p2_2, const glm::vec3& p2_3);
float area_of_triangle_intersection(const glm::vec3& p1_1, const glm::vec3& p2_1, const glm::vec3&  p3_1 , const glm::vec3& p1_2 , const glm::vec3& p2_2 , const glm::vec3& p2_3  );
void get_coefficients_from_plane(const Plane& plane , float& A , float& B , float& C , float& D);
std::vector<int> getNumberFromString(std::string s);
Eigen::VectorXd stdVectorToEigenVectorXd(const std::vector<float>& std_vec);
float distancePointToLine(const glm::vec3& point, const glm::vec3& linePoint1, const glm::vec3& linePoint2);


float compute_triangle_area(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);

glm::vec3 CoreType_conv_raylib_glm_vec3(Vector3 vec3);
Vector3 CoreType_conv_glm_raylib_vec3(glm::vec3 vec3);
// permuation
float permutation_return_smallest_dif(Eigen::VectorXf vec1, Eigen::VectorXf vec2, int N);


//templates
template <typename T>
int CoreType_return_smallest_k(std::vector<T>& arr, int k_th_smallest)
{

	// Ensure n is within the valid range
	if (k_th_smallest < 1 || k_th_smallest > arr.size()) {
		std::cerr << "Invalid input: n is out of range." << std::endl;
		return -1; // Error code
	}

	// Use nth_element to find the n-th smallest element (0-based index)
	std::nth_element(arr.begin(), arr.begin() + (k_th_smallest - 1), arr.end());

	// Return the n-th smallest element (1-based index, so subtract 1)
	return arr[k_th_smallest - 1];
}

template <typename T>
void CoreType_sort_by_value(std::vector<std::pair<T, int>>& vec) {
	// Sorting the vector of pairs based on the first element (value)
	std::sort(vec.begin(), vec.end(), [](const std::pair<T, int>& a, const std::pair<T, int>& b) {
		return a.first < b.first;
		});
}


