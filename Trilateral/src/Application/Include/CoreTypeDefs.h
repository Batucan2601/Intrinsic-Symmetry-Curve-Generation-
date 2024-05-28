#pragma once 
#include "glm/glm.hpp"
#include "Mesh.h"
#include <eigen/Eigen/Dense>

typedef struct {
	glm::vec3 point;
	glm::vec3 normal;
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
Mesh generate_mesh_from_plane( Plane* plane , glm::vec3 * m  );
float get_point_status_from_plane(Plane* plane, glm::vec3* point);
glm::vec3 project_point_to_plane(Plane* plane, glm::vec3* point);
glm::vec3 symmetry_point_from_plane(Plane* plane, glm::vec3* point);
Plane generate_plane_from_two_vectors(const glm::vec3& vec1 , const glm::vec3& vec2);
Plane generate_plane_from_three_points(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);
bool is_triangles_intersect(const glm::vec3& p1_1, const glm::vec3& p2_1, const glm::vec3& p3_1, const glm::vec3& p1_2, const glm::vec3& p2_2, const glm::vec3& p2_3);
float area_of_triangle_intersection(const glm::vec3& p1_1, const glm::vec3& p2_1, const glm::vec3&  p3_1 , const glm::vec3& p1_2 , const glm::vec3& p2_2 , const glm::vec3& p2_3  );
void get_coefficients_from_plane(const Plane& plane , float& A , float& B , float& C , float& D);
std::vector<int> getNumberFromString(std::string s);
Eigen::VectorXd stdVectorToEigenVectorXd(const std::vector<float>& std_vec);
