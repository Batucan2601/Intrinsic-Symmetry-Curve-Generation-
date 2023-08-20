#pragma once 
#include "glm/glm.hpp"
#include "Mesh.h"
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
};

enum ComparisonMethod
{
	absoulute_dif = 0,
	quadratic_dif = 1,
};
Mesh generate_mesh_from_plane( Plane* plane , glm::vec3 * m  );
float get_point_status_from_plane(Plane* plane, glm::vec3* point);