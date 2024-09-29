#pragma once 
#include "../Include/TrilateralMesh.h"

typedef struct
{
	glm::vec3 origin;
	glm::vec3 direction;

}TrilateralRay;
bool ray_triangle_intersection(TrilateralRay& ray, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3& hit_point);
