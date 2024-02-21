#include "../Include/Mesh.h"

typedef struct
{
	glm::vec3 origin;
	glm::vec3 direction;

}Ray;
bool ray_triangle_intersection(Ray& ray, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
