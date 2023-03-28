#pragma once 
#include "glm/glm.hpp"
#include "Mesh.h"
typedef struct {
	glm::vec3 point;
	glm::vec3 normal;
}Plane;

Mesh generate_mesh_from_plane( Plane* plane , glm::vec3 * m  );
float get_point_status_from_plane(Plane* plane, glm::vec3* point);