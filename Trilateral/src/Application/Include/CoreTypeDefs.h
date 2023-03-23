#include "glm/glm.hpp"
#include "MeshFactory.h"

typedef struct {
	glm::vec3 point;
	glm::vec3 normal;
}Plane;

Mesh generate_mesh_from_plane( Plane * plane , glm::vec3 * m  );