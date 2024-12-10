#include <vector>
#include "glm/glm.hpp"
#include "TrilateralMesh.h"
struct Curvature
{
	std::vector<glm::vec3> points; 
};

Curvature CurvatureGeneration_generate(TrilateralMesh* m);