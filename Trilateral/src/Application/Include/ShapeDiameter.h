#include "../Include/Mesh.h"

const int CONE_RADIUS = 2;
const int CONE_DISTANCE = 1;
const int NUMBER_OF_RAYS = 10;


void ShapeDiameter_calculate(Mesh* mesh, std::vector<unsigned int> indices, std::vector<float>& shape_diameter);