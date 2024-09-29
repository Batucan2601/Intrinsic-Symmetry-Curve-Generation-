#pragma once 
#include "../Include/TrilateralMesh.h"

const int CONE_RADIUS = 3;
const int CONE_RADIUS_MULTIPLIER = 100; //for random number generation

const int CONE_HEIGHT = 20;
const int NUMBER_OF_RAYS = 30;


void ShapeDiameter_calculate(TrilateralMesh* mesh, std::vector<unsigned int> indices, std::vector<float>& shape_diameter);