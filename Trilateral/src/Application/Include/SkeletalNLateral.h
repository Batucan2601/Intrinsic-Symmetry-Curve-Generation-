#pragma once
#include "glm/glm.hpp"
#include <map>
#include "Skeleton.h"
#include "Mesh.h"



struct SkeletalNLateral
{
	Skeleton skeleton;
	int N;
	std::vector<unsigned int> point_indices;
	std::vector<std::vector<double>>  geodesic_distances;
	SkeletalNLateral(Mesh& mesh, const std::vector<unsigned int>& point_indices, int N);

};

