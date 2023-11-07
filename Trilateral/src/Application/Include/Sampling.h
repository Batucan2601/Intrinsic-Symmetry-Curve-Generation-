#pragma once 
#include "MeshFactory.h"
#include <vector>

//this function returns a vector of integers
// where size of it is equal to no_of_samples in int
std::vector<unsigned int>  furthest_point_sampling(Mesh* m, int no_of_samples);
std::vector<unsigned int>  furthest_point_sampling_on_partial_points(Mesh* m, int no_of_samples, std::vector<unsigned int>& partial_points);
std::vector<unsigned int>  random_symmetry_indices_sampling(Mesh* m, int no_of_samples);
