#include "../Include/SkeletalNLateral.h"



//only geodesic distances for skeletal end points 
SkeletalNLateral::SkeletalNLateral(Mesh& mesh, const std::vector<unsigned int>& point_indices, int N)
{
	// 1 - copy data 
	this->point_indices = point_indices;
	this->N = N;

	// 2- calculate distances
}