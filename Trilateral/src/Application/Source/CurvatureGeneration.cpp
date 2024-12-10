#include "../Include/CurvatureGeneration.h"
#include "../Include/Geodesic.h"
Curvature CurvatureGeneration_generate(TrilateralMesh* m)
{
	// 1 - get the resemblance points mid points
	std::vector<int> curv_points;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int p1 = m->calculated_symmetry_pairs[i].first;
		int p2 = m->calculated_symmetry_pairs[i].second;
		std::vector<int> point_list = Geodesic_between_two_points(*m, p1 , p2);
		curv_points.push_back(point_list[point_list.size() / 2]);
	}
	// 2 - the algorithm is the following
	// do dijkstra for each
	// sum up to other curv points
	// start from the biggest summed point
	float longest_sum = -INFINITY;
	int curv_index = -1;
	std::vector<float> best_distances;
	for (size_t i = 0; i < curv_points.size(); i++)
	{
		int index = curv_points[i];
		std::vector<float> distances = Geodesic_dijkstra(*m,index);
		float sum = 0;
		for (size_t j = 0; j  < curv_points.size(); j ++)
		{
			sum += distances[curv_points[j]];
		}
		if (sum > longest_sum)
		{
			longest_sum = sum ;
			curv_index = i; 
			best_distances = distances;
		}
	}
	std::vector<std::pair<float , unsigned int>> best_distances_selected; 
	for (size_t i = 0; i < curv_points.size(); i++)
	{
		int index = curv_points[i];
		float dist = best_distances[index];
		best_distances_selected.push_back(std::make_pair(dist, index));
	}
	//sort
	std::sort(best_distances_selected.begin(), best_distances_selected.end());

	std::vector<glm::vec3> result;
	for (size_t i = 0; i < best_distances_selected.size()-1; i++)
	{
		int index1 = best_distances_selected[i].second;
		int index2 = best_distances_selected[i + 1].second;
		std::vector<int> points_between =  Geodesic_between_two_points(*m, index1, index2);
		for (size_t i = 0; i < points_between.size(); i++)
		{
			result.push_back(m->vertices[points_between[i]]);
		}
	}


	Curvature curv;
	curv.points = result;
	return curv; 
}