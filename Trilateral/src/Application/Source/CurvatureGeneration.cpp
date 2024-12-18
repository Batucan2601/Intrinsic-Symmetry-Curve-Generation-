#include "../Include/CurvatureGeneration.h"
#include "../Include/Geodesic.h"
#define _USE_MATH_DEFINES
#include "math.h"
static std::vector<unsigned int> get_mid_points(TrilateralMesh* m)
{
	// 1 - get the resemblance points mid points
	std::vector<unsigned int> curv_points;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int p1 = m->calculated_symmetry_pairs[i].first;
		int p2 = m->calculated_symmetry_pairs[i].second;
		std::vector<int> point_list = Geodesic_between_two_points(*m, p1, p2);
		float total_length = 0;
		for (size_t j = 0; j < point_list.size() - 1; j++)
		{
			int index1 = point_list[j];
			int index2 = point_list[j + 1];
			total_length += glm::distance(m->vertices[index1], m->vertices[index2]);
		}
		int halfway_index = -1;
		float dist = 0;
		for (size_t j = 0; j < point_list.size(); j++)
		{
			int index1 = point_list[j];
			int index2 = point_list[j + 1];

			dist += glm::distance(m->vertices[index1], m->vertices[index2]);
			if (dist >= total_length / 2)
			{
				halfway_index = j;
				break;
			}
		}
		curv_points.push_back(point_list[halfway_index]);
	}
	return curv_points;
}
Curvature CurvatureGeneration_generate(TrilateralMesh* m , float merge_distance_param)
{
	Curvature curv;
	// 1 - get the resemblance points mid points
	std::vector<unsigned int> curv_points = get_mid_points(m);
	std::vector<unsigned int>unwanted_pairs; 
	std::vector<unsigned int>wanted_pairs; 
	for (size_t i = 0; i < curv_points.size(); i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, curv_points[i]);
		for (size_t j = 0; j < curv_points.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			if (distances_i[curv_points[j]] < merge_distance_param)
			{
				//check if in wanted pairs
				bool in_wanted_pairs = false;
				for (size_t k = 0; k < wanted_pairs.size(); k++)
				{
					if (wanted_pairs[k] == i)
					{
						in_wanted_pairs = true;
						break; 
					}
				}
				if (!in_wanted_pairs)
				{
					unwanted_pairs.push_back(j);
					wanted_pairs.push_back(i);
				}
			}
		}
	}
	std::sort(unwanted_pairs.begin(), unwanted_pairs.end());
	for (size_t i = 0; i < unwanted_pairs.size(); i++)
	{
		curv_points.erase(curv_points.begin(), curv_points.begin() + unwanted_pairs[unwanted_pairs.size()-i-1]);
	}

	//now get closest geodesic distances
	for (size_t i = 0; i < curv_points.size(); i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, curv_points[i]);
		std::vector<float> distances; 
		for (size_t j = 0; j < curv_points.size(); j++)
		{
			if (i == j)
			{
				distances.push_back(INFINITY);
			}
			distances.push_back(distances_i[curv_points[j]]);
		}
		auto min_element = std::min_element(distances.begin(), distances.end());
		auto min_value = std::min(distances.begin(), distances.end());
		std::vector<int> path = Geodesic_between_two_points(*m, curv_points[i],  curv_points[*min_element]);
		for (size_t i = 0; i < path.size(); i++)
		{
			curv.points.push_back(m->vertices[path[i]]);
		}
	}
	for (size_t i = 0; i < curv_points.size(); i++)
	{
		curv.curvature_main_points.push_back(curv_points[i]);
	}
	return curv; 

}
float CurvatureGeneration_curvature_quality(TrilateralMesh* m, Curvature& curv)
{
	std::vector<unsigned int> curv_points;
	float distortion = 0;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int p1 = m->calculated_symmetry_pairs[i].first;
		int p2 = m->calculated_symmetry_pairs[i].second;
		distortion += std::abs(m->areas[p1] - m->areas[p2]);
	}
	float curv_dist = 0;
	for (size_t i = 0; i < curv.points.size() - 1; i++)
	{
		glm::distance(curv.points[i], curv.points[i+1]);
	}
	curv_dist /= curv.points.size();
	float ct = 0.5;
	float quality = (curv_dist - ct) / (1 - ct );
	curv.curvature_quality = quality;
	return quality; 
}
std::vector<Curve> CurvatureGeneration_generate_curve_paths(TrilateralMesh* m)
{
	std::vector<Curve> curves; 
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		Curve c;
		int p1 = m->calculated_symmetry_pairs[i].first;
		int p2 = m->calculated_symmetry_pairs[i].second;
		std::vector< int> path  = Geodesic_between_two_points(*m, p1, p2);
		for (size_t i = 0; i < path.size(); i++)
		{
			c.curve_path.push_back(path[i]);
		}
		CurvatureGeneration_get_curve_length(m, c);
		curves.push_back(c);
	}

	return curves;
}

float  CurvatureGeneration_get_curve_length(TrilateralMesh* m, Curve& curv)
{
	std::vector<int> path = Geodesic_between_two_points(*m, curv.curve_path[0], curv.curve_path[curv.curve_path.size() - 1]);
	float dist = 0;
	for (size_t i = 0; i < path.size()-1; i++)
	{
		int p1 = path[i];
		int p2 = path[i+1];
		float point_dist = glm::distance(m->vertices[p1], m->vertices[p2]);
		dist += point_dist; 
	}
	curv.length = dist; 
	return dist;
}