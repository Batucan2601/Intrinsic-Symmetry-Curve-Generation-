#include "../Include/CurvatureGeneration.h"
#include "../Include/Geodesic.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "../Include/Ray.h"
#include <raymath.h>
#include "../Include/CoreTypeDefs.h"
//eliminates the points with same index
unsigned int get_single_mid_point(TrilateralMesh* m, unsigned int p1, unsigned int p2)
{
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
	return halfway_index;
}
static float get_quality_sum(Curvature& c)
{
	float sum = 0;
	for (size_t i = 0; i < c.curve_points.size(); i++)
	{
		sum += c.curve_points[i].quality;
	}
	return sum;
}
static void remove_points_with_same_index(TrilateralMesh* m , Curvature& curvature)
{
	int N = curvature.curve_points.size();
	std::vector<std::pair<unsigned int, unsigned int >>  midpoint_with_indices;
	for (size_t i = 0; i < N; i++)
	{
		midpoint_with_indices.push_back(std::make_pair(curvature.curve_points[i].mid_point , i));
	}
	std::sort(midpoint_with_indices.begin(), midpoint_with_indices.end());
	//change every curve item
	Curvature temp = curvature;
	for (size_t i = 0; i < N; i++)
	{
		int index = midpoint_with_indices[i].second;
		curvature.curve_points[i] = temp.curve_points[index];
	}

	//remove the elements now
	std::vector<unsigned int> vertices_to_be_removed;
	for (size_t i = 0; i < N-1; i++)
	{
		unsigned int midpoint = curvature.curve_points[i].mid_point;
		unsigned int midpoint_i = curvature.curve_points[i+1].mid_point;
		if (midpoint == midpoint_i) //same
		{
			vertices_to_be_removed.push_back(i);
		}
	}
	std::sort(vertices_to_be_removed.begin(), vertices_to_be_removed.end(), std::greater<unsigned int>());
	for (size_t i = 0; i < vertices_to_be_removed.size(); i++)
	{
		curvature.curve_points.erase(curvature.curve_points.begin() + vertices_to_be_removed[i]);
	}


	
}
static void build_curvature(TrilateralMesh* m , Curvature& curvature)
{
	curvature.points.clear();
	//start generating curves
	std::vector<float> midpoint_fron_sum;
	for (size_t i = 0; i < curvature.curve_points.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, curvature.curve_points[i].mid_point);
		float sum = 0;
		for (size_t j = 0; j < curvature.curve_points.size(); j++)
		{
			sum += distances[curvature.curve_points[j].mid_point];
		}
		midpoint_fron_sum.push_back(sum);
	}
	auto best_midpoint_from_sum = std::max_element(midpoint_fron_sum.begin(), midpoint_fron_sum.end());
	int best_midpoint_from_index = (std::distance(midpoint_fron_sum.begin(), best_midpoint_from_sum));
	std::vector<float> distances = Geodesic_dijkstra(*m, curvature.curve_points[best_midpoint_from_index].mid_point);
	std::vector<std::pair<float, unsigned int>> distances_from_best;
	for (size_t i = 0; i < curvature.curve_points.size(); i++)
	{
		distances_from_best.push_back(std::make_pair(distances[curvature.curve_points[i].mid_point]
		, curvature.curve_points[i].mid_point));
	}
	std::sort(distances_from_best.begin(), distances_from_best.end());
	for (size_t i = 0; i < distances_from_best.size() - 1; i++)
	{
		int index1 = distances_from_best[i].second;
		int index2 = distances_from_best[i + 1].second;
		std::vector<int> path = Geodesic_between_two_points(*m, index1, index2);
		for (size_t j = 0; j < path.size(); j++)
		{
			curvature.points.push_back(m->vertices[path[j]]);
		}
	}
	CurvatureGeneration_curvature_quality(m, curvature);

	
}
static void update_curvature(TrilateralMesh* m, Curvature& c, 
unsigned int midpoint1, unsigned int midpoint_other, std::vector<unsigned int> agd_vertices , float hks_param , int histogram_size)
{
	float worst_quality = INFINITY;
	int  worst_quality_index = -1;
	for (size_t i = 0; i < c.curve_points.size(); i++)
	{
		float quality = c.curve_points[i].quality;
		if (quality < worst_quality)
		{
			worst_quality = quality;
			worst_quality_index = i;
		}
	}
	// this pair is wrong,
	int index1 = c.curve_points[worst_quality_index].correspondence.first;
	//int index2 = res_pair.second; 
	float current_quality_sum = get_quality_sum(c);
	float new_quality_sum = 0;
	int smallest_point_index = -1;
	int midpoint_index = -1;
	NLateralDescriptor best_desc;
	NLateralDescriptor best_desc_counterpart;
	while (current_quality_sum > new_quality_sum)
	{
		unsigned int selected_midpoint = midpoint_other;
		float smallest_hist_dif = INFINITY;
		int index1_pair = -1;
		glm::vec3 index1_normal = m->normals[index1];
		for (size_t i = 0; i < agd_vertices.size(); i++)
		{
			if (index1 == agd_vertices[i])
			{
				continue;
			}
			if (std::abs(m->normalized_heat_kernel_signature[index1] - m->normalized_heat_kernel_signature[agd_vertices[i]]) > hks_param)
			{
				continue; 
			}
			index1_pair = agd_vertices[i];
			//also check the hks
			//generate two trialteral
			std::vector<unsigned int> indices;
			indices.push_back(index1);
			indices.push_back(selected_midpoint);
			indices.push_back(index1_pair);

			NLateralDescriptor desc = NLateral_generate_descriptor(m, indices);
			desc.create_histogram_area(m, histogram_size);
			desc.create_histogram_HKS(m, histogram_size);
			desc.area_histogram.normalize(1);
			desc.hks_histogram.normalize(1);
			//swap 
			indices.clear();
			indices.push_back(index1_pair);
			indices.push_back(selected_midpoint);
			indices.push_back(index1);
			NLateralDescriptor desc1 = NLateral_generate_descriptor(m, indices);
			desc1.create_histogram_area(m, histogram_size);
			desc1.create_histogram_HKS(m, histogram_size);

			desc1.area_histogram.normalize(1);
			desc1.hks_histogram.normalize(1);
			float dif_area = Histogram_ChiSquareDistance(desc.area_histogram, desc1.area_histogram);
			float dif_hks = Histogram_ChiSquareDistance(desc.hks_histogram, desc1.hks_histogram);
			float total_dif = (dif_area * dif_area) + (dif_hks * dif_hks);

			if (total_dif < smallest_hist_dif)
			{
				smallest_hist_dif = total_dif;
				smallest_point_index = index1_pair;
				best_desc = desc;
				best_desc_counterpart = desc1;
			}
		}
		//add the new pair and build the system again
		c.curve_points.erase(c.curve_points.begin() + worst_quality_index);
		if (index1_pair == -1)
		{
			continue; 
		}
		CurvePoints p; 
		p.correspondence = std::make_pair(index1 ,index1_pair );
		unsigned int new_mid_point = get_single_mid_point(m, index1, index1_pair);
		p.mid_point = new_mid_point;

		c.curve_points.push_back(p);
		build_curvature(m, c);
		new_quality_sum =  get_quality_sum(c);
		if (new_quality_sum > current_quality_sum)
		{
			continue;
		}
		//else
		//remove this curvepoint it wont match with anything
		int erased_index = -1;
		for (size_t i = 0; i < c.curve_points.size(); i++)
		{
			if (c.curve_points[i].correspondence.first == index1 &&
				c.curve_points[i].correspondence.second == index1_pair &&
				c.curve_points[i].mid_point == new_mid_point)
			{
				erased_index = i;
				break;
			}
		}
		c.curve_points.erase(c.curve_points.begin() + erased_index);
		build_curvature(m, c);
		// this pair is wrong,
		float worst_quality = INFINITY;
		int  worst_quality_index = -1;
		for (size_t i = 0; i < c.curve_points.size(); i++)
		{
			float quality = c.curve_points[i].quality;
			if (quality < worst_quality)
			{
				worst_quality = quality;
				worst_quality_index = i;
			}
		}
		index1 = worst_quality_index;
	}
	Nlateral_display_desc(m, best_desc);
	Nlateral_display_desc(m, best_desc_counterpart);
	m->color_points(std::vector<unsigned int>(index1, smallest_point_index), WHITE);
}
static Curvature get_mid_points(TrilateralMesh* m)
{
	Curvature c;
	// 1 - get the resemblance points mid points
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
		CurvePoints curveP; 
		curveP.correspondence  = (std::make_pair(p1, p2));
		curveP.mid_point = point_list[halfway_index];
		c.curve_points.push_back(curveP);
	}
	return c;
}
Curvature CurvatureGeneration_generate(TrilateralMesh* m , float merge_distance_param, std::vector<unsigned int>& agd_indices,
float hks_param)
{
	Curvature curv_front;
	Curvature curv_back;
	Curvature main_curv; 
	unsigned int midpoint_front, midpoint_back;
	// 1 - get the resemblance points mid points
	main_curv = get_mid_points(m);
	CurvatureGeneration_mid_point_w_AGD(m, midpoint_front, midpoint_back);
	for (size_t i = 0; i < main_curv.curve_points.size(); i++)
	{
		unsigned int index = main_curv.curve_points[i].mid_point;
		glm::vec3 normal = m->normals[index];

		float dot_front = glm::dot(normal, m->normals[midpoint_front]);
		float dot_back = glm::dot(normal,  m->normals[midpoint_back]);
		if (dot_front >= dot_back)
		{
			CurvePoints p = main_curv.curve_points[i];
			curv_front.curve_points.push_back(p);
		}
		else
		{
			CurvePoints p = main_curv.curve_points[i];
			curv_back.curve_points.push_back(p);
		}
	}
	remove_points_with_same_index(m, curv_front);
	build_curvature(m, curv_front);
	//remove_points_with_same_index(m, curv_back);
	//build_curvature(m, curv_back);
	//for front
	update_curvature(m, curv_front, midpoint_front, midpoint_back,agd_indices, hks_param, 10);

	return curv_front;
}
// the quality is the measure of 
// the following
// the actual path over  
//  the path p-1 and p+1
void CurvatureGeneration_curvature_quality(TrilateralMesh* m, Curvature& curv)
{
	for (size_t i = 0; i < curv.curve_points.size(); i++)
	{
		curv.curve_points[i].quality = 0;
	}
	curv.curve_points[0].quality = 1.0f;
	for (size_t i = 1; i < curv.curve_points.size()-1; i++)
	{
		float geodesic_distance = 0;
		float algorithm_distance = 0;
		int index_p1 = curv.curve_points[i-1].mid_point;
		int index_p2 = curv.curve_points[i].mid_point;
		int index_p3 = curv.curve_points[i+1].mid_point;
		std::vector<int> paths = Geodesic_between_two_points(*m, index_p1 , index_p3);
		for (size_t i = 0; i < paths.size()-1; i++)
		{
			int index_i = paths[i];
			int index_i1 = paths[i+1];
			float dist = glm::distance(m->vertices[index_i] , m->vertices[index_i1]);
			geodesic_distance += dist; 
		}
		std::vector<int> paths_1_2 = Geodesic_between_two_points(*m, index_p1, index_p2);
		std::vector<int> paths_2_3 = Geodesic_between_two_points(*m, index_p2, index_p3);
		for (size_t i = 0; i < paths_1_2.size() - 1; i++)
		{
			int index_i = paths_1_2[i];
			int index_i1 = paths_1_2[i + 1];
			float dist = glm::distance(m->vertices[index_i], m->vertices[index_i1]);
			algorithm_distance += dist;
		}
		for (size_t i = 0; i < paths_2_3.size() - 1; i++)
		{
			int index_i = paths_2_3[i];
			int index_i1 = paths_2_3[i + 1];
			float dist = glm::distance(m->vertices[index_i], m->vertices[index_i1]);
			algorithm_distance += dist;
		}
		float quality_of_point = geodesic_distance / algorithm_distance;
		curv.curve_points[i].quality= quality_of_point;
	}
	curv.curve_points[curv.curve_points.size()-1].quality = 1.0f;
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

//least agd and its counterpart
void CurvatureGeneration_mid_point_w_AGD(TrilateralMesh* m , unsigned int& p1 , unsigned int& p2 )
{
	//find agd for all
	int N = m->vertices.size();
	float smallest_val = INFINITY;
	unsigned int smallest_index = -1;
	for (size_t i = 0; i < N; i++)
	{
		std::vector<float> distances =Geodesic_dijkstra(*m, i);
		float sum = 0;
		for (size_t j = 0; j < N ; j++)
		{
			sum += distances[j];
		}
		if (sum < smallest_val)
		{
			smallest_val = sum;
			smallest_index = i;
		}
	}
	// that is the smallest index 
	//find the counterpart by using a ray
	TrilateralRay ray;
	ray.origin = m->vertices[smallest_index];
	glm::vec3 vec1(m->normals_display[smallest_index * 12], m->normals_display[smallest_index*12 + 1], m->normals_display[smallest_index*12 + 2]);
	glm::vec3 vec2(m->normals_display[smallest_index*12 + 6], m->normals_display[smallest_index*12 + 7], m->normals_display[smallest_index*12 + 8]);
	glm::vec3 dir = vec2 - vec1;
	dir = dir * -1.0f;
	ray.direction = dir;
	float smallest_dist = INFINITY;
	unsigned int smallest_hit_index = -1; 
	for (size_t i = 0; i < m->triangles.size(); i+=3 )
	{
		glm::vec3 hit_point;
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		bool is_hit = ray_triangle_intersection(ray , m->vertices[index1] , m->vertices[index2] , m->vertices[index3], hit_point );
		if (is_hit)
		{
			float dist = glm::distance(m->vertices[smallest_index], hit_point);
			if (smallest_dist > dist)
			{
				dist = smallest_dist;
				//get the closest dist on triangle 
				float disti = glm::distance(hit_point , m->vertices[index1]);
				float disti1 = glm::distance(hit_point, m->vertices[index2]);
				float disti2 = glm::distance(hit_point, m->vertices[index3]);
				if (disti < disti1 && disti < disti2)
				{
					smallest_hit_index = index1;
				}
				else if (disti1 < disti && disti1 < disti2)
				{
					smallest_hit_index = index2; 
				}
				else
				{
					smallest_hit_index = index3; 
				}
			}
		}

	}

	m->color_all(BLACK);
	m->color_points(std::vector<unsigned int>{smallest_index, smallest_hit_index}, RED);
	
	p1 = smallest_index;
	p2 = smallest_hit_index;
	
}