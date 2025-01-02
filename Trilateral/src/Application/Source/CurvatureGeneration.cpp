#include "../Include/CurvatureGeneration.h"
#include "../Include/Geodesic.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "../Include/Ray.h"
#include <raymath.h>
#include "../Include/CoreTypeDefs.h"

static void build_curvature(TrilateralMesh* m, Curvature& curvature);
static void remove_points_with_same_index(TrilateralMesh* m, Curvature& curvature);
static float get_quality(TrilateralMesh* m, Curvature& c, int index);

static unsigned int current_midpoint; 
static unsigned int current_midpoint_inverse; 

float Curvature::get_avg_quality(TrilateralMesh* m  )
{
	float total = 0;
	for (size_t i = 1; i < this->curve_points.size()-1; i++)
	{
		total += get_quality(m, *this, i);
	}
	return total / this->curve_points.size();
}
std::vector<unsigned int> Curvature::get_strong_points()
{
	std::vector<unsigned int> strong_list;
	for (size_t i = 0; i < this->curve_points.size(); i++)
	{
		for (size_t j = 0; j < this->curve_points.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			unsigned int midpoint_i = this->curve_points[i].mid_point;
			unsigned int midpoint_j = this->curve_points[j].mid_point;
			if (midpoint_i == midpoint_j)
			{
				//check already added
				bool is_already_added = false;
				for (size_t k = 0; k < strong_list.size(); k++)
				{
					if (strong_list[k] == midpoint_i)
					{
						is_already_added = true; 
						break; 
					}
				}
				if (!is_already_added)
				{
					strong_list.push_back(midpoint_i);
				}
			}
		}
	}
	return strong_list;
}
void Curvature::add_strong_list(std::vector<unsigned int> strong_list)
{
	for (size_t i = 0; i < strong_list.size(); i++)
	{
		CurvePoints p;
		p.mid_point = strong_list[i];
		p.is_strong = true; 
		this->curve_points.push_back(p);
	}
}

void Curvature::generate_curve_quality(TrilateralMesh* m )
{
	this->curve_quality.clear();
	this->curve_quality.push_back(1.0f);

	int first_index = this->curve_points[0].mid_point;
	std::vector<float> distances_from_first = Geodesic_dijkstra(*m, first_index);
	float distance_so_far = 0;
	for (size_t i = 1; i < this->curve_points.size(); i++)
	{
		int index = this->curve_points[i].mid_point;
		float nominator = distances_from_first[index];
		std::vector<float> distances_from_i = Geodesic_dijkstra(*m, index);
		float distance = distances_from_i[this->curve_points[i-1].mid_point];
		distance_so_far += distance;
		float quality = nominator / distance_so_far; 
		this->curve_quality.push_back(quality);
	}

}
static float get_quality(TrilateralMesh* m, Curvature& c , int index  )
{
	unsigned int index_m1 = c.curve_points[index - 1].mid_point;
	unsigned int index_i = c.curve_points[index ].mid_point;
	unsigned int index_p1 = c.curve_points[index + 1].mid_point;
	std::vector<float> distances_p1 = Geodesic_dijkstra(*m, index_p1);
	std::vector<float> distances_index = Geodesic_dijkstra(*m, index_i);

	float dist_on_curve = distances_p1[index_i] +distances_index[index_m1];
	float true_dist = distances_p1[index_m1];

	return  true_dist / dist_on_curve;
	
}

static bool remove_index_w_correspon(TrilateralMesh* m, Curvature& c, unsigned int i)
{
	Curvature curve_temp = c; 
	CurvePoints temp = curve_temp.curve_points[i];
	curve_temp.curve_points.erase(curve_temp.curve_points.begin() + i);
	if (!temp.is_strong)
	{
		for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
		{
			if (temp.correspondence.first == m->calculated_symmetry_pairs[i].first &&
				temp.correspondence.second == m->calculated_symmetry_pairs[i].second)
			{
				m->calculated_symmetry_pairs.erase(m->calculated_symmetry_pairs.begin() + i);
				break;
			}
		}
	}
	build_curvature(m, curve_temp);

	int no_of_normal = 0;
	int no_of_inv_normal = 0;
	for (size_t i = 0; i < curve_temp.paths.size(); i++)
	{
		for (size_t j = 0; j < curve_temp.paths[i].size(); j++)
		{
			int index = curve_temp.paths[i][j];
			glm::vec3 normal = m->normals[index];

			float dot = glm::dot(m->normals[curve_temp.midpoint_index], m->normals[index]);
			float dot_inv = glm::dot(m->normals[curve_temp.midpoint_inv_index], m->normals[index]);
			if (dot > dot_inv)
			{
				no_of_normal++;
			}
			else
			{
				no_of_inv_normal++;
			}
		}
	}

	if (no_of_normal > no_of_inv_normal)
	{
		c = curve_temp;
		return true; 
	}
	return false;
}
static CurvePoints get_worst_curv_point(Curvature& c ,int& index)
{
	// 1 - get the vertex with worst quality. 
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
	index = worst_quality_index;
	return c.curve_points[worst_quality_index];
}
static Curvature connect_front_and_back(TrilateralMesh* m, Curvature& front, Curvature& back)
{
	//detemrine the ends
	CurvePoints curv_front_start = front.curve_points[0];
	CurvePoints curv_front_end = front.curve_points[front.curve_points.size()-1];
	CurvePoints curv_back_start = back.curve_points[0];
	CurvePoints curv_back_end = back.curve_points[back.curve_points.size()-1];

	std::vector<float> distances_from_front_start = Geodesic_dijkstra(*m , curv_front_start.mid_point);
	bool is_front_start_closer_to_back_start = false;
	if (distances_from_front_start[curv_back_end.mid_point] > distances_from_front_start[curv_back_start.mid_point])
	{
		is_front_start_closer_to_back_start = !is_front_start_closer_to_back_start;
	}
	Curvature full_curvature = front; 
	if (is_front_start_closer_to_back_start)
	{
		for (size_t i = 0; i < back.curve_points.size(); i++)
		{
			full_curvature.curve_points.push_back(back.curve_points[back.curve_points.size() - 1 - i]);

		}
	}
	else
	{
		for (size_t i = 0; i < back.curve_points.size(); i++)
		{
			full_curvature.curve_points.push_back(back.curve_points[i]);
		}
	}
	//remove_points_with_same_index(m, full_curvature);
	full_curvature.paths.clear();
	//shall generate paths by hand 
	for (size_t i = 0; i < full_curvature.curve_points.size()-1; i++)
	{
		std::vector<int> paths = Geodesic_between_two_points(*m,
		full_curvature.curve_points[i].mid_point, full_curvature.curve_points[i+1].mid_point);
		std::vector<unsigned int> paths_unsigned;
		for (size_t i = 0; i < paths.size(); i++)
		{
			paths_unsigned.push_back(paths[i]);
		}
		full_curvature.paths.push_back(paths_unsigned);
	}
	//first and last
	std::vector<int> paths = Geodesic_between_two_points(*m,
	full_curvature.curve_points[0].mid_point, full_curvature.curve_points[full_curvature.curve_points.size()-1].mid_point);
	std::vector<unsigned int> paths_unsigned;
	for (size_t i = 0; i < paths.size(); i++)
	{
		paths_unsigned.push_back(paths[i]);
	}
	full_curvature.paths.push_back(paths_unsigned);
	
	full_curvature.generate_curve_quality(m);
	return full_curvature;
}
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
			halfway_index = point_list[j];
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
static float get_quality_sum_avg(Curvature& c)
{
	float sum = 0;
	for (size_t i = 0; i < c.curve_points.size(); i++)
	{
		sum += c.curve_points[i].quality;
	}
	return sum / c.curve_points.size();
}

static CurvePoints get_smallest_quality(Curvature& c)
{
	float smallest = INFINITY;
	int smallest_index = -1;
	for (size_t i = 0; i < c.curve_points.size(); i++)
	{
		float quality = c.curve_points[i].quality;
		if (quality < smallest)
		{
			smallest = quality; 
			smallest_index = i;
		}
	}
	return c.curve_points[smallest_index];
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

	//determine the point's strength
	for (size_t i = 0; i < N; i++)
	{
		curvature.curve_points[i].strength = 0;
	}
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				continue; 
			}
			unsigned int midpoint = curvature.curve_points[i].mid_point;
			unsigned int midpoint_i = curvature.curve_points[j].mid_point;
			if (midpoint == midpoint_i) //same
			{
				curvature.curve_points[i].strength += 1; 
			}
		}
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
bool distance_sort_func(std::pair<float, CurvePoints>& p1 , std::pair<float, CurvePoints>& p2 )
{
	return p2.first > p1.first;
}
static void build_curvature(TrilateralMesh* m , Curvature& curvature)
{
	curvature.paths.clear();
	//start generating curves
	std::vector<float> midpoint_from_sum;
	for (size_t i = 0; i < curvature.curve_points.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, curvature.curve_points[i].mid_point);
		float sum = 0;
		for (size_t j = 0; j < curvature.curve_points.size(); j++)
		{
			sum += distances[curvature.curve_points[j].mid_point];
		}
		midpoint_from_sum.push_back(sum);
	}
	auto best_midpoint_from_sum = std::max_element(midpoint_from_sum.begin(), midpoint_from_sum.end());
	int best_midpoint_from_index = (std::distance(midpoint_from_sum.begin(), best_midpoint_from_sum));
	std::vector<float> distances = Geodesic_dijkstra(*m, curvature.curve_points[best_midpoint_from_index].mid_point);
	std::vector<std::pair<float,CurvePoints>> distances_from_best;
	for (size_t i = 0; i < curvature.curve_points.size(); i++)
	{
		distances_from_best.push_back(std::make_pair(distances[curvature.curve_points[i].mid_point]
		, curvature.curve_points[i]));
	}
	curvature.curve_points.clear();
	std::sort(distances_from_best.begin(), distances_from_best.end(), distance_sort_func);
	for (size_t i = 0; i < distances_from_best.size() - 1; i++)
	{
		int index1 = distances_from_best[i].second.mid_point;
		int index2 = distances_from_best[i + 1].second.mid_point;
		std::vector<int> path = Geodesic_between_two_points(*m, index1, index2);
		std::vector<unsigned int> path_points; 
		for (size_t j = 0; j < path.size(); j++)
		{
			path_points.push_back(path[j]);
		}
		curvature.paths.push_back(path_points);
	}
	for (size_t i = 0; i < distances_from_best.size(); i++)
	{
		curvature.curve_points.push_back(distances_from_best[i].second);

	}

	m->color_points(std::vector<unsigned int>{curvature.curve_points[best_midpoint_from_index].mid_point} , BLUE );
	CurvatureGeneration_curvature_quality(m, curvature);
	curvature.generate_curve_quality(m);
	
}
static void update_curvature(TrilateralMesh* m, Curvature& c, 
unsigned int midpoint_selected, unsigned int midpoint_other, std::vector<unsigned int> agd_vertices , float hks_param , int histogram_size)
{
	// 1 - get the vertex with worst quality. 
	int worst_quality_index;
	CurvePoints least_c = get_worst_curv_point(c , worst_quality_index);
	float worst_quality = least_c.quality;
	
	NLateralDescriptor best_desc;
	NLateralDescriptor best_desc_counterpart;
	std::pair<int, int> result = { -1,-1 };
	unsigned int selected_midpoint = midpoint_other;
	float smallest_hist_dif = INFINITY;
	std::vector<float> distances = Geodesic_dijkstra(*m, selected_midpoint);
	// this pair is wrong,
	CurvePoints worst_curve_point = c.curve_points[worst_quality_index];

	// 2 - erase the low quality pair
	c.curve_points.erase(c.curve_points.begin() + worst_quality_index);
	//erase from mesh
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		std::pair<unsigned int, unsigned int> pair = m->calculated_symmetry_pairs[i];
		if (pair.first == worst_curve_point.correspondence.first && pair.second == worst_curve_point.correspondence.second)
		{
			m->calculated_symmetry_pairs.erase(m->calculated_symmetry_pairs.begin() + i);
			break;
		}
	}

	unsigned int index1 = worst_curve_point.correspondence.first;
	for (size_t i = 0; i < agd_vertices.size(); i++)
	{
		if (index1 == agd_vertices[i] || agd_vertices[i] == worst_curve_point.correspondence.second)
		{
			continue;
		}
		if (std::abs(m->normalized_heat_kernel_signature[index1] - m->normalized_heat_kernel_signature[agd_vertices[i]]) > hks_param)
		{
			continue; 
		}

		//check their distances from midpoint
		float dist_to_midpoint = distances[index1] / distances[agd_vertices[i]];
		if (dist_to_midpoint > 1)
		{
			dist_to_midpoint = 1 / dist_to_midpoint;
		}
		if(dist_to_midpoint < 0.8 )
		{
			continue;
		}
			
		unsigned int mid_index = get_single_mid_point(m, index1, agd_vertices[i]);
		float dot_front = glm::dot(m->normals[midpoint_selected], m->normals[mid_index]);
		float dot_back = glm::dot(m->normals[midpoint_other], m->normals[mid_index]);
		if (dot_back > dot_front)
		{
			continue; 
		}
		//also check the hks
		//generate two trialteral
		std::vector<unsigned int> indices;
		indices.push_back(index1);
		indices.push_back(selected_midpoint);
		indices.push_back(agd_vertices[i]);

		NLateralDescriptor desc = NLateral_generate_descriptor(m, indices);
		desc.create_histogram_area(m, histogram_size,0);
		desc.create_histogram_HKS(m, histogram_size,0);
		desc.area_histogram[0].normalize(1);
		desc.hks_histogram[0].normalize(1);
		//swap 
		indices.clear();
			
		indices.push_back(agd_vertices[i]);
		indices.push_back(selected_midpoint);
		indices.push_back(index1);
			
		NLateralDescriptor desc1 = NLateral_generate_descriptor(m, indices);
		desc1.create_histogram_area(m, histogram_size,0);
		desc1.create_histogram_HKS(m, histogram_size,0);
		desc1.area_histogram[0].normalize(1);
		desc1.hks_histogram[0].normalize(1);

		float dif_area = Histogram_ChiSquareDistance(desc.area_histogram[0], desc1.area_histogram[0]);
		float dif_hks = Histogram_ChiSquareDistance(desc.hks_histogram[0], desc1.hks_histogram[0]);
		float total_dif = (dif_area * dif_area) + (dif_hks * dif_hks);

		if (total_dif < smallest_hist_dif)
		{
			result.first = index1;
			result.second = agd_vertices[i];
			smallest_hist_dif = total_dif;
			best_desc = desc;
			best_desc_counterpart = desc1;
		}
	}
	//if not unmatched 
	if (result.second != -1)
	{
		CurvePoints p;
		p.correspondence = std::make_pair(result.first, result.second);
		unsigned int new_mid_point = get_single_mid_point(m, result.first, result.second);
		p.mid_point = new_mid_point;
		c.curve_points.push_back(p);
	}
	build_curvature(m, c);

	// match second point now
	index1 = worst_curve_point.correspondence.second;
	smallest_hist_dif = INFINITY;
	for (size_t i = 0; i < agd_vertices.size(); i++)
	{
		if (index1 == agd_vertices[i] || agd_vertices[i] == worst_curve_point.correspondence.first)
		{
			continue;
		}
		if (std::abs(m->normalized_heat_kernel_signature[index1] - m->normalized_heat_kernel_signature[agd_vertices[i]]) > hks_param)
		{
			continue;
		}
		//check their distances from midpoint
		float dist_to_midpoint = distances[index1] / distances[agd_vertices[i]];
		if (dist_to_midpoint > 1)
		{
			dist_to_midpoint = 1 / dist_to_midpoint;
		}
		if (dist_to_midpoint < 0.8)
		{
			continue;
		}
		unsigned int mid_index = get_single_mid_point(m, index1, agd_vertices[i]);
		float dot_front = glm::dot(m->normals[midpoint_selected], m->normals[mid_index]);
		float dot_back = glm::dot(m->normals[midpoint_other], m->normals[mid_index]);
		if (dot_back > dot_front)
		{
			continue;
		}
		//also check the hks
		//generate two trialteral
		std::vector<unsigned int> indices;
		indices.push_back(index1);
		indices.push_back(selected_midpoint);
		indices.push_back(agd_vertices[i]);

		NLateralDescriptor desc = NLateral_generate_descriptor(m, indices);
		desc.create_histogram_area(m, histogram_size,0);
		desc.create_histogram_HKS(m, histogram_size,0);
		desc.area_histogram[0].normalize(1);
		desc.hks_histogram[0].normalize(1);
		//swap 
		indices.clear();

		indices.push_back(agd_vertices[i]);
		indices.push_back(selected_midpoint);
		indices.push_back(index1);

		NLateralDescriptor desc1 = NLateral_generate_descriptor(m, indices);
		desc1.create_histogram_area(m, histogram_size, 0);
		desc1.create_histogram_HKS(m, histogram_size,0);
		desc1.area_histogram[0].normalize(1);
		desc1.hks_histogram[0].normalize(1);

		float dif_area = Histogram_ChiSquareDistance(desc.area_histogram[0], desc1.area_histogram[0]);
		float dif_hks = Histogram_ChiSquareDistance(desc.hks_histogram[0], desc1.hks_histogram[0]);
		float total_dif = (dif_area * dif_area) + (dif_hks * dif_hks);

		if (total_dif < smallest_hist_dif)
		{
			result.first = index1;
			result.second = agd_vertices[i];
			smallest_hist_dif = total_dif;
			best_desc = desc;
			best_desc_counterpart = desc1;
		}
	}
	//if not unmatched 
	if (result.second != -1)
	{
		CurvePoints p;
		p.correspondence = std::make_pair(result.first, result.second);
		unsigned int new_mid_point = get_single_mid_point(m, result.first, result.second);
		p.mid_point = new_mid_point;
		c.curve_points.push_back(p);
	}
	build_curvature(m, c);

	//Nlateral_display_desc(m, best_desc);
	//Nlateral_display_desc(m, best_desc_counterpart);
	//m->color_points(std::vector<unsigned int>(result.first, result.second), WHITE);
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
Curvature CurvatureGeneration_generate(TrilateralMesh* m , std::vector<unsigned int>& agd_indices,
float hks_param , float quality_param)
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
	current_midpoint = midpoint_back;
	current_midpoint_inverse = midpoint_front;
	//front
	remove_points_with_same_index(m, curv_front);
	build_curvature(m, curv_front);
	/*while (get_smallest_quality(curv_front).quality < quality_param)
	{
		update_curvature(m, curv_front, midpoint_front, midpoint_back, agd_indices, hks_param, 10);
	}
	CurvatureGeneration_curvature_quality(m, curv_front);*/
	// back 
	remove_points_with_same_index(m, curv_back);
	build_curvature(m, curv_back);
	/*while (get_smallest_quality(curv_back).quality < quality_param)
	{
		update_curvature(m, curv_back, midpoint_back, midpoint_front, agd_indices, hks_param, 10);
	}*/
	CurvatureGeneration_curvature_quality(m, curv_back);

	//Curvature full_curvature = connect_front_and_back(m, curv_front, curv_back);
	m->color_points(std::vector<unsigned int>{ midpoint_front }, WHITE);
	m->color_points(std::vector<unsigned int>{ midpoint_back }, WHITE);

	


	return curv_back;
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
		std::vector<float> dist_1 = Geodesic_dijkstra(*m, index_p1);
		geodesic_distance = dist_1[index_p3]; // dsit 1 - 3 
		std::vector<float> dist_2 = Geodesic_dijkstra(*m, index_p2);
		algorithm_distance += dist_1[index_p2];// dist 1 - 2 
		algorithm_distance += dist_2[index_p3]; // dist 2 -3 
		float quality_of_point = geodesic_distance / algorithm_distance;
		curv.curve_points[i].quality= quality_of_point;
		std::cout << " current point " << index_p2 << std::endl;
		std::cout << " dijksra dist " << geodesic_distance << std::endl;
		std::cout << " alg dist " << dist_1[index_p2]  << " " << dist_2[index_p3] << std::endl;
		std::cout << "alg dist sum " << algorithm_distance << std::endl;
		std::cout << "divided  " << geodesic_distance / algorithm_distance << std::endl;
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

	p1 = smallest_index;
	p2 = smallest_hit_index;
	
}

void CurvatureGeneration_update(TrilateralMesh* m,Curvature& c, std::vector<unsigned int>& agd_indices,
	float hks_param, float quality_param)
{
	update_curvature(m, c, current_midpoint, current_midpoint_inverse,
		agd_indices, hks_param, 10);
}

void CurvatureGeneration_update_w_quality(TrilateralMesh* m, Curvature& c, std::vector<unsigned int>& agd_indices,
	float hks_param, float quality_param)
{
	while (get_smallest_quality(c).quality < quality_param)
	{
		update_curvature(m, c, current_midpoint, current_midpoint_inverse,
			agd_indices, hks_param, 10);
	}

}
static void get_base_points_from_full_curvature(TrilateralMesh* m , Curvature& c , unsigned int& index1,
unsigned int& index2 )
{
	//find the biggest distanced to everyone and it's most distanced counterpart
	std::vector<float> midpoint_from_sum;
	for (size_t i = 0; i < c.curve_points.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, c.curve_points[i].mid_point);
		float sum = 0;
		for (size_t j = 0; j < c.curve_points.size(); j++)
		{
			sum += distances[c.curve_points[j].mid_point];
		}
		midpoint_from_sum.push_back(sum);
	}
	auto best_midpoint_from_sum = std::max_element(midpoint_from_sum.begin(), midpoint_from_sum.end());
	int best_midpoint_from_index = (std::distance(midpoint_from_sum.begin(), best_midpoint_from_sum));
	std::vector<float> distances = Geodesic_dijkstra(*m, c.curve_points[best_midpoint_from_index].mid_point);

	float biggest = -INFINITY;
	int biggest_index = -1;
	for (size_t i = 0; i < c.curve_points.size(); i++)
	{
		int index = c.curve_points[i].mid_point;
		float dist = distances[index];
		if (dist > 1e-12 && dist > biggest)
		{
			biggest = dist;
			biggest_index = i;
		}
	}

	index1 =  best_midpoint_from_index;
	index2 = biggest_index;
}
Curvature CurvatureGeneration_generate_update_full_curv(TrilateralMesh* m, Curvature& c, std::vector<unsigned int>& agd_vertices,
float hks_param , int histogram_size )
{
	unsigned int base_index1;
	unsigned int base_index2;
	get_base_points_from_full_curvature(m, c, base_index1, base_index2);
	//get_full_curve_quality( m ,c , base_index1 , base_index2);
	return c; 
}
Curvature CurvatureGeneration_generate_full_curv(TrilateralMesh* m, std::vector<unsigned int>& agd_indices,
	float hks_param, float quality_param)
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
		float dot_back = glm::dot(normal, m->normals[midpoint_back]);
		if (dot_front >= dot_back)
		{
			CurvePoints p = main_curv.curve_points[i];
			p.strength = 0;
			curv_front.curve_points.push_back(p);
		}
		else
		{
			CurvePoints p = main_curv.curve_points[i];
			p.strength = 0;
			curv_back.curve_points.push_back(p);
		}
	}

	//get strong indices
	std::vector<unsigned int> strong_list = curv_back.get_strong_points();
	std::vector<unsigned int> strong_list_f = curv_front.get_strong_points();
	if (strong_list_f.size() > 0)
	{
		strong_list.insert(strong_list.end(), strong_list_f.begin(), strong_list_f.end());
	}


	//front
	curv_front.midpoint_index = midpoint_front;
	curv_front.midpoint_inv_index = midpoint_back;
	curv_front.add_strong_list(strong_list);

	remove_points_with_same_index(m, curv_front);
	build_curvature(m, curv_front);
	//CurvatureGeneration_update_w_quality(m, curv_front, agd_indices, hks_param, quality_param);
	// back 
	curv_back.midpoint_index = midpoint_back;
	curv_back.midpoint_inv_index = midpoint_front;
	curv_back.add_strong_list(strong_list);

	remove_points_with_same_index(m, curv_back);
	build_curvature(m, curv_back);
	//CurvatureGeneration_update_w_quality(m, curv_back, agd_indices, hks_param, quality_param);

	while(CurvatureGeneration_curve_smoothing(m, curv_front, quality_param)); 
	while(CurvatureGeneration_curve_smoothing(m, curv_back, quality_param)); 

	Curvature full_curvature = connect_front_and_back(m, curv_front, curv_back);
	m->color_points(std::vector<unsigned int>{ midpoint_front }, WHITE);
	m->color_points(std::vector<unsigned int>{ midpoint_back }, WHITE);

	m->color_points(std::vector<unsigned int>{ strong_list }, RED);

	return full_curvature;
}

void CurvatureGeneration_laplacian_smoothing(TrilateralMesh* m, Curvature& c, float quality_param )
{
	float worst_quality = INFINITY; 
	int worst_quality_index = -1; 
	for (size_t i = 1; i < c.curve_points.size()-1; i++)
	{
		int index_i_m = c.curve_points[i-1].mid_point; 
		int index_i_p = c.curve_points[i+1].mid_point;
		int index = c.curve_points[i].mid_point;
		
		// find midpoint 
		unsigned int new_index = Geodesic_find_midpoint(m, index_i_m, index_i_p);
		c.curve_points[i].mid_point = index; 
		float quality_i = get_quality(m, c, i);
		if (quality_i  < quality_param)
		{
			//remove point and correspondence
			if (quality_i < worst_quality  )
			{
				worst_quality = quality_i;
				worst_quality_index = i; 
			}
		}
	} 
	if (worst_quality_index != -1 )
	{
		remove_index_w_correspon(m, c, worst_quality_index);
	}

}
static bool check_parameters(TrilateralMesh* m, Curvature& c, unsigned int index1, unsigned int index2,
	float hks_param, float distance_to_midpoint_param  )
{
	float hks_dif = std::abs(m->normalized_heat_kernel_signature[index1] - m->normalized_heat_kernel_signature[index2]);
	bool is_hks = hks_dif < hks_param; 

	std::vector<float> distances = Geodesic_dijkstra(*m,c.midpoint_index);
	float dist1 = distances[index1];
	float dist2 = distances[index2];
	float ratio = dist1 / dist2;
	if (ratio > 1)
	{
		ratio = 1 / ratio;
	}
	if (ratio > distance_to_midpoint_param && is_hks )
	{
		return true;
	}
	return false;

}
void CurvatureGeneration_add_new_matching(TrilateralMesh* m, Curvature& c, 
std::vector<unsigned int>& agd_indices , float quality_param , float hks_param , float distance_to_midpoint_param)
{
	std::vector<unsigned int> currently_avaliable_indices;
	//select currently avaliable indices
	for (size_t i = 0; i < agd_indices.size(); i++)
	{
		unsigned int curent_index = agd_indices[i]; 
		bool is_index_exist = false; 
		for (size_t j = 0; j < m->calculated_symmetry_pairs.size(); j++)
		{
			unsigned int index1 = m->calculated_symmetry_pairs[j].first;
			unsigned int index2 = m->calculated_symmetry_pairs[j].second;
			if (index1 == curent_index || index2 == curent_index)
			{
				is_index_exist = true; 
				break; 
			}
		}
		if (!is_index_exist)
		{
			currently_avaliable_indices.push_back(agd_indices[i]);
		}
	}
	glm::vec3 normal_midpoint = m->normals[c.midpoint_index];
	glm::vec3 normal_inv_midpoint = m->normals[c.midpoint_inv_index];

	float best_quality =0 ;
	std::pair<unsigned int,unsigned int> best_quality_indices = std::make_pair(0,0);
	for (size_t i = 0; i < currently_avaliable_indices.size(); i++)
	{
		int index_i = currently_avaliable_indices[i];
		for (size_t j = 0; j < currently_avaliable_indices.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			int index_j = currently_avaliable_indices[j];
			// 1- get midpoint
			unsigned int midpoint_i_j = get_single_mid_point(m, index_i, index_j);

			//2 - compare normal
			glm::vec3 normal_midpoint_i_j = m->normals[midpoint_i_j];
			float dot_normal = glm::dot(normal_midpoint, normal_midpoint_i_j);
			float dot_normal_inv = glm::dot(normal_inv_midpoint, normal_midpoint_i_j);
			if (dot_normal_inv > dot_normal)
			{
				continue; // failed may belong to other curve 
			}
			Curvature temp = c; 

			bool is_ok = check_parameters(m, c, index_i, index_j, hks_param, distance_to_midpoint_param);

			if (!is_ok)
			{
				continue; 
			}
			//add the new point
			CurvePoints p; 
			p.correspondence = std::make_pair(index_i, index_j);
			p.mid_point = midpoint_i_j;
			
			temp.curve_points.push_back(p);

			build_curvature(m, temp);

			float quality = temp.get_avg_quality(m);

			if (quality > best_quality)
			{
				best_quality = quality;
				best_quality_indices = std::make_pair(index_i, index_j);
			}
		}
	}

	//add the new pair
	unsigned int index_i = best_quality_indices.first;
	unsigned int index_j = best_quality_indices.second;
	if (index_i == 0 && index_j == 0)
	{
		std::cout << " no build " << std::endl;
		return;
	}
	unsigned int midpoint_i_j = get_single_mid_point(m, index_i, index_j);

	CurvePoints p;
	p.mid_point = midpoint_i_j;
	p.correspondence = std::make_pair(index_i, index_j);

	c.curve_points.push_back(p);
	m->calculated_symmetry_pairs.push_back(p.correspondence);

	build_curvature(m,c);
}

bool CurvatureGeneration_curve_smoothing(TrilateralMesh* m, Curvature& c, float quality_dif_param)
{
	if (c.curve_points.size() <= 2)
	{
		return false;
	}
	int index_i = -1;
	for (size_t i = 0; i < c.curve_points.size()-1; i++)
	{
		if (c.curve_quality[i] - c.curve_quality[i + 1] > quality_dif_param)
		{
			index_i = i;
			break;
		}
	}
	if (index_i != -1)
	{
		bool is_ok = remove_index_w_correspon(m, c, index_i);
		if (is_ok)
		{
			return true;
		}
		return false;
	}
	return false;
}