#include "../Include/Voronoi.h"
#include "../Include/Geodesic.h"
#include <stack>
#include <random>
#include "../Include/VarianceMinimizingTransportPlan.h"
#include "../Include/ROI.h"

Voronoi::Voronoi()
{
	return;
}
Voronoi::Voronoi(TrilateralMesh* m, unsigned int p1, unsigned int p2 , float param )
{
	this->m = m;
	this->p1 = p1;
	this->p2 = p2;
	//we have to build indices 
	// part 1 - voronoi boundary
	std::vector<float> distances_1 = Geodesic_dijkstra(*m, p1);
	std::vector<float> distances_2 = Geodesic_dijkstra(*m, p2);

	// lets normalzie boundary
	auto dist_1_max_auto = std::max_element(distances_1.begin(), distances_1.end());
	auto dist_2_max_auto = std::max_element(distances_2.begin(), distances_2.end());
	float dist_1_max = *dist_1_max_auto;
	float dist_2_max = *dist_2_max_auto;
	for (size_t i = 0; i < distances_1.size(); i++)
	{
		distances_1[i] = distances_1[i] / dist_1_max;
		distances_2[i] = distances_2[i] / dist_2_max;
	}

	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		float dif = std::abs(distances_1[i] - distances_2[i]);
		if (dif < param)
		{
			indices.push_back(i);
		}
	}
}
void Voronoi::generate_voronoi_parts()
{
	int no_of_points_visited = 0;

	this->status = std::vector<int>(this->m->vertices.size() , 0 );
	for (size_t i = 0; i < this->indices.size(); i++)
	{
		this->status[this->indices[i]] = -1;
		no_of_points_visited++;
	}
	int current_partition_no = 1; 
	while (no_of_points_visited < this->m->vertices.size())
	{
		//pick a non curve point
		unsigned int random_pick = 0;
		while (this->status[random_pick] != 0) //get a visited index
		{
			random_pick++;
			if (random_pick == this->status.size())
			{
				return;
			}
		}
		std::vector<unsigned int> visited_vertices;
		std::stack<int> stack;
		stack.push(random_pick);
		while (!stack.empty())
		{
			int index = stack.top();
			stack.pop(); //vertex index popped from stack
			if (this->status[index] == 0) //not visited
			{
				this->status[index] = current_partition_no; // now te vertex has been visited
				visited_vertices.push_back(index);
				// this region of loop assumes index is not edge, therefore add the adjacencies
				for (size_t i = 0; i < m->adjacenies[index].size(); i++) //process pairs 
				{
					stack.push(m->adjacenies[index][i].first);
				}
			}
			if (this->status[index] == -1) //do nothing 
			{
				;
			}
		}
		current_partition_no++;
	}
	
	this->no_of_partition = current_partition_no;
}

void Voronoi::area_compare()
{
	float area_red = 0;
	float area_green = 0;
	for (size_t i = 0; i < this->status.size(); i++)
	{
		if (this->status[i] == 1)
		{
			area_red = area_red + this->m->areas[i];
		}
		else if (this->status[i] == 2)
		{
			area_green = area_green + this->m->areas[i];
		}

	}
	std::cout << " area comp " << std::min(area_green, area_red) / std::max(area_green, area_red) << std::endl;
}
void Voronoi::color()
{
	std::vector<unsigned int> voronoi_part; 
	std::vector<unsigned int> half1;
	std::vector<unsigned int> half2;
	std::vector<unsigned int> half3;
	std::vector<unsigned int> half4;
	std::vector<unsigned int> half5;
	for (  int i = 0; i < this->status.size();  i++)
	{
		/*if (this->status[i] == -1)
		{
			voronoi_part.push_back(i);
		}*/
		if (this->status[i] == 1)
		{
			half1.push_back(i);
		}
		else if (this->status[i] == 2)
		{
			half2.push_back(i);
		}
		else if (this->status[i] == 3)
		{
			half3.push_back(i);
		}
		else if (this->status[i] == 4)
		{
			half4.push_back(i);
		}
		else if (this->status[i] == 5)
		{
			half5.push_back(i);
		}
	}
	this->m->color_points(this->indices, YELLOW);
	this->m->color_points(half1, GREEN);
	this->m->color_points(half2, RED);
	this->m->color_points(half3, BLUE);
	this->m->color_points(half4, BLACK);
	this->m->color_points(half5, PINK);
}
float Voronoi::distance_to_index(unsigned int index)
{
	std::vector<float> distances = Geodesic_dijkstra(*this->m, index);
	
	float min = INFINITY;
	for (size_t i = 0; i < this->indices.size(); i++)
	{
		int index = this->indices[i];
		if (distances[index] < min )
		{
			min = distances[index];
		}
	}
	return min; 
}
float NLateral_generate_descriptors_with_random_voronoi_points(
	TrilateralMesh* m, unsigned int p1, unsigned int p2 , float voronoi_param,float fuziness, int hist_no
, int no_of_point_param)
{
	// 1 - create Voronoi
	Voronoi v(m, p1, p2, voronoi_param);
	// Create a non-deterministic random device for seeding
	std::random_device rd;
	// Use the Mersenne Twister engine seeded with rd()
	std::mt19937 gen(rd());
	// Define a uniform integer distribution from 0 to 99
	std::uniform_int_distribution<> dist(0, v.indices.size()-1);
	
	
	if (no_of_point_param > v.indices.size())
	{
		no_of_point_param = v.indices.size();
	}

	std::vector<int > indices(no_of_point_param , -1);
	unsigned int mid_index = Geodesic_get_midpoint_from_path(m, p1, p2);
	indices[0] =  mid_index;
	float max_dist = 0;
	float min_dist = INFINITY;
	int max_index = 0;
	std::vector<float> dist_from_mid = Geodesic_dijkstra(*m, mid_index);
	for (size_t i = 0; i < v.indices.size(); i++)
	{
		int index = v.indices[i];
		if (max_dist < dist_from_mid[index])
		{
			max_dist = dist_from_mid[index];
			max_index = index; 
		}
		if (min_dist > dist_from_mid[index])
		{
			min_dist = dist_from_mid[index];
		}
	}
	indices[no_of_point_param - 1] = max_index;
	float step = ( max_dist - min_dist ) / no_of_point_param;
	for (size_t i = 0; i < v.indices.size(); i++)
	{	
		int index = v.indices[i];
		float dist = dist_from_mid[index];
		
		int  which_index = (dist - min_dist) / step; 
		if (which_index < no_of_point_param)
		{
			if (indices[which_index] == -1 )
			{
				indices[which_index] = index;
			}
		}
		
	}
	for (size_t i = 0; i < indices.size(); i++)
	{
		if (indices[i] == -1)
		{
			indices[i] = dist(gen);
		}
	}
	float total = 0;
	for (size_t i = 0; i < no_of_point_param; i++)
	{
		unsigned int index = indices[i];
		NLateralDescriptor desc1 = NLateral_generate_symmetric_descriptor(m, p1, index, hist_no, fuziness);
		NLateralDescriptor desc2 = NLateral_generate_symmetric_descriptor(m, p2, index, hist_no, fuziness);
		float dif = VarianceMin_compare(m, desc1, desc2, true, hist_no, 1);
		total = total + dif; 
	}

	return total;
}
Voronoi Voronoi_algorithm_in_action(TrilateralMesh* m, float voronoi_param, float hks_param, float sdf_param,float dist_param, int hist_no, float fuziness, int no_of_points,
	std::vector<unsigned int>& agd_indices)
{
	Voronoi v = Voronoi_get_closest_voronoi(m, voronoi_param);
	//Voronoi_prune_voronoi(m, v, voronoi_param);
	v = Voronoi_destroy_wrong_matches_and_recalculate(m, voronoi_param, v);
	v = Voronoi_check_pair_closeness_and_recalculate(m, voronoi_param, dist_param, v);

	v = Voronoi_add_points_and_recalculate(m, v, voronoi_param, hks_param, sdf_param, dist_param, hist_no, fuziness, no_of_points, agd_indices);
	v = Voronoi_destroy_wrong_matches_and_recalculate(m, voronoi_param, v);
	//v = Voronoi_get_closest_voronoi(m, voronoi_param);
	return v; 
}
Voronoi Voronoi_algorithm_in_action(TrilateralMesh* m,Voronoi& v, float voronoi_param, float hks_param, float sdf_param, float dist_param, int hist_no, float fuziness, int no_of_points,
	std::vector<unsigned int>& agd_indices)
{

	v.generate_voronoi_parts();
	v = Voronoi_destroy_wrong_matches_and_recalculate(m, voronoi_param, v);
	v = Voronoi_check_pair_closeness_and_recalculate(m, voronoi_param, dist_param, v);

	v = Voronoi_add_points_and_recalculate(m, v, voronoi_param, hks_param, sdf_param, dist_param, hist_no, fuziness, no_of_points, agd_indices);
	v = Voronoi_destroy_wrong_matches_and_recalculate(m, voronoi_param, v);
	return v;
}


Voronoi Voronoi_add_points_and_recalculate(TrilateralMesh* m, Voronoi& v,float voronoi_param, float hks_param, float sdf_param,float dist_param, int hist_no, float fuziness, int no_of_points
,std::vector<unsigned int>& agd_points)
{
	std::vector<bool> is_in_matching(agd_points.size(), false);
	int selected_agd_index = -1;
	for (size_t i = 0; i < no_of_points; i++)
	{
		for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
		{
			int index1 = m->calculated_symmetry_pairs[i].first;
			int index2 = m->calculated_symmetry_pairs[i].second;
			for (size_t j = 0; j < agd_points.size(); j++)
			{
				if (agd_points[j] == index1)
				{
					is_in_matching[j] = true;
				}
				if (agd_points[j] == index2)
				{
					is_in_matching[j] = true;
				}
			}

		}
		//select a never selected index
		for (size_t i = 0; i < agd_points.size(); i++)
		{
			if (!is_in_matching[i])
			{
				selected_agd_index = agd_points[i];
				is_in_matching[i] = true; 
				break;
			}
		}
		if (selected_agd_index == -1)
		{
			break; 
		}
		Voronoi_add_point(m, v, selected_agd_index, voronoi_param,  hks_param, sdf_param,dist_param, m->midpoint, hist_no, fuziness , agd_points);
	}
	//v = Voronoi_get_closest_voronoi(m, voronoi_param);
	return v; 
}
void Voronoi_add_point(TrilateralMesh* m, Voronoi& voronoi , int selected_agd_index,  float voronoi_param , float hks_param , float sdf_param , float dist_param, unsigned int midpoint , int hist_no, float fuziness
, std::vector<unsigned int>& agd_points)
{
	std::vector<bool> is_in_matching(agd_points.size(), false);
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		for (size_t j = 0; j < agd_points.size(); j++)
		{
			if (agd_points[j] == index1)
			{
				is_in_matching[j] = true;
			}
			if (agd_points[j] == index2)
			{
				is_in_matching[j] = true;
			}
		}

	}

	float hks = m->normalized_heat_kernel_signature[selected_agd_index];
	std::vector<float> distances = Geodesic_dijkstra(*m, m->midpoint);
	//select a counterpart
	std::vector<float> ratios; 
	std::vector<std::pair<float, std::pair<unsigned int,unsigned int>> > results; 
	//for (size_t i = 0; i < agd_points.size(); i++)
	{
		for (size_t j = 0; j < agd_points.size(); j++)
		{
			unsigned int index_j = agd_points[j];
			if (is_in_matching[j] == true || voronoi.status[index_j] == -1 )
			{
				continue;
			}
			if (voronoi.status[index_j] == voronoi.status[selected_agd_index])
			{
				continue; 
			}
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[selected_agd_index] - m->normalized_heat_kernel_signature[index_j]);
			bool is_opposite = voronoi.status[selected_agd_index] == voronoi.status[index_j];
			if (hks_dif < hks_param  && is_opposite)
			{
				continue;
			}
			std::vector<float> distances_index_j = Geodesic_dijkstra(*m, index_j);
			std::vector<float> distances_agd_index = Geodesic_dijkstra(*m , selected_agd_index);
			for (size_t k = 0; k < voronoi.indices.size(); k++) // get 10 distances from indcies
			{
				float dist1 = distances_index_j[voronoi.indices[k]];
				float dist2 = distances_agd_index[voronoi.indices[k]];
				float ratio = std::min(dist1, dist2) / std::max(dist1, dist2);
				ratios.push_back(ratio);
			}
			/*float dif = NLateral_generate_descriptors_with_random_voronoi_points(m, random_index,
				agd_points[j], voronoi_param, fuziness, hist_no, 10);*/
			float dif = 0;
			for (size_t i = 0; i < ratios.size(); i++)
			{
				dif += ratios[i];
			}
			dif = dif / ratios.size();
			if (dif > dist_param)
			{
				results.push_back(std::make_pair(dif, std::make_pair(selected_agd_index, index_j)));
			}
		}
	
		
	}
	if (results.size() == 0)
	{
		return; 
	}
	auto min_index_auto = std::min_element( results.begin(), results.end());
	std::pair<unsigned int, unsigned int> min_index = (*min_index_auto).second;

	std::cout << " newly added " << min_index.first << " " << min_index.second << std::endl;
	m->calculated_symmetry_pairs.push_back(std::make_pair(min_index.first, min_index.second));

}
Voronoi Voronoi_get_closest_voronoi(TrilateralMesh* m, float voronoi_param)
{

	std::vector<std::pair<float,unsigned int>> midpoints;
	std::vector<float> total_distances(m->calculated_symmetry_pairs.size(),0);
	std::vector<float> weights(m->calculated_symmetry_pairs.size(), 0);
	std::vector<float> distances_from_mid = 	Geodesic_dijkstra(*m, m->midpoint);
	std::vector<float> voronoi_lengths(m->calculated_symmetry_pairs.size(), 0);

	std::vector<std::vector<int>> matching_paths; 

	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		std::vector<int> path = Geodesic_between_two_points(*m, index1, index2);
		matching_paths.push_back(path);
	}
	float biggest_size = 0;
	float smallest_size = INFINITY;
	for (size_t i = 0; i < matching_paths.size(); i++)
	{
		if (matching_paths[i].size() > biggest_size)
		{
			biggest_size = matching_paths[i].size();
		}
		if (matching_paths[i].size() < smallest_size)
		{
			smallest_size = matching_paths[i].size();
		}
	}
	for (size_t i = 0; i < weights.size(); i++)
	{
		float length_normalized = (matching_paths[i].size() - smallest_size) / (biggest_size - smallest_size);
		weights[i] = length_normalized;
		weights[i] = std::pow(weights[i],1);
	}

	//trying new length 
	std::vector<float> distances_from_midpoint = Geodesic_dijkstra(*m, m->midpoint);
	std::vector<float> midpoint_distances_from(m->calculated_symmetry_pairs.size(), 0);
	for (size_t i = 0; i < weights.size(); i++)
	{
		unsigned int midpoint = Geodesic_find_midpoint(m,m->calculated_symmetry_pairs[i].first, m->calculated_symmetry_pairs[i].second);
		float distance = distances_from_midpoint[midpoint];
		midpoint_distances_from[i] = distance; 
	}
	
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		unsigned int midpoint = Geodesic_find_midpoint(m ,index1,index2 );
		midpoints.push_back(std::make_pair(0, midpoint));
	}
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		Voronoi v(m, index1, index2, voronoi_param);
		v.generate_voronoi_parts();
		v.connect_boundary();
		v.generate_voronoi_parts();
		float part1 = 0;
		float part2 = 0;
		bool is_bad_voronoi = false; 

		//areas of both sides
		for (size_t i = 0; i < v.status.size(); i++)
		{
			if (v.status[i] == 1)
			{
				part1 = part1 + m->areas[i];
			}
			if (v.status[i] == 2)
			{
				part2 = part2 + m->areas[i];
			}
		}
		float area_ratio = std::min(part1, part2) / std::max(part1, part2);
		std::cout << i << " area ratio  " << area_ratio << std::endl;
		if (is_bad_voronoi || area_ratio < 0.9)
		{
			total_distances[i] = INFINITY;
			continue; 
		}

		voronoi_lengths[i] = v.indices.size();
		for (size_t j = 0; j < m->calculated_symmetry_pairs.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			unsigned int midpoint = midpoints[j].second;
			std::vector<float> distances = Geodesic_dijkstra(*m, midpoint);
			float min_dist = INFINITY;
			float dist_to_voronoi_1 = v.distance_to_index(m->calculated_symmetry_pairs[j].first);
			float dist_to_voronoi_2 = v.distance_to_index(m->calculated_symmetry_pairs[j].second);
			float ratio = std::min(dist_to_voronoi_1, dist_to_voronoi_2) / std::max(dist_to_voronoi_1, dist_to_voronoi_2);

			//check the best path voronoi intersects this
			std::vector<float> distances1 = Geodesic_dijkstra(*m, m->calculated_symmetry_pairs[j].first);
			std::vector<float> distances2 = Geodesic_dijkstra(*m, m->calculated_symmetry_pairs[j].second);

			float best_ratio = 0; 
			for (size_t k = 0; k < v.indices.size(); k++)
			{
				for (size_t t = 0; t < matching_paths[j].size(); t++)
				{
					if (v.indices[k] == matching_paths[j][t])
					{
						float dist_to_index1 = distances1[v.indices[k]];
						float dist_to_index2 = distances2[v.indices[k]];
						float ratio = std::min(dist_to_index1 , dist_to_index2) / std::max(dist_to_index1, dist_to_index2);
						if (ratio > best_ratio)
						{
							best_ratio = ratio; 
						}
					}
				}
			}

			for (size_t k = 0; k < v.indices.size(); k++)
			{
				int voronoi_index = v.indices[k];
				float voronoi_dist = distances[voronoi_index];
				if (voronoi_dist < min_dist)
				{
					min_dist =  voronoi_dist;
				}
			}
			//total_distances[i] += min_dist * weights[j];
			total_distances[i] += best_ratio * (1 - area_ratio); //* weights[j];
		}
		
	}

	auto least_elem_auto = std::min_element(total_distances.begin(), total_distances.end());
	int minIndex = std::distance(total_distances.begin(), least_elem_auto);

	for (size_t i = 0; i < total_distances.size(); i++)
	{
		std::cout << i << " " << total_distances[i] << std::endl;
	}
	std::cout << " selected voronoi pair " << m->calculated_symmetry_pairs[minIndex].first << " "
    << m->calculated_symmetry_pairs[minIndex].second << " " << std::endl;

	int index1 = m->calculated_symmetry_pairs[minIndex].first;
	int index2 = m->calculated_symmetry_pairs[minIndex].second;
	Voronoi v(m, index1, index2, voronoi_param);
	v.generate_voronoi_parts();
	v.connect_boundary();
	v.generate_voronoi_parts();
	v.color();

	// area comparison
	v.area_compare();
	return v;
}

void Voronoi_prune_voronoi(TrilateralMesh* m, Voronoi& voronoi, float voronoi_param)
{
	if (voronoi.indices.size() > m->vertices.size() / 2)
	{
		return; 
	}


	unsigned int midpoint = Geodesic_get_midpoint_from_path(m, voronoi.p1, voronoi.p2);
	unsigned  int voronoi_midpoint_inverse = Geodesic_send_ray_get_counterpart(m, midpoint);
	std::vector<float> distances_from_mid = Geodesic_dijkstra(*m, midpoint);

	float distance_to_mesh_midpoint = distances_from_mid[m->midpoint];

	std::vector<unsigned int> to_removed; 
	for (int i = voronoi.indices.size()-1; i >= 0  ; i--)
	{
		int index = voronoi.indices[i];
		if (distances_from_mid[index] > distance_to_mesh_midpoint)
		{
			to_removed.push_back(i);
		}
	}
	std::vector<unsigned int> re_colored;
	for (size_t i = 0; i < to_removed.size(); i++)
	{
		re_colored.push_back(voronoi.indices[to_removed[i]]);
		voronoi.indices.erase(voronoi.indices.begin() + to_removed[i]);
	}
	//midpoints
	std::vector<unsigned int> correspondence_midpoints;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		unsigned int midpoint = Geodesic_get_midpoint_from_path(m, index1, index2);
		correspondence_midpoints.push_back(midpoint);
	}

	std::vector<unsigned int> voronoi_front;
	std::vector<unsigned int> voronoi_back;
	unsigned int normal_midpoint_inverse = Geodesic_send_ray_get_counterpart(m, m->midpoint);
	glm::vec3 normal_midpoint_inv = m->normals[normal_midpoint_inverse];
	glm::vec3 normal_midpoint = m->normals[m->midpoint];
	for (size_t i = 0; i < voronoi.indices.size(); i++)
	{
		int index = voronoi.indices[i];
		glm::vec3 normal = m->normals[index];
		float dotfront = glm::dot(normal, normal_midpoint);
		float dotback = glm::dot(normal, normal_midpoint_inv);
		if (dotfront > dotback)
		{
			voronoi_front.push_back(index);
		}
		else
		{
			voronoi_back.push_back(index);
		}
	}

	float farthest_distance_front = 0;
	int farthest_midpoint_front = -1;
	for (size_t i = 0; i < correspondence_midpoints.size(); i++)
	{
		glm::vec3 normal = m->normals[correspondence_midpoints[i]];
		float dotfront = glm::dot(normal, normal_midpoint);
		float dotback = glm::dot(normal, normal_midpoint_inv);
		if (!(dotfront > dotback))
		{
			continue; 
		}
		float distance = 0;
		std::vector<float> distances = Geodesic_dijkstra(*m, correspondence_midpoints[i]);
		for (size_t j = 0; j < voronoi_front.size(); j++)
		{
			int index = voronoi_front[j];
			distance = distance + distances[index];
		}

		if (distance > farthest_distance_front)
		{
			farthest_distance_front = distance;
			farthest_midpoint_front = i;
		}
	}
	/*std::vector<int> path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_front], m->midpoint);
	std::vector<unsigned int> path_us;
	for (size_t i = 0; i < path.size(); i++)
	{
		path_us.push_back(path[i]);
	}
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
	} */

	float farthest_distance_back = 0;
	int farthest_midpoint_back = -1;
	for (size_t i = 0; i < correspondence_midpoints.size(); i++)
	{
		glm::vec3 normal = m->normals[correspondence_midpoints[i]];
		float dotfront = glm::dot(normal, normal_midpoint);
		float dotback = glm::dot(normal, normal_midpoint_inv);
		if (!(dotback > dotfront))
		{
			continue;
		}

		float distance = 0;
		std::vector<float> distances = Geodesic_dijkstra(*m, correspondence_midpoints[i]);
		for (size_t j = 0; j < voronoi_back.size(); j++)
		{
			int index = voronoi_back[j];
			distance = distance + distances[index];
		}

		if (distance > farthest_distance_back)
		{
			farthest_distance_back = distance;
			farthest_midpoint_back = i;
		}
	}
	m->color_all(WHITE);
	m->color_points(re_colored, WHITE);
	m->color_points(voronoi_back, MAGENTA);
	m->color_points(voronoi_front, BLUE);
	return;
	float closest_dist = INFINITY;
	int closest_index = -1;
	std::vector<float> distances_from_inverse_voronoi = Geodesic_dijkstra(*m , voronoi_midpoint_inverse); 
	for (size_t i = 0; i < voronoi.indices.size(); i++)
	{
		int index = voronoi.indices[i];
		float dist = distances_from_inverse_voronoi[index];
		if (closest_dist > dist)
		{
			closest_dist = dist; 
			closest_index = voronoi.indices[i];
		}
	}
	voronoi_midpoint_inverse = closest_index;

	std::vector<float> distances_from_voroni_midpoint = Geodesic_dijkstra(*m, midpoint);


	m->color_all(WHITE);
	std::vector<unsigned int > new_path = voronoi.get_closest_path(midpoint, voronoi_midpoint_inverse);
	m->color_points(new_path, YELLOW);
	bool is_midpoint_front_bad = true; 
	//if (distances_from_mid[correspondence_midpoints[farthest_midpoint_back]] > distances_from_mid[normal_midpoint_inverse])
	{
		std::vector<int> path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_back], normal_midpoint_inverse);

		m->color_points(conv_int_to_unsigned(path), DARKBLUE);
		for (size_t i = 0; i < path.size(); i++)
		{
			voronoi.indices.push_back(path[i]);
			new_path.push_back(path[i]);
		}
		is_midpoint_front_bad = false; 
	}

	bool is_midpoint_back_bad = true;
	//if (distances_from_mid[correspondence_midpoints[farthest_midpoint_front]] > distances_from_mid[m->midpoint])
	{
		std::vector<int> path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_front], m->midpoint);
		for (size_t i = 0; i < path.size(); i++)
		{
			voronoi.indices.push_back(path[i]);
			new_path.push_back(path[i]);
		}
		m->color_points(conv_int_to_unsigned(path), PURPLE);
		is_midpoint_back_bad = false;
	}
	
	std::vector<int> path = Geodesic_between_two_points(*m, voronoi_midpoint_inverse, m->midpoint);
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
		new_path.push_back(path[i]);
		m->color_points(conv_int_to_unsigned(path), GREEN);
	}
	path = Geodesic_between_two_points(*m, midpoint, normal_midpoint_inverse);
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
		new_path.push_back(path[i]);
		m->color_points(conv_int_to_unsigned(path), RED);
	}
	path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_front], normal_midpoint_inverse);
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
		new_path.push_back(path[i]);
		m->color_points(conv_int_to_unsigned(path), BLUE);
	}
	voronoi.indices.clear();
	voronoi.indices = new_path;
	voronoi.generate_voronoi_parts();
	voronoi.color(); 
	voronoi.area_compare();

	//m->color_points(new_path, ORANGE);

	/*std::vector<int> path1 = Geodesic_between_two_points(*m, midpoint, correspondence_midpoints[farthest_midpoint_front]);
	std::vector<int> path2 = Geodesic_between_two_points(*m, midpoint, correspondence_midpoints[farthest_midpoint_back]);
	std::vector<int> path3 = Geodesic_between_two_points(*m, midpoint, m->midpoint);
	std::vector<int> path4 = Geodesic_between_two_points(*m, midpoint, normal_midpoint_inverse);


	std::vector<int> path2 = Geodesic_between_two_points(*m, m->midpoint, correspondence_midpoints[farthest_midpoint_front]);
	std::vector<int> path3 = Geodesic_between_two_points(*m, midpoint, correspondence_midpoints[farthest_midpoint_front]);*/
	/*std::vector<int> path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_back], normal_midpoint_inverse);
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
	}
	path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_front], midpoint);
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
	}
	path = Geodesic_between_two_points(*m, correspondence_midpoints[farthest_midpoint_front], correspondence_midpoints[farthest_midpoint_back]);
	for (size_t i = 0; i < path.size(); i++)
	{
		voronoi.indices.push_back(path[i]);
	} 

	voronoi.generate_voronoi_parts(); */
	//voronoi.connect_boundary();
	//voronoi.color();


	//m->color_points(path_us, ORANGE);
	//m->color_points(path_us, ORANGE);


	//m->color_points(voronoi.indices, YELLOW);
}
Voronoi Voronoi_destroy_wrong_matches_and_recalculate(TrilateralMesh* m, float voronoi_param, Voronoi& v)
{
	std::vector<unsigned int> to_removed; 
	for (size_t i = m->calculated_symmetry_pairs.size()-1; i >0 ; i--)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		if (v.status[index1] == v.status[index2] && v.status[index1] != -1)
		{
			to_removed.push_back(i);
		}
	}
	for (size_t i = 0; i < to_removed.size(); i++)
	{
		m->calculated_symmetry_pairs.erase(m->calculated_symmetry_pairs.begin() + to_removed[i]);
	}
	//v = Voronoi_get_closest_voronoi(m, voronoi_param);
	return v; 
}
Voronoi Voronoi_check_pair_closeness_and_recalculate(TrilateralMesh* m, float voronoi_param, float dist_to_voronoi_param, Voronoi& v)
{
	std::vector<unsigned int > to_removed;
	for (size_t i = m->calculated_symmetry_pairs.size() - 1; i > 0; i--)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;

		std::vector<float> ratios; 

		std::vector<float> distances1 = Geodesic_dijkstra(*m, index1);
		std::vector<float> distances2 = Geodesic_dijkstra(*m, index2);

		for (size_t j = 0; j <  v.indices.size(); j++) // get 10 distances from indcies
		{
			float dist1 = distances1[v.indices[j]];
			float dist2 = distances2[v.indices[j]];
			float ratio = std::min(dist1, dist2) / std::max(dist1, dist2);
			ratios.push_back(ratio);
			
		}
		float avg_ratio = 0;
		for (size_t j = 0; j < ratios.size(); j++)
		{
			avg_ratio = avg_ratio + ratios[j];
		}
		avg_ratio = avg_ratio / ratios.size();
		if (avg_ratio < dist_to_voronoi_param)
		{
			to_removed.push_back(i);
		}
	}
	for (size_t i = 0; i < to_removed.size(); i++)
	{
		m->calculated_symmetry_pairs.erase(m->calculated_symmetry_pairs.begin() + to_removed[i]);
	}
	//v = Voronoi_get_closest_voronoi(m, voronoi_param);
	return v;
}

void Voronoi_color_every_pairs_voronoi(TrilateralMesh* m, float voronoi_param)
{
	std::vector<Voronoi> vec_v; 
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		Voronoi v(m, index1, index2, voronoi_param);
		v.generate_voronoi_parts();
		vec_v.push_back(v);

	}
	std::vector<unsigned int> points_mentioned(m->vertices.size(), 0);
	for (size_t i = 0; i < vec_v.size(); i++)
	{
		for (size_t j = 0; j < vec_v[i].indices.size(); j++)
		{
			int index = vec_v[i].indices[j];
			points_mentioned[index]++;
		}
	}

	float max = *std::max_element(points_mentioned.begin() , points_mentioned.end());
	float mean = 0;
	for (size_t i = 0; i < points_mentioned.size(); i++)
	{
		mean = mean + points_mentioned[i];
	}
	int size_without_0 = 0;
	for (size_t i = 0; i < points_mentioned.size(); i++)
	{
		if (points_mentioned[i] > 0)
		{
			size_without_0 += points_mentioned[i];
		}
	}

	mean = mean / size_without_0; // / N 

	float variance = 0;
	for (size_t i = 0; i < points_mentioned.size(); i++)
	{
		float dif = points_mentioned[i] - mean; 
		variance = dif * dif; 
	}
	variance = variance / points_mentioned.size();

	float std_dev = std::sqrt(variance);
	std_dev = 1; 

	std::vector<unsigned int> colored_vec; 
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (points_mentioned[i] >= mean - 1 && points_mentioned[i] <= mean + 1)
		{
			colored_vec.push_back(i);
		}
	}
	m->color_all(WHITE);
	m->color_points(colored_vec, BLUE);
}

Voronoi Voronoi_show_voronoi(TrilateralMesh* m,  unsigned int pair_no , float voronoi_param)
{
	int index1 = m->calculated_symmetry_pairs[pair_no].first;
	int index2 = m->calculated_symmetry_pairs[pair_no].second;
	Voronoi v(m, index1, index2, voronoi_param);
	v.generate_voronoi_parts();
	v.color();
	std::vector<unsigned int> pairs; 
	pairs.push_back(index1);
	pairs.push_back(index2);
	std::vector<float> distance = Geodesic_dijkstra(*m, index1);
	std::cout << " distance " << distance[index2] << std::endl;
	std::cout << index1 << " " << index2 << std::endl; 
	m->color_points(pairs, PURPLE);
	v.area_compare();
	return v;
}
int Voronoi_get_closest_unvisited_index(TrilateralMesh* m, Voronoi& v, std::vector<int>& is_visited , int index)
{
	std::vector<float> distances = Geodesic_dijkstra(*m, index);
	float min = INFINITY; 
	int min_index = -1;
	for (size_t i = 0; i < v.indices.size(); i++)
	{
		int index = v.indices[i];
		float dist = distances[index];
		if (dist < min && !is_visited[index])
		{
			min = dist; 
			min_index = index; 
		}
	}
	return min_index;
}
void Voronoi::connect_boundary()
{
	//select a random point
	//start a breadth first search
	int no_of_connected_points = 0;
	std::vector<bool> is_on_voronoi(m->vertices.size(), false);
	for (size_t i = 0; i < this->indices.size(); i++)
	{
		is_on_voronoi[this->indices[i]] = true;
	}
	std::stack<int> stack;  // a stack consisting of indices
	// get the adjacencies
	std::vector<std::vector<std::pair<int, float>>> mesh_adjacencies = m->adjacenies;
	//lastly get a int array with  size of vertices in order to check if the vertex has been visited ( -1 edge , 0 not visisted , 1 visited) 
	std::vector<int>is_visited(m->vertices.size(), 0);
	//push our point to stack
	std::vector<unsigned int> end_result; 
	bool last_time = false; 
	unsigned int start_index = this->indices[0];
	while (1)
	{
		stack.push(start_index);
		while (!stack.empty())
		{
			int index = stack.top();
			stack.pop(); //vertex index popped from stack
			if (is_visited[index] == 0) //not visited
			{
				end_result.push_back(index);
				is_visited[index] = 1; // now te vertex has been visited
				no_of_connected_points++;
				// this region of loop assumes index is not edge, therefore add the adjacencies
				for (size_t i = 0; i < mesh_adjacencies[index].size(); i++) //process pairs 
				{
					int adjacent_index = mesh_adjacencies[index][i].first;
					if (is_on_voronoi[adjacent_index])
					{
						stack.push(adjacent_index);
					}
				}
			}
			if (is_visited[index] == -1) //do nothing 
			{
				;
			}
		}
		if (last_time)
		{
			break; 
		}
		start_index = Voronoi_get_closest_unvisited_index(this->m, *this, is_visited, end_result[end_result.size() - 1]);
		if (start_index == -1) //means the end
		{
			last_time = true;
			start_index = this->indices[0];
			for (size_t i = 0; i < is_visited.size(); i++)
			{
				is_visited[0] = false; 
			}

		}
		std::vector<int> path = Geodesic_between_two_points(*this->m, start_index, end_result[end_result.size() - 1]);
		for (size_t i = 0; i < path.size(); i++)
		{
			end_result.push_back(path[i]);
		}

		//make end result unique
		std::sort(end_result.begin(), end_result.end());
		for (size_t i = end_result.size() - 1; i > 0; i--)
		{
			if (end_result[i] == end_result[i - 1])
				end_result.erase(end_result.begin() + i);
		}

	}
	this->indices.clear();

	//make end result unique
	std::sort(end_result.begin(), end_result.end());
	for (size_t i = end_result.size()-1; i > 0  ; i--)
	{
		if(end_result[i] == end_result[i-1])
		end_result.erase(end_result.begin() + i);
	}
	for (size_t i = 0; i < end_result.size(); i++)
	{
		this->indices.push_back(end_result[i]);
	}
	m->color_points(end_result, ORANGE);
}

//index1 , index2 mesh indices not voronoi indices 
std::vector<unsigned int> Voronoi::get_closest_path(unsigned int index1, unsigned int index2)
{
	std::vector<int> vector_of_indices(this->m->vertices.size(), -1);
	std::vector<int> vector_of_indices_inverse(this->indices.size(), -1);
	TrilateralMesh temp;
	for (size_t i = 0; i < this->indices.size(); i++)
	{
		temp.vertices.push_back(this->m->vertices[this->indices[i]]);
		vector_of_indices[this->indices[i]] = i;
		vector_of_indices_inverse[i] = this->indices[i];
	}
	temp.adjacenies = std::vector<std::vector<std::pair<int, float>>>(this->indices.size());
	for (size_t i = 0; i < this->indices.size(); i++)
	{
		int index = this->indices[i];
		std::vector<std::pair<int, float>> adjacency = this->m->adjacenies[index];
		std::vector<std::pair<int, float>> new_adjacency;
		for (size_t j = 0; j < adjacency.size(); j++)
		{
			int adjacent_index = adjacency[j].first;
			float distance = adjacency[j].second;
			if (vector_of_indices[adjacent_index] != -1)
			{
				std::pair<int, float> pair = std::make_pair(vector_of_indices[adjacent_index],distance);
				new_adjacency.push_back(pair);
			}
		}
		temp.adjacenies[i] = new_adjacency;
	}

	std::vector<unsigned int> path = conv_int_to_unsigned(Geodesic_between_two_points(temp, vector_of_indices[index1], vector_of_indices[index2]));

	std::vector<unsigned int> converted_path;
	for (size_t i = 0; i < path.size(); i++)
	{
		int index = path[i];
		int converted_back = vector_of_indices_inverse[index];
		converted_path.push_back(converted_back);
	}
	return converted_path;
}

void Voronoi_get_two_pair_and_generate_a_path(TrilateralMesh* m, Voronoi& voronoi, int voronoi_p1, float fuzziness , int voronoi_division_no)
{
	
	std::vector<float> distances = Geodesic_dijkstra(*m, voronoi_p1);
	std::vector<float> distances_among_voronoi;
	for (size_t i = 0; i < voronoi.indices.size(); i++)
	{
		int index = voronoi.indices[i];
		distances_among_voronoi.push_back(distances[index]);
	}
	float step = (*std::max_element(distances_among_voronoi.begin(), distances_among_voronoi.end()) - *std::min_element(distances_among_voronoi.begin(), distances_among_voronoi.end())) / 5; // default 
	std::cout << " step " << step << std::endl; 
	std::cout << " max " << *std::max_element(distances_among_voronoi.begin(), distances_among_voronoi.end()) << std::endl;
	std::cout << " min " << *std::min_element(distances_among_voronoi.begin(), distances_among_voronoi.end()) << std::endl;

	int index = -1; 
	for (size_t i = 0; i < voronoi.indices.size(); i++)
	{
		float distance = distances_among_voronoi[i];
		std::cout << "distance " << distance << std::endl;
		distance = distance - *std::min_element(distances_among_voronoi.begin(), distances_among_voronoi.end()); //normalized
		distance = distance / step; 
		if ( distance < voronoi_division_no + 1 && distance > voronoi_division_no - 1)
		{
			index = voronoi.indices[i];
			break; 
		}
	}

	float biggest_dist = 0; //useless
	std::vector<unsigned int> indices = { (unsigned int)voronoi_p1 , (unsigned int)index };
	NLateralDescriptor desc = NLateral_generate_descriptor_w_midpoints(m, indices, fuzziness, biggest_dist);
	desc.create_histogram_HKS(m , 5 ,voronoi_p1 );
	//m->color_points(desc.triangles_inside, BLUE);


}

void Voronoi_connect_midpoints(TrilateralMesh* m, Voronoi& voronoi)
{
	std::vector<unsigned int> midpoints;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		unsigned int midpoint = Geodesic_get_midpoint_from_path(m, index1, index2);
		midpoints.push_back(midpoint);
	}

	std::vector<int> paths; 
	//start from 0'th
	std::vector<bool> is_used(midpoints.size(), false);
	int no_of_used = 1; 
	is_used[0] = true; 
	unsigned int selected_index = 0;
	int last_index = -1;
	while (no_of_used < midpoints.size())
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, midpoints[selected_index]);
		int closest = -1; 
		float closest_dist = INFINITY; 
		for (size_t i = 0; i < midpoints.size(); i++)
		{
			if (!is_used[i])
			{
				if (closest_dist > distances[midpoints[i]])
				{
					closest_dist = distances[midpoints[i]];
					closest = i;
					last_index = i; 
				}
			}
		}
		std::vector<int> path = Geodesic_between_two_points(*m, midpoints[selected_index], midpoints[closest]);
		is_used[closest] = true; 
		no_of_used++;
		selected_index = closest;
		for (size_t i = 0; i < path.size(); i++)
		{
			paths.push_back(path[i]);
		}
	}
	std::vector<int> last_path = Geodesic_between_two_points(*m, midpoints[0], midpoints[last_index]);
	for (size_t i = 0; i < last_path.size(); i++)
	{
		paths.push_back(last_path[i]);
	}

	std::vector<bool> selected_points(m->vertices.size(), false);
	for (size_t i = 0; i < paths.size(); i++)
	{
		selected_points[paths[i]] = true; 
	}
	unsigned int seed_point = 0;
	while(selected_points[seed_point])
	{
		seed_point++;
	}
	std::vector<int> is_visited(m->vertices.size() , OUTSIDE );
	for (size_t i = 0; i < paths.size(); i++)
	{
		is_visited[paths[i]] = EDGE;
	}
	std::vector<unsigned int> path1 = breadth_first_search(m,seed_point, is_visited);
	//get other part
	for (size_t i = 0; i < path1.size(); i++)
	{
		selected_points[path1[i]] = true;
	}
	seed_point = 0;
	while (selected_points[seed_point])
	{
		seed_point++;
	}
	std::vector<unsigned int> path2 = breadth_first_search(m, seed_point, is_visited);

	m->color_points(path1, RED);
	m->color_points(path2, GREEN);
	m->color_points(conv_int_to_unsigned(paths), YELLOW);


	float area_red = 0;
	float area_green = 0;
	for (size_t i = 0; i < path1.size(); i++)
	{
		int index = path1[i];
		area_red += m->areas[index];
	}
	for (size_t i = 0; i < path2.size(); i++)
	{
		int index = path2[i];
		area_green += m->areas[index];
	}
	std::cout << " area is == > " << std::min(area_green ,area_red) / std::max(area_green, area_red) << std::endl;
}
