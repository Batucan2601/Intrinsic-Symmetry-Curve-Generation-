#include "../Include/Voronoi.h"
#include "../Include/Geodesic.h"
#include <stack>
#include <random>
#include "../Include/VarianceMinimizingTransportPlan.h"

Voronoi::Voronoi()
{
	return;
}
Voronoi::Voronoi(TrilateralMesh* m, unsigned int p1, unsigned int p2 , float param )
{
	this->m = m;
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
		if (this->status[i] == -1 )
		{
			voronoi_part.push_back(i);
		}
		else if (this->status[i] == 1)
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

	this->m->color_points(voronoi_part, YELLOW);
	this->m->color_points(half1, GREEN);
	this->m->color_points(half2, RED);
	this->m->color_points(half3, BLUE);
	this->m->color_points(half4, ORANGE);
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
	v = Voronoi_destroy_wrong_matches_and_recalculate(m, voronoi_param, v);
	v = Voronoi_check_pair_closeness_and_recalculate(m, voronoi_param, dist_param, v);

	v = Voronoi_add_points_and_recalculate(m, v, voronoi_param, hks_param, sdf_param, dist_param, hist_no, fuziness, no_of_points, agd_indices);
	v = Voronoi_get_closest_voronoi(m, voronoi_param);
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
	v = Voronoi_get_closest_voronoi(m, voronoi_param);
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
			for (size_t k = 0; k < voronoi.indices.size(); k++) // get 10 distances from indcies
			{
				std::vector<float> distances = Geodesic_dijkstra(*m, voronoi.indices[k]);
				float dist1 = distances[selected_agd_index];
				float dist2 = distances[index_j];
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

	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		int midpoint = Geodesic_get_midpoint_from_path(m, index1, index2);
		midpoints.push_back(std::make_pair(0,midpoint));
		//weights 
		std::vector<float> distances_from_mid = Geodesic_dijkstra(*m, midpoint);
		//weights.push_back((distances_from_mid[index1] + distances_from_mid[index2]) / 2);
		
		float dif_midpoint = std::abs(distances_from_mid[index1] + distances_from_mid[index2]);
		float dif_hks = std::abs(m->normalized_heat_kernel_signature[index1] - m->normalized_heat_kernel_signature[index2]);
		float sigma = 0.02; // random
		//float weight = std::exp(-(dif_midpoint * dif_midpoint) / (2 * sigma * sigma));
		float weight = std::exp(-(dif_midpoint * dif_midpoint) / (2 * sigma * sigma));
		//weights[i] = ( weight);
		Voronoi v(m, index1, index2, voronoi_param);
		std::vector<float> dist_index1 = Geodesic_dijkstra(*m, index1);
		std::vector<float> dist_index2 = Geodesic_dijkstra(*m, index2);
		weight = 0; 
		for (size_t i = 0; i < v.indices.size(); i++)
		{
			float dist1 = dist_index1[v.indices[i]];
			float dist2 = dist_index2[v.indices[i]];
			weight = weight + std::min(dist1,dist2) / std::max(dist1,dist2);
		}
		weight = weight / v.indices.size();
		weights[i] = weight ;
	}
	float max = *std::max_element(weights.begin() , weights.end());
	/*for (size_t i = 0; i < weights.size(); i++)
	{
		weights[i] = weights[i] / max;
	} */
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		Voronoi v(m, index1, index2, voronoi_param);
		for (size_t j = 0; j < m->calculated_symmetry_pairs.size(); j++)
		{
			unsigned int midpoint = midpoints[j].second;
			std::vector<float> distances = Geodesic_dijkstra(*m, midpoint);
			float min_dist = INFINITY;
			float dist_to_voronoi_1 = v.distance_to_index(m->calculated_symmetry_pairs[j].first);
			float dist_to_voronoi_2 = v.distance_to_index(m->calculated_symmetry_pairs[j].second);
			float ratio = std::min(dist_to_voronoi_1, dist_to_voronoi_2) / std::max(dist_to_voronoi_1, dist_to_voronoi_2);
			/*float sigma = 0.02; // random
			float weight = std::exp(-(1.0/ratio * 1.0/ratio) / (2 * sigma * sigma));
			weights[j] = weight;*/
			for (size_t k = 0; k < v.indices.size(); k++)
			{
				int voronoi_index = v.indices[k];
				float voronoi_dist = distances[voronoi_index];
				if (voronoi_dist < min_dist)
				{
					min_dist =  voronoi_dist;
				}
			}
			total_distances[i] += min_dist * std::pow(weights[j] , 4 );
		}
	}
	auto least_elem_auto = std::min_element(total_distances.begin(), total_distances.end());
	int minIndex = std::distance(total_distances.begin(), least_elem_auto);
	std::cout << " selected voronoi pair " << m->calculated_symmetry_pairs[minIndex].first << " "
    << m->calculated_symmetry_pairs[minIndex].second << " " << std::endl;

	int index1 = m->calculated_symmetry_pairs[minIndex].first;
	int index2 = m->calculated_symmetry_pairs[minIndex].second;
	Voronoi v(m, index1, index2, voronoi_param);
	v.generate_voronoi_parts();
	v.connect_boundary();
	v.generate_voronoi_parts();
	v.color();
	return v;
}
Voronoi Voronoi_destroy_wrong_matches_and_recalculate(TrilateralMesh* m, float voronoi_param, Voronoi& v)
{
	std::vector<unsigned int> to_removed; 
	for (size_t i = m->calculated_symmetry_pairs.size()-1; i >0 ; i--)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		if (v.status[index1] == v.status[index2] || v.status[index2] == -1 || v.status[index1] == -1)
		{
			to_removed.push_back(i);
		}
	}
	for (size_t i = 0; i < to_removed.size(); i++)
	{
		m->calculated_symmetry_pairs.erase(m->calculated_symmetry_pairs.begin() + to_removed[i]);
	}
	v = Voronoi_get_closest_voronoi(m, voronoi_param);
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
		for (size_t j = 0; j <  v.indices.size(); j++) // get 10 distances from indcies
		{
			std::vector<float> distances = Geodesic_dijkstra(*m, v.indices[j]);
			float dist1 = distances[index1];
			float dist2 = distances[index2];
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
	v = Voronoi_get_closest_voronoi(m, voronoi_param);
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

void Voronoi_show_voronoi(TrilateralMesh* m,  unsigned int pair_no , float voronoi_param)
{
	int index1 = m->calculated_symmetry_pairs[pair_no].first;
	int index2 = m->calculated_symmetry_pairs[pair_no].second;
	Voronoi v(m, index1, index2, voronoi_param);
	v.generate_voronoi_parts();
	v.color();
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
