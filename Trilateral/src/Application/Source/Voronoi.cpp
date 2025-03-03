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
Voronoi Voronoi_algorithm_in_action(TrilateralMesh* m, float voronoi_param, float hks_param, float distance_param, int hist_no, float fuziness, int no_of_points,
	std::vector<unsigned int>& agd_indices)
{
	Voronoi v = Voronoi_get_closest_voronoi(m, voronoi_param);
	v = Voronoi_destroy_wrong_matches_and_recalculate(m, voronoi_param, v);
	v = Voronoi_check_pair_closeness_and_recalculate(m, voronoi_param, distance_param, v);

	v = Voronoi_add_points_and_recalculate(m, v, voronoi_param, hks_param, distance_param, hist_no, fuziness, no_of_points, agd_indices);
	return v; 
}
Voronoi Voronoi_add_points_and_recalculate(TrilateralMesh* m, Voronoi& v,float voronoi_param, float hks_param, float distance_param, int hist_no, float fuziness, int no_of_points
,std::vector<unsigned int>& agd_points)
{
	for (size_t i = 0; i < no_of_points; i++)
	{
		Voronoi_add_point(m, v, voronoi_param,  hks_param, distance_param, m->midpoint, hist_no, fuziness , agd_points);
	}
	v = Voronoi_get_closest_voronoi(m, voronoi_param);
	return v; 
}
void Voronoi_add_point(TrilateralMesh* m, Voronoi& voronoi , float voronoi_param , float hks_param , float dist_param , unsigned int midpoint , int hist_no, float fuziness
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
	//select a never selected index
	unsigned int random_index = 0; 
	for (size_t i = 0; i < agd_points.size(); i++)
	{
		if (!is_in_matching[i])
		{
			random_index = agd_points[i];
			break;
		}
	}

	float hks = m->normalized_heat_kernel_signature[random_index];
	std::vector<float> distances = Geodesic_dijkstra(*m, m->midpoint);
	//select a counterpart
	std::vector<std::pair<float, unsigned int> > results; 
	for (size_t i = 0; i < agd_points.size(); i++)
	{
		unsigned int index = agd_points[i];
		if (is_in_matching[i] == true || voronoi.status[index] == -1 || i == random_index)
		{
			continue; 
		}
		float hks_i = m->normalized_heat_kernel_signature[index];
		float hks_dif = std::abs(hks - hks_i);
		float ratio = std::min(distances[random_index], distances[index]) / std::max(distances[random_index], distances[index]);
		if (hks_dif < hks_param && ratio < dist_param)
		{
			continue; 
		}
		NLateralDescriptor desc_j_i = NLateral_generate_symmetric_descriptor(m, index, i, hist_no, fuziness);
		NLateralDescriptor desc_i_j = NLateral_generate_symmetric_descriptor(m, i, index, hist_no, fuziness);
		float dif = VarianceMin_compare(m, desc_i_j, desc_j_i, true, hist_no, 1);
		results.push_back(std::make_pair(dif,index));
	}
	auto min_index_auto = std::min_element( results.begin(), results.end());
	unsigned int min_index = std::distance(results.begin() , min_index_auto );

	m->calculated_symmetry_pairs.push_back(std::make_pair(random_index, min_index));

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
		weights[i] = 1 ;
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
			for (size_t k = 0; k < v.indices.size(); k++)
			{
				int voronoi_index = v.indices[k];
				float voronoi_dist = distances[voronoi_index];
				if (voronoi_dist < min_dist)
				{
					min_dist =  voronoi_dist;
				}
			}
			total_distances[i] += min_dist * weights[j];
		}
	}
	auto least_elem_auto = std::min_element(total_distances.begin(), total_distances.end());
	int minIndex = std::distance(total_distances.begin(), least_elem_auto);
	std::cout << " min element" << m->calculated_symmetry_pairs[minIndex].first << " "
    << m->calculated_symmetry_pairs[minIndex].second << " " << std::endl;

	int index1 = m->calculated_symmetry_pairs[minIndex].first;
	int index2 = m->calculated_symmetry_pairs[minIndex].second;
	Voronoi v(m, index1, index2, voronoi_param);
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
		if (v.status[index1] == v.status[index2])
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
		float dist1 = v.distance_to_index(index1);
		float dist2 = v.distance_to_index(index2);
		
		float ratio = std::min(dist1,dist2) / std::max(dist1,dist2);
		if (ratio < dist_to_voronoi_param)
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