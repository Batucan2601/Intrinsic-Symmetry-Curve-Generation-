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
	this->status = std::vector<int>(this->m->vertices.size() , 0 );
	for (size_t i = 0; i < this->indices.size(); i++)
	{
		this->status[this->indices[i]] = -1;
	}
	//pick a non curve point
	unsigned int random_pick = 0;
	while (this->status[random_pick] == -1)
	{
		random_pick++;
		if (random_pick == this->status.size())
		{
			return;
		}
	}
	std::vector<unsigned int> visited_vertices;
	//push our point to stack
	std::stack<int> stack;

	stack.push(random_pick);
	while (!stack.empty())
	{
		int index = stack.top();
		stack.pop(); //vertex index popped from stack
		if (this->status[index] == 0) //not visited
		{
			this->status[index] = 1; // now te vertex has been visited
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

	//now some are 1 some are 0 
}

void Voronoi::color()
{
	std::vector<unsigned int> voronoi_part; 
	std::vector<unsigned int> half1;
	std::vector<unsigned int> half2;
	for (  int i = 0; i < this->status.size();  i++)
	{
		if (this->status[i] == -1 )
		{
			voronoi_part.push_back(i);
		}
		else if (this->status[i] == 0)
		{
			half1.push_back(i);
		}
		else if (this->status[i] == 1)
		{
			half2.push_back(i);
		}
	}

	this->m->color_points(voronoi_part, YELLOW);
	this->m->color_points(half1, GREEN);
	this->m->color_points(half2, RED);
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
	float total = 0;
	for (size_t i = 0; i < no_of_point_param; i++)
	{
		unsigned int index = dist(gen);
		index = v.indices[index];
		NLateralDescriptor desc1 = NLateral_generate_symmetric_descriptor(m, p1, index, hist_no, fuziness);
		NLateralDescriptor desc2 = NLateral_generate_symmetric_descriptor(m, p2, index, hist_no, fuziness);
		float dif = VarianceMin_compare(m, desc1, desc2, true, hist_no, 1);
		total = total + dif; 
	}

	return total;
}

Voronoi Voronoi_get_closest_voronoi(TrilateralMesh* m, float voronoi_param)
{
	std::vector<std::pair<float,unsigned int>> midpoints;
	std::vector<float> total_distances(m->calculated_symmetry_pairs.size(),0);
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		int midpoint = Geodesic_get_midpoint_from_path(m, index1, index2);
		midpoints.push_back(std::make_pair(0,midpoint));
	}
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
					min_dist = voronoi_dist;
				}
			}
			total_distances[i] += min_dist;
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