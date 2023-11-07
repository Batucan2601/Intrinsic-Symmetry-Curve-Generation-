#include "../Include/Sampling.h"
#include <stdlib.h>     /* srand, rand */
std::vector<unsigned int>  furthest_point_sampling(Mesh* m, int no_of_samples)
{
	float* distance = new float[m->vertices.size()];
	int* sampled = new int[no_of_samples];
	//initialize distance
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		distance[i] = 1e6; 
	}
	for (size_t i = 0; i < no_of_samples; i++)
	{
		sampled[i] = -1;
	}
	srand(time(NULL));
	//sampling starts
	int sample_idx = rand() % no_of_samples;
	sampled[0] = sample_idx;
	for (size_t i = 1; i < no_of_samples; i++)
	{
		//update distances
		std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m, sample_idx);
		for (size_t j = 0; j < distance_matrix_p1.size(); j++)
		{
			if (distance_matrix_p1[j] < distance[j])
			{
				distance[j] = distance_matrix_p1[j];
			}
		}
		//sample next
		// get max
		float maxValue = -1;
		int maxIndex = -1;
		for (size_t j = 0; j < distance_matrix_p1.size(); j++)
		{
			if (maxValue < distance[j])
			{
				maxValue = distance[j];
				maxIndex = j; 
			}
		}
		sample_idx = maxIndex;
		sampled[i] = sample_idx;
	}
	std::vector<unsigned int> sampled_id_vector; 
	for (size_t i = 0; i < no_of_samples; i++)
	{
		sampled_id_vector.push_back(sampled[i]);
	}
	delete[] distance;
	delete[] sampled;


	// redrawing part
	std::vector<glm::vec3> new_color_buffer;
	for (size_t i = 0; i < m->colors.size(); i++)
	{
		bool is_sampled = false; 
		for (size_t j = 0; j < no_of_samples; j++)
		{
			if (sampled_id_vector[j] == i)
			{
				is_sampled = true; 
				break;
			}
		}
		if (is_sampled ) 
		{
			new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
		}
		else
		{
			new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
		}
		
	}
	m->colors = new_color_buffer;

	return sampled_id_vector;
}

std::vector<unsigned int>  random_symmetry_indices_sampling(Mesh* m, int no_of_samples)
{
	srand(time(NULL));
	std::vector<unsigned int> random_sym_pairs;
	std::vector<bool> is_sample_exists(m->vertices.size(), 0);
	for (size_t i = 0; i < no_of_samples / 2 ; i++)
	{

		int random_no;
		while (true)
		{
			random_no = rand() % (m->vertices.size());
			if (!is_sample_exists[random_no])
			{
				is_sample_exists[random_no] = true; 
				break; 
			}
		}
		for (size_t j = 0; j < m->symmetry_pairs.size(); j++)
		{
			if (m->symmetry_pairs[j].first == random_no )
			{
				random_sym_pairs.push_back(m->symmetry_pairs[j].first);
				random_sym_pairs.push_back(m->symmetry_pairs[j].second);
				is_sample_exists[m->symmetry_pairs[j].second] = true; 
			}
		}
	}

	m->colors.clear();
	// color
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (is_sample_exists[i])
		{
			m->colors.push_back(glm::vec3(1.0f, 1.0f, 0.0f));
		}
		else
		{
			m->colors.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
		}
	}

	return random_sym_pairs; 
}

//furthest points sampling in specific part of a mesh decided by vector partial points 
std::vector<unsigned int>  furthest_point_sampling_on_partial_points(Mesh* m, int no_of_samples, std::vector<unsigned int>& partial_points)
{
	// generate a new mesh from existing mesh
	Mesh partial_mesh = *m;
	// create a boolean vector in order to select partial indices for O(1) selection
	std::vector<bool> is_point_on_partial_mesh(m->vertices.size() , false);
	for (size_t i = 0; i < partial_points.size(); i++)
	{
		is_point_on_partial_mesh[partial_points[i]] = true;
	}
	partial_mesh.vertices.clear();
	partial_mesh.adjacenies.clear();
	partial_mesh.colors.clear();
	for (size_t i = 0; i < partial_points.size(); i++)
	{
		partial_mesh.vertices.push_back(m->vertices[partial_points[i]]);
	}
	for (size_t i = 0; i < partial_points.size(); i++)
	{
		for (size_t j = 0; j < m->adjacenies[partial_points[i]].size(); j++)
		{
			if (is_point_on_partial_mesh[m->adjacenies[partial_points[i]][j].first ])
			{
				std::pair<int, float> adjacency;
				
				//adjacency.first = 
				//partial_mesh.adjacenies.push_back()
			}
		}
	}
}