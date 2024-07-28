#include "../Include/Sampling.h"
#include "../Include/Geodesic.h"
#include <stdlib.h>     /* srand, rand */
std::vector<unsigned int>  furthest_point_sampling(Mesh* m, int no_of_samples, bool is_points_colored )
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
		std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, sample_idx);
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
	if (is_points_colored)
	{
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
			if (is_sampled)
			{
				new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
			}
			else
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			}

		}
		m->colors = new_color_buffer;
	}
	

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
std::vector<unsigned int>  furthest_point_sampling_on_partial_points(Mesh* m, int no_of_samples, std::vector<unsigned int>& map_partial_to_original)
{
	// generate a new mesh from existing mesh
	Mesh partial_mesh = *m;
	// create a boolean vector in order to select partial indices for O(1) selection
	std::vector<bool> is_point_on_partial_mesh(m->vertices.size() , false);
	for (size_t i = 0; i < map_partial_to_original.size(); i++)
	{
		is_point_on_partial_mesh[map_partial_to_original[i]] = true;
	}
	std::vector<unsigned int> map_original_to_partial(m->vertices.size(), -1);
	for (size_t i = 0; i < map_partial_to_original.size(); i++)
	{
		map_original_to_partial[map_partial_to_original[i]] = i;
	}

	partial_mesh.vertices.clear();
	partial_mesh.adjacenies.clear();
	partial_mesh.colors.clear();
	for (size_t i = 0; i < map_partial_to_original.size(); i++)
	{
		partial_mesh.vertices.push_back(m->vertices[map_partial_to_original[i]]);
	}
	for (size_t i = 0; i < map_partial_to_original.size(); i++)
	{
		partial_mesh.adjacenies.push_back(std::vector<std::pair<int, float>>());
		for (size_t j = 0; j < m->adjacenies[map_partial_to_original[i]].size(); j++)
		{
			if (is_point_on_partial_mesh[m->adjacenies[map_partial_to_original[i]][j].first ])
			{
				std::pair<int, float> adjacency;
				adjacency.first = map_original_to_partial[m->adjacenies[map_partial_to_original[i]][j].first];
				adjacency.second = m->adjacenies[map_partial_to_original[i]][j].second;
				partial_mesh.adjacenies[i].push_back(adjacency);
			}
		}
		partial_mesh.colors.push_back(glm::vec3(0, 0, 0));
	}

	//now with the new mesh, do furthespoint search 
	std::vector<unsigned int>  partial_mesh_fps_points =  furthest_point_sampling(&partial_mesh, no_of_samples , false);
	//convert the points back 
	std::vector<unsigned int> fps_points_corrected;
	for (size_t i = 0; i < partial_mesh_fps_points.size(); i++)
	{
		unsigned int fps_index_for_partial_point = map_partial_to_original[partial_mesh_fps_points[i]];
		
		fps_points_corrected.push_back(fps_index_for_partial_point);
	}
	return fps_points_corrected;
}