#include "../Include/Sampling.h"
#include <stdlib.h>     /* srand, rand */
std::vector<int>  furthest_point_sampling(Mesh* m, int no_of_samples)
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
	std::vector<int> sampled_id_vector; 
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
			}
		}
		if (is_sampled ) 
		{
			new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
		}
		else
		{
			new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
		}
		
	}
	m->colors = new_color_buffer;

	return sampled_id_vector;
}