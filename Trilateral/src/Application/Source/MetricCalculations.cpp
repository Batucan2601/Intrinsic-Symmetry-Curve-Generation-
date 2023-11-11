#include "../Include/MetricCalculations.h"
#include "../Include/MeshFactory.h"

//v6'nin ground-truth simetrik noktasi v66 ise ve senin methodun v6->v77'ye gonderdiyse o zaman geodesic(v66, v77) costun olacak;
float get_geodesic_cost(Mesh* m, unsigned int point_index1, unsigned int calculated_index1_correspondence)
{
	unsigned int ground_truth = -1;
	//get symmetry pair
	for (size_t i = 0; i < m->symmetry_pairs.size(); i++)
	{
		if (m->symmetry_pairs[i].first == point_index1)
		{
			ground_truth = m->symmetry_pairs[i].second;
			break;
		}
		else if (m->symmetry_pairs[i].second == point_index1)
		{
			ground_truth = m->symmetry_pairs[i].first;
			break;
		}
	}
	std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m, ground_truth);

	float cost = distance_matrix_p1[calculated_index1_correspondence]; 

	return cost;

}
float get_geodesic_cost_with_list(Mesh* m, std::vector<unsigned int> point_indices, std::vector<unsigned int> calculated_index_correspondence_list)
{
	unsigned int N = point_indices.size();
	float error = 0;
	std::vector<unsigned int> ground_truth_indices; 
	for (size_t i = 0; i < N; i++)
	{
		float single_error = get_geodesic_cost( m , point_indices[i], calculated_index_correspondence_list[i]);
		error += single_error;
	}
	return error/N ; // return average error
}

