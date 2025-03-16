#include "../Include/Geodesic.h"
#include <unordered_set>
#include <queue>
#include <unordered_map>
#include "../Include/CoreTypeDefs.h"
#include "raylib.h"
#include "../Include/Ray.h"
#include <stack>
#include <random>
#include "../Include/Sampling.h"
// Function to find n-ring neighbors for a specific point
static std::vector<unsigned int> findNRingNeighbors(TrilateralMesh* m, int startPointId, int n) {
	std::unordered_set<int> visited; // Track visited points
	std::queue<std::pair<int, int>> q; // Pair of (current point ID, depth level)
	std::unordered_set<int> result; // Store neighbors within n-rings

	q.push({ startPointId, 0 });
	visited.insert(startPointId);

	while (!q.empty()) {
		std::pair<int,int>p = q.front();
		int depth = p.second;
		q.pop();

		if (depth > 0) {
			result.insert(p.first);
		}

		if (depth < n) {

			for (int i = 0; i < m->adjacenies[p.first].size(); i++ ) {
				int neighbour = m->adjacenies[p.first][i].first;
				if (visited.find(neighbour) == visited.end()) {
					visited.insert(neighbour);
					q.push({ neighbour, depth + 1 });
				}
			}
		}
	}
	std::vector<unsigned int> return_vec; 
	for (auto& num : result) {
		return_vec.push_back(num);
	}
	return return_vec;
}

std::vector<float> Geodesic_dijkstra(TrilateralMesh& m, int point_index)
{
	// return for prdecessors
	std::vector<float> matrix_vec;
	//init fibonacci heap
	FibonacciHeap<std::pair<float, int >> min_fib_heap;
	//init array 
	float* matrix = new float[m.vertices.size()];
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if (point_index == i)
		{
			matrix[i] = 0.0f;
		}
		else
		{
			//INFINITY
			matrix[i] = (float)INFINITY;
		}
	}
	//init extra stuff
	int vertex_visited = 0; //total umber of visits 
	unsigned int* is_vertex_visited = new unsigned int[m.vertices.size()];
	int* predecessor = new int[m.vertices.size()]; //predecessor vertices for all 
	//float* will_be_visited_list = new float[m.vertices.size()];
	// std::vector<std::pair<int , float>> shortest_path_set; //not yet included but calculated 
	//fill array 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_vertex_visited[i] = 0; // not visited
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessor[i] = -1; // not visited
	}

	//dijkstra
	std::pair<float, int > first_index;
	first_index.first = 0.0f; //distance 
	first_index.second = point_index;

	std::pair<float, int > current_pair;
	min_fib_heap.insert(first_index);
	while (!min_fib_heap.isEmpty())
	{
		current_pair = min_fib_heap.removeMinimum();
		// get the adjacent vertices
		int current_index = current_pair.second;
		// if it is already visited pass
		if (is_vertex_visited[current_index] == 1)
		{
			continue;
		}
		is_vertex_visited[current_index] = 1;
		float distance = current_pair.first;
		for (size_t i = 0; i < m.adjacenies[current_index].size(); i++)
		{
			int adjacent_index = m.adjacenies[current_index][i].first;
			float adjacent_distance = m.adjacenies[current_index][i].second;
			if (is_vertex_visited[adjacent_index] == 0)
			{
				//relaxation
				if (matrix[adjacent_index] > distance + m.adjacenies[current_index][i].second)
				{
					//update
					matrix[adjacent_index] = distance + m.adjacenies[current_index][i].second;
					predecessor[adjacent_index] = current_index;
					//push 
					std::pair<float, int> new_vertex;
					new_vertex.first = matrix[adjacent_index];
					new_vertex.second = adjacent_index;
					min_fib_heap.insert(new_vertex);
				}
			}
		}
	}
	// copy the predecessor 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		matrix_vec.push_back(matrix[i]);
	}
	delete[] matrix;
	delete[] is_vertex_visited;
	delete[] predecessor;
	return matrix_vec;

}
std::vector<int> Geodesic_dijkstra_predecessors(TrilateralMesh& m, int point_index)
{
	// return for prdecessors
	std::vector<int> predecessors;
	//init fibonacci heap
	FibonacciHeap<std::pair<float, int >> min_fib_heap;
	//init array 
	float* matrix = new float[m.vertices.size()];
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if (point_index == i)
		{
			matrix[i] = 0.0f;
		}
		else
		{
			//INFINITY
			matrix[i] = (float)INFINITY;
		}
	}
	//init extra stuff
	int vertex_visited = 0; //total umber of visits 
	unsigned int* is_vertex_visited = new unsigned int[m.vertices.size()];
	int* predecessor = new int[m.vertices.size()]; //predecessor vertices for all 
	//float* will_be_visited_list = new float[m.vertices.size()];
	// std::vector<std::pair<int , float>> shortest_path_set; //not yet included but calculated 
	//fill array 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_vertex_visited[i] = 0; // not visited
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessor[i] = -1; // not visited
	}

	//dijkstra
	std::pair<float, int > first_index;
	first_index.first = 0.0f; //distance 
	first_index.second = point_index;

	std::pair<float, int > current_pair;
	min_fib_heap.insert(first_index);
	while (!min_fib_heap.isEmpty())
	{
		current_pair = min_fib_heap.removeMinimum();
		// get the adjacent vertices
		int current_index = current_pair.second;
		// if it is already visited pass
		if (is_vertex_visited[current_index] == 1)
		{
			continue;
		}
		is_vertex_visited[current_index] = 1;
		float distance = current_pair.first;
		for (size_t i = 0; i < m.adjacenies[current_index].size(); i++)
		{
			int adjacent_index = m.adjacenies[current_index][i].first;
			float adjacent_distance = m.adjacenies[current_index][i].second;
			if (is_vertex_visited[adjacent_index] == 0)
			{
				//relaxation
				if (matrix[adjacent_index] > distance + m.adjacenies[current_index][i].second)
				{
					//update
					matrix[adjacent_index] = distance + m.adjacenies[current_index][i].second;
					predecessor[adjacent_index] = current_index;
					//push 
					std::pair<float, int> new_vertex;
					new_vertex.first = matrix[adjacent_index];
					new_vertex.second = adjacent_index;
					min_fib_heap.insert(new_vertex);
				}
			}
		}
	}
	// copy the predecessor 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessors.push_back(predecessor[i]);
	}
	delete[] matrix;
	delete[] is_vertex_visited;
	delete[] predecessor;
	return predecessors;

}
std::vector<int> Geodesic_between_two_points(TrilateralMesh& m, int p1_index, int p2_index)
{
	std::vector<int> predecessors = Geodesic_dijkstra_predecessors(m, p1_index);
	std::vector<int> consec_indices;
	int pred = p2_index;
	while (pred != p1_index)
	{
		consec_indices.push_back(pred);
		pred = predecessors[pred];
	}
	consec_indices.insert(consec_indices.end(), p1_index);
	return consec_indices;
}

bool Geodesic_proximity(TrilateralMesh& m , NLateralDescriptor& desc1, NLateralDescriptor& desc2, float proximity)
{
	std::vector<float> distances = Geodesic_dijkstra(m, desc1.indices[0]);
	return distances[desc2.indices[0]] < proximity;
}

std::vector<unsigned int> Geodesic_avg_dijkstra_modified(TrilateralMesh* m, float sweep_percentage, int N_ring, bool is_color , float& biggest_dijkstra)
{
	std::vector<std::pair<float, unsigned int>> gradient_indices;
	int vertex_size = m->vertices.size();
	std::vector<float> avg_geodesic_distances(vertex_size);
	float biggest = -INFINITY;
	biggest_dijkstra = -INFINITY;
	for (size_t i = 0; i < vertex_size; i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, i);
		//sum the distances
		float sum = 0;
		for (size_t j = 0; j < vertex_size; j++)
		{
			sum += distances_i[j];
			if (biggest_dijkstra < distances_i[j])
			{
				biggest_dijkstra = distances_i[j];
			}
		}
		avg_geodesic_distances[i] = sum;
		if (sum > biggest)
		{
			biggest = sum;
		}
	}
	//average the one ring
	std::vector<float> temp_avg_values = avg_geodesic_distances;
	for (size_t i = 0; i < vertex_size; i++)
	{
		float avg = temp_avg_values[i];
		std::vector<unsigned int> neighbours = findNRingNeighbors(m, i, N_ring);
		for (size_t j = 0; j < neighbours.size(); j++)
		{
			int neighbour_index = neighbours[j];
			avg += temp_avg_values[neighbour_index];
		}
		avg = avg / (neighbours.size() + 1);
		avg_geodesic_distances[i] = avg;
	}
	//sample where gradient is 0
	for (size_t i = 0; i < vertex_size; i++)
	{
		bool is_maxima = true;
		bool is_minima = true;

		std::vector<unsigned int> neighbours;
		neighbours = findNRingNeighbors(m, i, N_ring);
		for (size_t j = 0; j < neighbours.size(); j++)
		{
			int index = neighbours[j];
			if (index == i)
			{
				continue;
			}
			if (avg_geodesic_distances[index] > avg_geodesic_distances[i])
			{
				is_maxima = false;
			}
			if (avg_geodesic_distances[index] < avg_geodesic_distances[i])
			{
				is_minima = false;
			}
		}
		if (is_maxima || is_minima)
		{
			gradient_indices.push_back(std::make_pair(avg_geodesic_distances[i], i));
		}
	}
	std::vector< unsigned int> extremums;
	for (size_t i = 0; i < gradient_indices.size(); i++)
	{
		int index_front = gradient_indices[i].second;
		bool is_index_front = true;
		for (size_t j = 0; j < extremums.size(); j++)
		{
			int ext_index = extremums[j];
			std::vector<float> distances_ext = Geodesic_dijkstra(*m, ext_index);
			if (distances_ext[index_front] < (biggest_dijkstra * sweep_percentage))
			{
				is_index_front = false;
			}

		}
		if (is_index_front)
		{
			extremums.push_back(gradient_indices[i].second);
		}
	}
	std::vector<unsigned int> extremum_simplified; 
	for (size_t i = 0; i < extremums.size(); i++)
	{
		int index_i = extremums[i];
		float val_i = avg_geodesic_distances[index_i];
		float percentage = 5;
		bool is_close_found = false; 
		for (size_t j = 0; j < extremums.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			int index_j = extremums[j];
			float val_j = avg_geodesic_distances[index_j];
			float dif = std::abs(val_i - val_j);
			if (dif < (val_i * percentage / 100 ))
			{
				is_close_found = true; 
				break; 
			}
		}
		if (is_close_found)
		{
			extremum_simplified.push_back(index_i);
		}
	}
	extremums = extremum_simplified;
	if (is_color)
	{
		m->color_all(BLACK);
		//min max float 
		auto min = std::min_element(avg_geodesic_distances.begin() , avg_geodesic_distances.end());
		auto max = std::max_element(avg_geodesic_distances.begin() , avg_geodesic_distances.end());
		for (size_t i = 0; i < vertex_size; i++)
		{
			glm::vec3 color = CoreType_getColor(avg_geodesic_distances[i], *min, *max);
			m->raylib_mesh.colors[i * 4] = color[0] * 255;
			m->raylib_mesh.colors[i * 4 + 1] = color[1] * 255;
			m->raylib_mesh.colors[i * 4 + 2] = color[2] * 255;
			m->raylib_mesh.colors[i * 4 + 3] = 255;
		}
		/*for (size_t i = 0; i < extremums.size(); i++)
		{
			int index = extremums[i];
			m->raylib_mesh.colors[index * 4] = 255;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}*/
		m->update_raylib_mesh();
	}

	return extremums;

}
std::vector<unsigned int> Geodesic_avg_dijkstra_modified_to_points(TrilateralMesh* m, std::vector<unsigned int> points, int& no_of_points, 
float sweep_percentage, int N_ring, bool is_color)
{
	std::vector<std::pair<float, unsigned int>> gradient_indices;
	int vertex_size = m->vertices.size();
	std::vector<float> avg_geodesic_distances(vertex_size);
	float biggest = -INFINITY;
	float biggest_dijkstra = -INFINITY;
	for (size_t i = 0; i < vertex_size; i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, i);

		//sum the distances
		float sum = 0;
		for (size_t j = 0; j < points.size(); j++)
		{
			int index = points[j];
			sum += distances_i[index];
			if (biggest_dijkstra < distances_i[index])
			{
				biggest_dijkstra = distances_i[index];
			}
		}
		avg_geodesic_distances[i] = sum;
		if (sum > biggest)
		{
			biggest = sum;
		}

	}

	//average the one ring
	std::vector<float> temp_avg_values = avg_geodesic_distances;
	for (size_t i = 0; i < vertex_size; i++)
	{
		float avg = temp_avg_values[i];
		std::vector<unsigned int> neighbours = findNRingNeighbors(m, i, N_ring);
		for (size_t j = 0; j < neighbours.size(); j++)
		{
			int neighbour_index = neighbours[j];
			avg += temp_avg_values[neighbour_index];
		}
		avg = avg / (neighbours.size() + 1);
		avg_geodesic_distances[i] = avg;
	}
	//sample where gradient is 0
	for (size_t i = 0; i < vertex_size; i++)
	{
		bool is_maxima = true;
		bool is_minima = true;

		std::vector<unsigned int> neighbours;
		neighbours = findNRingNeighbors(m, i, N_ring);
		for (size_t j = 0; j < neighbours.size(); j++)
		{
			int index = neighbours[j];
			if (index == i)
			{
				continue;
			}
			if (avg_geodesic_distances[index] > avg_geodesic_distances[i])
			{
				is_maxima = false;
			}
			if (avg_geodesic_distances[index] < avg_geodesic_distances[i])
			{
				is_minima = false;
			}

		}
		/*int neighbour_size = m->adjacenies[i].size();
		for (size_t j = 0; j < neighbour_size; j++)
		{
			int neighbour_index = m->adjacenies[i][j].first;
			if (avg_geodesic_distances[neighbour_index] > avg_geodesic_distances[i])
			{
				is_maxima = false;
			}
			if (avg_geodesic_distances[neighbour_index] < avg_geodesic_distances[i])
			{
				is_minima = false;
			}
		}*/
		if (is_maxima || is_minima)
		{
			gradient_indices.push_back(std::make_pair(avg_geodesic_distances[i], i));
		}
	}
	std::vector< unsigned int> extremums;
	for (size_t i = 0; i < gradient_indices.size(); i++)
	{
		int index_front = gradient_indices[i].second;
		bool is_index_front = true;
		for (size_t j = 0; j < extremums.size(); j++)
		{
			int ext_index = extremums[j];
			std::vector<float> distances_ext = Geodesic_dijkstra(*m, ext_index);
			if (distances_ext[index_front] < (biggest_dijkstra * sweep_percentage))
			{
				is_index_front = false;
			}

		}
		if (is_index_front)
		{
			extremums.push_back(gradient_indices[i].second);
		}
	}
	std::sort(gradient_indices.begin(), gradient_indices.end());

	extremums.insert(extremums.end(), points.begin(), points.end());
	if (is_color)
	{
		for (size_t i = 0; i < vertex_size; i++)
		{
			m->raylib_mesh.colors[i * 4] = avg_geodesic_distances[i] / biggest * 255;
			m->raylib_mesh.colors[i * 4 + 1] = 0;
			m->raylib_mesh.colors[i * 4 + 2] = 0;
			m->raylib_mesh.colors[i * 4 + 3] = 255;
		}
		for (size_t i = 0; i < extremums.size(); i++)
		{
			int index = extremums[i];
			m->raylib_mesh.colors[index * 4] = 255;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
		m->update_raylib_mesh();
	}

	return extremums;
}
std::vector<unsigned int> Geodesic_avg_dijkstra(TrilateralMesh* m, int& c, float sweep_percentage, int N_ring,  bool is_color)
{
	std::vector<std::pair<float , unsigned int>> gradient_indices;
	int vertex_size = m->vertices.size();
	std::vector<float> avg_geodesic_distances( vertex_size);
	float biggest = -INFINITY;
	float biggest_dijkstra = -INFINITY; 
	for (size_t i = 0; i < vertex_size; i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, i);

		//sum the distances
		float sum = 0;
		for (size_t j = 0; j < vertex_size; j++)
		{
			sum += distances_i[j];
			if (biggest_dijkstra < distances_i[j])
			{
				biggest_dijkstra = distances_i[j];
			}
		}
		
		sum = sum * m->areas[i] / 3.0;
		avg_geodesic_distances[i] = sum;
		if (sum > biggest)
		{
			biggest = sum;
		}

	}
	
	//average the one ring
	std::vector<float> temp_avg_values = avg_geodesic_distances;
	for (size_t i = 0; i < vertex_size; i++)
	{
		float avg = temp_avg_values[i]; 
		std::vector<unsigned int> neighbours = findNRingNeighbors(m, i, N_ring);
		for (size_t j = 0; j < neighbours.size(); j++)
		{
			int neighbour_index = neighbours[j];
	
			avg += temp_avg_values[neighbour_index];
		}
		avg = avg / (neighbours.size() + 1);
		avg_geodesic_distances[i] = avg; 
	}
	int n_ring = 3;
	//sample where gradient is 0
	for (size_t i = 0; i < vertex_size; i++)
	{
		bool is_maxima = true; 
		bool is_minima = true; 
		
		std::vector<unsigned int> neighbours;
		neighbours = findNRingNeighbors(m, i, N_ring);
		for (size_t j = 0; j < neighbours.size(); j++)
		{
			int index = neighbours[j];
			if (index == i)
			{
				continue;
			}
			if (avg_geodesic_distances[index] > avg_geodesic_distances[i])
			{
				is_maxima = false;
			}
			if (avg_geodesic_distances[index] < avg_geodesic_distances[i])
			{
				is_minima = false;
			}

		}
		/*int neighbour_size = m->adjacenies[i].size();
		for (size_t j = 0; j < neighbour_size; j++)
		{
			int neighbour_index = m->adjacenies[i][j].first;
			if (avg_geodesic_distances[neighbour_index] > avg_geodesic_distances[i])
			{
				is_maxima = false; 
			}
			if (avg_geodesic_distances[neighbour_index] < avg_geodesic_distances[i])
			{
				is_minima = false; 
			}
		}*/
		if (is_maxima || is_minima )
		{
			gradient_indices.push_back(std::make_pair(avg_geodesic_distances[i], i));
		}
	}

	std::sort(gradient_indices.begin(), gradient_indices.end());
	


	std::vector< unsigned int> extremums;
	for (size_t i = 0; i < gradient_indices.size(); i++)
	{
		extremums.push_back(gradient_indices[i].second);
	}
	

	if (is_color)
	{
		for (size_t i = 0; i < vertex_size; i++)
		{
			m->raylib_mesh.colors[i * 4] = avg_geodesic_distances[i] / biggest * 255;
			m->raylib_mesh.colors[i * 4 + 1] = 0;
			m->raylib_mesh.colors[i * 4 + 2] = 0;
			m->raylib_mesh.colors[i * 4 + 3] = 255;
		}
		for (size_t i = 0; i < extremums.size(); i++)
		{
			int index = extremums[i];
			m->raylib_mesh.colors[index * 4] = 255;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
		m->update_raylib_mesh();
	}
	
	return extremums;
	
}
std::vector<unsigned int> Geodesic_min_dijkstra(TrilateralMesh*m , 
std::vector< unsigned int> agd_extremums, float sweep_percentage, float tau, bool is_color)
{
	int vertex_size = m->vertices.size();
	std::vector<std::pair<float, unsigned int>> min_geo_distances(vertex_size); 
	std::vector<std::pair<float, unsigned int>> gradient_indices;
	float maximum_length = -INFINITY; 
	float maximum_min_geodesic = -INFINITY; 
	for (size_t i = 0; i < vertex_size; i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, i);
		for (size_t j = 0; j < distances_i.size(); j++)
		{
			if (maximum_length < distances_i[j])
			{
				maximum_length = distances_i[j];
			}
		}
		float smallest = INFINITY; 
		for (size_t j = 0; j < agd_extremums.size(); j++)
		{
			int ext_index = agd_extremums[j];
			if (distances_i[ext_index] < smallest)
			{
				smallest = distances_i[ext_index];
			}
		}
		min_geo_distances[i] = std::make_pair(smallest , i); 
	}

	//for tau
	for (size_t i = 0; i < min_geo_distances.size(); i++)
	{
		if (maximum_min_geodesic < min_geo_distances[i].first)
		{
			maximum_min_geodesic = min_geo_distances[i].first;
		}
	}
	//sample where gradient is 0
	for (size_t i = 0; i < vertex_size; i++)
	{
		bool is_maxima = true;
		int neighbour_size = m->adjacenies[i].size();
		for (size_t j = 0; j < neighbour_size; j++)
		{
			int neighbour_index = m->adjacenies[i][j].first;
			if (min_geo_distances[neighbour_index] > min_geo_distances[i])
			{
				is_maxima = false;
			}
		}
		if (is_maxima && min_geo_distances[i].first > (tau * maximum_min_geodesic) )
		{
			gradient_indices.push_back(min_geo_distances[i]);
			std::cout << " min geo percentage " << min_geo_distances[i].first<< "  " << maximum_min_geodesic << " "
			<< min_geo_distances[i].first / maximum_min_geodesic << std::endl;
		}
	}
	std::sort(gradient_indices.begin(), gradient_indices.end());

	std::vector< unsigned int> extremums;
	for (size_t i = 0; i < gradient_indices.size() ; i++)
	{
		int index_front = gradient_indices[i].second;
		bool is_index_front = true;
		for (size_t j = 0; j < extremums.size(); j++)
		{
			int ext_index = extremums[j];
			std::vector<float> distances_ext = Geodesic_dijkstra(*m, ext_index);
			if (distances_ext[index_front] < (maximum_length * sweep_percentage))
			{
				is_index_front = false;
			}
		}
		if (is_index_front)
		{
			extremums.push_back(gradient_indices[i].second);
		}
	}

	if (is_color)
	{
		for (size_t i = 0; i < extremums.size(); i++)
		{
			int index = extremums[i];
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 0;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
		m->update_raylib_mesh();
	}
	
	extremums.insert(extremums.end(), agd_extremums.begin(), agd_extremums.end());
	return extremums;
}


void Geodesic_write_sampled_points(TrilateralMesh* m, std::vector<unsigned int>& agd_points)
{
	std::ofstream file("sampled_points.txt");
	for (size_t i = 0; i < agd_points.size(); i++)
	{
		file << " " << agd_points[i] << " " << "\n";
	}
	file.close();

}
void Geodesic_read_sampled_points(TrilateralMesh* m, std::vector<unsigned int>& sampled_points)
{
	std::ifstream file("sampled_points.txt");
	// Read the file line by line
	std::string line;
	while (std::getline(file, line)) {
		std::stringstream ss(line); // Use stringstream to parse the line
		int number;
		std::vector<int> nums;
		// Extract numbers from the line
		while (ss >> number) {
			nums.push_back(number);
		}
		sampled_points.push_back(nums[0]);
	}
}


std::vector<unsigned int> Geodesic_find_biggest_AGD(TrilateralMesh* m, float sweep_percentage , float stop_param  )
{
	std::vector<unsigned int> biggest_agd_indices; 
	int vertex_size = m->vertices.size();
	std::vector<float> avg_geodesic_distances(vertex_size);
	float biggest = -INFINITY;
	float biggest_dijkstra = -INFINITY;
	for (size_t i = 0; i < vertex_size; i++)
	{
		std::vector<float> distances_i = Geodesic_dijkstra(*m, i);
		//sum the distances
		float sum = 0;
		for (size_t j = 0; j < vertex_size; j++)
		{
			sum += distances_i[j];
			if (biggest_dijkstra < distances_i[j])
			{
				biggest_dijkstra = distances_i[j];
			}
		}
		sum = sum * m->areas[i] / 3.0;
		avg_geodesic_distances[i] = sum;
		if (sum > biggest)
		{
			biggest = sum;
		}
	}

	//while( true )
	for (size_t i = 0; i < stop_param; i++)
	{
		auto max_auto= std::max_element(avg_geodesic_distances.begin(), avg_geodesic_distances.end());
		int max_index = std::distance(avg_geodesic_distances.begin() , max_auto);
		float max_val = *max_auto;
		/*if (max_val < (biggest * stop_param))
		{
			break; 
		}*/
		std::vector<float> distances = Geodesic_dijkstra(*m , max_index);
		bool is_new_point_added = true; 
		for (size_t j = 0; j < biggest_agd_indices.size(); j++)
		{
			if (distances[biggest_agd_indices[j]] <  biggest_dijkstra *  sweep_percentage )
			{
				is_new_point_added = false;
				break; 
			}
		}
		if (is_new_point_added)
		{
			biggest_agd_indices.push_back(max_index);
		}
		avg_geodesic_distances.erase(avg_geodesic_distances.begin() + max_index);
	}


	return biggest_agd_indices;


}

unsigned int Geodesic_find_midpoint(TrilateralMesh* m, unsigned int index1, unsigned int index2)
{
	std::vector<int> paths =Geodesic_between_two_points(*m, index1, index2);
	float dist = 0;
	for (size_t i = 0; i < paths.size()-1; i++)
	{
		dist = dist + glm::distance(m->vertices[paths[i]], m->vertices[paths[i + 1]]);
	}
	float half_dist = 0; 
	for (size_t i = 0; i < paths.size()-1; i++)
	{
		half_dist = half_dist + glm::distance(m->vertices[paths[i]], m->vertices[paths[i + 1]]);
		if (half_dist >= dist / 2)
		{
			return paths[i+1];
		}
	}
}

void Geodesic_mid_point_w_AGD(TrilateralMesh* m, unsigned int& p1, unsigned int& p2, float& biggest_dijkstra)
{
	//find agd for all
	int N = m->vertices.size();
	float smallest_val = INFINITY;
	unsigned int smallest_index = -1;
	std::vector<float> minimum_of_maximums;
	for (size_t i = 0; i < N; i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, i);
		float sum = 0;
		float biggest = -INFINITY; 
		for (size_t j = 0; j < N; j++)
		{
			sum += distances[j];
			if (distances[j] > biggest)
			{
				biggest = distances[j];
			}
		}
		if (sum < smallest_val)
		{
			smallest_val = sum;
			smallest_index = i;
		}

		minimum_of_maximums.push_back(biggest);
	} 
	auto smallest_auto = std::min_element(minimum_of_maximums.begin(), minimum_of_maximums.end());
	smallest_index = std::distance(minimum_of_maximums.begin(), smallest_auto);
	/*glm::vec3 p(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		p = p + m->vertices[i];
	}
	p /= m->vertices.size();
	float closest = INFINITY; 
	int closest_index = -1;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		float dist = glm::distance(m->vertices[i], p);
		if (dist < closest)
		{
			closest = dist; 
			closest_index = i; 
		}
	}
	smallest_index = closest_index;*/
	auto best_biggest = std::max_element(minimum_of_maximums.begin(), minimum_of_maximums.end());
	biggest_dijkstra = *best_biggest;
	// that is the smallest index 
	//find the counterpart by using a ray
	p2 = Geodesic_send_ray_get_counterpart(m, smallest_index);
	p1 = smallest_index;

	std::vector<unsigned int> mid_points = {p1,p2};
	m->color_points(mid_points, BLACK);
}
std::vector<unsigned int> Geodesic_generate_secondary_curve_w_midpoints(TrilateralMesh* m, unsigned int& midpoint1 , unsigned int& midpoint2 )
{
	float biggest; 
	Geodesic_mid_point_w_AGD(m, midpoint1, midpoint2, biggest);
	std::vector<unsigned int> points = midpoint_sampling(m, 0.02,biggest,midpoint1 , midpoint2);
	
	std::vector<float> distances_from_midpoint = Geodesic_dijkstra(*m, midpoint1);

	std::vector<std::pair<float,int> > points_curvature; 
	for (size_t i = 0; i < points.size(); i++)
	{
		float dist_from_mid = distances_from_midpoint[i];
		if (dist_from_mid / biggest < 0.4)
		{
			points_curvature.push_back(std::make_pair( dist_from_mid, points[i] ) );
		}
	}
	std::vector<unsigned int> paths; 
	std::sort(points_curvature.begin(), points_curvature.end());
	for (size_t i = 0; i < points_curvature.size()-1; i++)
	{
		std::vector<int> path = Geodesic_between_two_points(*m, points_curvature[i].second, points_curvature[i + 1].second);
		paths.insert(paths.end(), path.begin(), path.end());
	}
	//m->color_points(paths, YELLOW);
	return paths;
}

bool Geodesic_path_intersection(TrilateralMesh* m, std::vector<unsigned int>& path1, std::vector<unsigned int>& path2, unsigned int& no_of_times)
{
	bool is_hit = false; 
	no_of_times = 0; 
	for (size_t i = 0; i < path1.size(); i++)
	{
		int index_i1 = path1[i];
		for (size_t j = 0; j < path2.size(); j++)
		{
			int index_j1 = path2[j];
			if (index_i1 == index_j1)
			{
				is_hit = true; 
				no_of_times++; 
			}
		}
	}
	return is_hit; 
}

std::vector<unsigned int> conv_int_to_unsigned(std::vector<int> vec)
{
	std::vector<unsigned int> u_vec;
	for (size_t i = 0; i < vec.size(); i++)
	{
		u_vec.push_back(vec[i]);
	}
	return u_vec;
}

// give the first point and fect the second 
std::vector<unsigned int> Geodesic_generate_secondary_curve(TrilateralMesh* m, unsigned int& midpoint1, unsigned int& midpoint2)
{
	midpoint2 = Geodesic_send_ray_get_counterpart(m, midpoint1);
	std::vector<unsigned int> path = conv_int_to_unsigned(Geodesic_between_two_points(*m, midpoint1, midpoint2));
	unsigned int midpoint_of_midpoints;
	float total_length = 0;
	for (size_t i = 0; i < path.size() - 1; i++)
	{
		total_length += glm::distance(m->vertices[path[i]], m->vertices[path[i + 1]]);
	}
	float half_length = 0;
	for (size_t i = 0; i < path.size() - 1; i++)
	{
		half_length += glm::distance(m->vertices[path[i]], m->vertices[path[i + 1]]);
		if (half_length >= total_length / 2)
		{
			midpoint_of_midpoints = path[i];
			break;
		}
	}
	unsigned int midpoint_of_midpoint_2 = Geodesic_send_ray_get_counterpart(m, midpoint_of_midpoints);
	std::vector<unsigned int> path_m1_to_m = conv_int_to_unsigned(Geodesic_between_two_points(*m, midpoint1, midpoint_of_midpoint_2));
	std::vector<unsigned int> path_m_to_m2 = conv_int_to_unsigned(Geodesic_between_two_points(*m, midpoint_of_midpoint_2, midpoint2));

	for (size_t i = path_m1_to_m.size() - 1; i > 0; i--)
	{
		path.push_back(path_m1_to_m[i]);
	}
	for (size_t i = path_m_to_m2.size() - 1; i > 0; i--)
	{
		path.push_back(path_m_to_m2[i]);
	}
	m->color_points(path, YELLOW);

	std::vector<unsigned int> midpoints = { midpoint1 };
	m->color_points(midpoints, RED);
	midpoints.clear();
	midpoints.push_back(midpoint2);
	m->color_points(midpoints, GREEN);
	midpoints.clear();
	midpoints.push_back(midpoint_of_midpoint_2);
	m->color_points(midpoints, BLUE);

	return path;
}
std::vector<unsigned int> Geodesic_generate_multiple_secondary_curve(TrilateralMesh* m, unsigned int& midpoint1, unsigned int& midpoint2)
{

	std::vector<unsigned int> mid_curve = Geodesic_generate_secondary_curve_w_midpoints(m, midpoint1, midpoint2);

	glm::vec3 dif = m->vertices[mid_curve[1]] - m->vertices[mid_curve[0]];
	dif = glm::normalize(dif);

	glm::vec3 new_dir = glm::normalize( glm::cross(m->normals[midpoint2], dif));
	glm::vec3 mid_point1 = m->vertices[midpoint2];

	// lets do bfs 
	std::stack<int> stack;  // a stack consisting of indices
	// get the adjacencies
	std::vector<std::vector<std::pair<int, float>>> mesh_adjacencies = m->adjacenies;
	//lastly get a int array with  size of vertices in order to check if the vertex has been visited ( -1 edge , 0 not visisted , 1 visited) 
	std::vector<unsigned int> visited_points;
	std::vector<int> is_visited(m->vertices.size(), 0);
	for (size_t i = 0; i < mid_curve.size(); i++)
	{
		is_visited[mid_curve[i]] = -1;
	}
	//push our point to stack
	stack.push(0); // random index 0 chosen 
	while (!stack.empty())
	{
		int index = stack.top();
		stack.pop(); //vertex index popped from stack
		visited_points.push_back(index);
		if (is_visited[index] == 0) //not visited
		{
			is_visited[index] = 1; // now te vertex has been visited

			// this region of loop assumes index is not edge, therefore add the adjacencies
			for (size_t i = 0; i < mesh_adjacencies[index].size(); i++) //process pairs 
			{
				stack.push(mesh_adjacencies[index][i].first);
			}
		}
		if (is_visited[index] == -1) //do nothing 
		{
			;
		}
	}

	//now get a new average on stack
	glm::vec3 avg(0, 0, 0);
	for (size_t i = 0; i < visited_points.size(); i++)
	{
		avg = avg + m->vertices[visited_points[i]];
	}
	
	avg = avg / (float)visited_points.size();
	unsigned int closest_index = mesh_get_closest_index(m, avg);

	std::vector<unsigned int> vertices = { closest_index };
	m->color_points(visited_points, RED);
	m->color_points(vertices, BLACK);

	return mid_curve;
}

unsigned int Geodesic_send_ray_get_counterpart(TrilateralMesh* m, unsigned int& midpoint1)
{
	TrilateralRay ray;
	ray.origin = m->vertices[midpoint1];
	glm::vec3 vec1(m->normals_display[midpoint1 * 12], m->normals_display[midpoint1 * 12 + 1], m->normals_display[midpoint1 * 12 + 2]);
	glm::vec3 vec2(m->normals_display[midpoint1 * 12 + 6], m->normals_display[midpoint1 * 12 + 7], m->normals_display[midpoint1 * 12 + 8]);
	glm::vec3 dir = vec2 - vec1;
	dir = dir * -1.0f;
	dir = glm::normalize(dir);
	ray.direction = dir;
	float smallest_dist = INFINITY;
	unsigned int smallest_hit_index = -1;
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		glm::vec3 hit_point;
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		bool is_hit = ray_triangle_intersection(ray, m->vertices[index1], m->vertices[index2], m->vertices[index3], hit_point);
		bool is_same_point = (index1 == midpoint1) || (index2 == midpoint1) || (index3 == midpoint1);
		if (is_hit && !is_same_point)
		{
			float dist = glm::distance(m->vertices[midpoint1], hit_point);
			if (smallest_dist > dist)
			{
				dist = smallest_dist;
				//get the closest dist on triangle 
				float disti = glm::distance(hit_point, m->vertices[index1]);
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
	std::cout << " ray cast point " << smallest_hit_index<<  std::endl;
	return smallest_hit_index; 
}

void Geodesic_color_path(TrilateralMesh * m , unsigned int p1 , unsigned int p2 )
{
	std::vector<int> path = Geodesic_between_two_points(*m, p1, p2);
	m->color_points(conv_int_to_unsigned(path), YELLOW);
}

void Geodesic_color_according_to_midpoints(TrilateralMesh* m)
{
	unsigned int mid1, mid2;
	float dijkstra_biggest;
	Geodesic_mid_point_w_AGD(m, mid1, mid2, dijkstra_biggest);
	std::vector<unsigned int> closer_to_1;
	std::vector<unsigned int> closer_to_2;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		float dot1 = glm::dot(m->normals[i], m->normals[mid1]);
		float dot2 = glm::dot(m->normals[i], m->normals[mid2]);
		if (dot1 > dot2)
		{
			closer_to_1.push_back(i);
		}
		else
		{
			closer_to_2.push_back(i);
		}
	}

	m->color_points(closer_to_1, RED);
	m->color_points(closer_to_2, BLUE);

	std::vector<unsigned int> mid = { mid1,  mid2 };
	m->color_points(mid, BLACK);
}
unsigned int Geodesic_get_midpoint_from_path(TrilateralMesh* m, unsigned int p1, unsigned int p2)
{
	// get mid point of newly added 
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
	for (size_t j = 0; j < point_list.size()-1; j++)
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