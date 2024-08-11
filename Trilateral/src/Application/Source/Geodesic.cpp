#include "../Include/Geodesic.h"
std::vector<float> Geodesic_dijkstra(Mesh& m, int point_index)
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
std::vector<int> Geodesic_dijkstra_predecessors(Mesh& m, int point_index)
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
std::vector<int> Geodesic_between_two_points(Mesh& m, int p1_index, int p2_index)
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