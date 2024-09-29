#pragma once
/*
#include "TrilateralMesh.h"
static void buffer_mesh(TrilateralMesh& m)
{
	//first convert to float

	//then buffer 
	//convert to float buffer before loading
	std::vector<float> temp_array;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		temp_array.push_back(m.vertices[i].x);
		temp_array.push_back(m.vertices[i].y);
		temp_array.push_back(m.vertices[i].z);
		//default zero
		temp_array.push_back(0.0f);
		temp_array.push_back(0.0f);
		temp_array.push_back(0.0f);
	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.vertex_indices.size() * sizeof(int), &m.vertex_indices[0], GL_STATIC_DRAW);
	;

}
static void buffer_mesh(TrilateralMesh& m, std::vector<int> path)
{
	//first convert to float

	//then buffer 
	//convert to float buffer before loading
	std::vector<float> temp_array;
	bool is_color_changed = false;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_color_changed = false;
		temp_array.push_back(m.vertices[i].x);
		temp_array.push_back(m.vertices[i].y);
		temp_array.push_back(m.vertices[i].z);
		for (size_t j = 0; j < path.size(); j++)
		{
			if (path[j] == i)
			{
				is_color_changed = true;
				temp_array.push_back(1.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}
		//default zero
		if (!is_color_changed)
		{
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}

	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);
	;

}
static void draw_buffer(TrilateralMesh& m)
{
	glDrawElements(GL_LINES, m.vertex_indices.size(), GL_UNSIGNED_INT, 0);
}
#pragma region assigment 1 
// find the distance and path between all points
static std::vector<int> compute_geodesic_distances_array(TrilateralMesh& m, int point_index)
{
	//init array 
	float* matrix = new float[m.vertices.size()];
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if (point_index == i)
		{
			matrix[i] = 0;
		}
		else
		{
			//INFINITY
			matrix[i] = (float)INFINITY;
		}
	}
	//init extra stuff
	int vertex_visited = 0; //total umber of visits 
	unsigned int* if_vertex_visited = new unsigned int[m.vertices.size()];
	int* predecessor = new int[m.vertices.size()]; //predecessor vertices for all 
	//float* will_be_visited_list = new float[m.vertices.size()];
	// std::vector<std::pair<int , float>> shortest_path_set; //not yet included but calculated 
	//fill array 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if_vertex_visited[i] = 0; // not visited
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessor[i] = -1; // not visited
	}
	// dijkstra
	if_vertex_visited[point_index] = 1;
	vertex_visited = vertex_visited + 1;
	int current_vertex_no = point_index;
	predecessor[point_index] = point_index;
	while (vertex_visited < m.vertices.size())
	{
		//1 - look around of the current index and update them
		for (size_t i = 0; i < m.adjacenies[current_vertex_no].size(); i++)
		{
			int adjacent_vertex = m.adjacenies[current_vertex_no][i].first;
			if (if_vertex_visited[adjacent_vertex] == 0) //if not visited 
			{
				if (matrix[current_vertex_no] + m.adjacenies[current_vertex_no][i].second < matrix[adjacent_vertex])
				{
					matrix[adjacent_vertex] = matrix[current_vertex_no] + m.adjacenies[current_vertex_no][i].second;
					predecessor[adjacent_vertex] = current_vertex_no;
				}
			}
		}
		// 2 - after update select the minimum non visited vertex from matrix and if_visited array
		int minimum_size = (float)INFINITY;
		for (size_t i = 0; i < m.vertices.size(); i++)
		{
			if (if_vertex_visited[i] == 0 && matrix[i] != (float)INFINITY) //update 
			{
				if (matrix[i] < minimum_size)
				{
					minimum_size = matrix[i];
					current_vertex_no = i;
				}

			}
		}
		vertex_visited += 1;
		if_vertex_visited[current_vertex_no] = 1;

	}

	std::vector<int> predecessors;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessors.push_back(predecessor[i]);
	}
	//destroy arrays
	delete matrix;
	//delete extra stuff 
	delete if_vertex_visited;
	delete predecessor;

	return predecessors;
}
static std::vector<float> compute_geodesic_distances_array_distance(TrilateralMesh& m, int point_index)
{
	//init array 
	float* matrix = new float[m.vertices.size()];
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if (point_index == i)
		{
			matrix[i] = 0;
		}
		else
		{
			//INFINITY
			matrix[i] = (float)INFINITY;
		}
	}
	//init extra stuff
	int vertex_visited = 0; //total umber of visits 
	unsigned int* if_vertex_visited = new unsigned int[m.vertices.size()];
	int* predecessor = new int[m.vertices.size()]; //predecessor vertices for all 
	//float* will_be_visited_list = new float[m.vertices.size()];
	// std::vector<std::pair<int , float>> shortest_path_set; //not yet included but calculated 
	//fill array 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if_vertex_visited[i] = 0; // not visited
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessor[i] = -1; // not visited
	}
	// dijkstra
	if_vertex_visited[point_index] = 1;
	vertex_visited = vertex_visited + 1;
	int current_vertex_no = point_index;
	predecessor[point_index] = point_index;
	while (vertex_visited < m.vertices.size())
	{
		//1 - look around of the current index and update them
		for (size_t i = 0; i < m.adjacenies[current_vertex_no].size(); i++)
		{
			int adjacent_vertex = m.adjacenies[current_vertex_no][i].first;
			if (if_vertex_visited[adjacent_vertex] == 0) //if not visited 
			{
				if (matrix[current_vertex_no] + m.adjacenies[current_vertex_no][i].second < matrix[adjacent_vertex])
				{
					matrix[adjacent_vertex] = matrix[current_vertex_no] + m.adjacenies[current_vertex_no][i].second;
					predecessor[adjacent_vertex] = current_vertex_no;
				}
			}
		}
		// 2 - after update select the minimum non visited vertex from matrix and if_visited array
		int minimum_size = (float)INFINITY;
		for (size_t i = 0; i < m.vertices.size(); i++)
		{
			if (if_vertex_visited[i] == 0 && matrix[i] != (float)INFINITY) //update 
			{
				if (matrix[i] < minimum_size)
				{
					minimum_size = matrix[i];
					current_vertex_no = i;
				}

			}
		}
		vertex_visited += 1;
		if_vertex_visited[current_vertex_no] = 1;

	}

	std::vector<int> predecessors;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessors.push_back(predecessor[i]);
	}
	std::vector<float> distances;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		distances.push_back(matrix[i]);
	}
	//destroy arrays
	delete matrix;
	//delete extra stuff 
	delete if_vertex_visited;
	delete predecessor;


	return distances;
}
static std::vector<int> compute_geodesic_distances_min_heap(TrilateralMesh& m, int point_index)
{
	//create min heap
	std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>> > min_heap;
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
	min_heap.push(first_index);
	while (!min_heap.empty())
	{
		current_pair = min_heap.top();
		min_heap.pop();
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
					min_heap.push(new_vertex);
				}
			}
		}
	}
	// copy the predecessor 
	std::vector<int> predecessors;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		predecessors.push_back(predecessor[i]);
	}
	delete matrix;
	delete is_vertex_visited;
	delete predecessor;

	return predecessors;
}
static std::vector<float> compute_geodesic_distances_min_heap_distances(TrilateralMesh& m, int point_index)
{
	//create min heap
	std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>> > min_heap;
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
	min_heap.push(first_index);
	while (!min_heap.empty())
	{
		current_pair = min_heap.top();
		min_heap.pop();
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
					min_heap.push(new_vertex);
				}
			}
		}
	}
	// copy the predecessor 
	std::vector<float> distances;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		distances.push_back(matrix[i]);
	}
	delete matrix;
	delete is_vertex_visited;
	delete predecessor;
	return distances;
}
static std::vector<int> Geodesic_dijkstra_predecessors(TrilateralMesh& m, int point_index)
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
static std::vector<float> Geodesic_dijkstra(TrilateralMesh& m, int point_index)
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
static std::vector<int> draw_with_heap_implementation(TrilateralMesh& m, int p1_index, int p2_index)
{
	std::vector<int> predecessors = compute_geodesic_distances_min_heap(m, p1_index);
	std::vector<int> consec_indices;
	int pred = p2_index;
	while (pred != p1_index)
	{
		consec_indices.push_back(pred);
		pred = predecessors[pred];
	}
	return consec_indices;
}
static std::vector<int> draw_with_array_implementation(TrilateralMesh& m, int p1_index, int p2_index)
{
	std::vector<int> predecessors = compute_geodesic_distances_array(m, p1_index);
	std::vector<int> consec_indices;
	int pred = p2_index;
	while (pred != p1_index)
	{
		consec_indices.push_back(pred);
		pred = predecessors[pred];
	}
	return consec_indices;
}

static void compute_all(TrilateralMesh& m)
{

	// part 1 - compute for arrays and output
	std::ofstream file;
	file.open("matrixM.txt", std::ofstream::out | std::ofstream::app);
	std::vector<float> all_output_float;
	auto start_time = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_array_distance(m, i);
		for (size_t j = 0; j < m.vertices.size(); j++)
		{
			all_output_float.push_back(distances[j]);
		}
	}
	auto end_time = std::chrono::high_resolution_clock::now();

	auto passed_time = end_time - start_time;
	std::cout << "time passed with array " << passed_time / std::chrono::milliseconds(1) << std::endl;
	for (size_t i = 0; i < all_output_float.size(); i++)
	{
		file << all_output_float[i] << "  ";
	}
	file.close();

	// part 2 - compute for min heap
	file.open("minheapM.txt", std::ofstream::out | std::ofstream::app);
	all_output_float.clear();
	start_time = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_min_heap_distances(m, i);
		for (size_t j = 0; j < m.vertices.size(); j++)
		{
			all_output_float.push_back(distances[j]);
		}
	}
	end_time = std::chrono::high_resolution_clock::now();
	passed_time = end_time - start_time;
	std::cout << "time passed with min heap " << passed_time / std::chrono::milliseconds(1) << std::endl;
	for (size_t i = 0; i < all_output_float.size(); i++)
	{
		file << all_output_float[i] << "  ";
	}
	file.close();


	// part 3 - compute for fibonacci heap
	file.open("fibheapM.txt", std::ofstream::out | std::ofstream::app);
	all_output_float.clear();
	start_time = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_min_heap_distances(m, i);
		for (size_t j = 0; j < m.vertices.size(); j++)
		{
			all_output_float.push_back(distances[j]);
		}
	}
	end_time = std::chrono::high_resolution_clock::now();
	passed_time = end_time - start_time;
	std::cout << "time passed with fib heap " << passed_time / std::chrono::milliseconds(1) << std::endl;
	for (size_t i = 0; i < all_output_float.size(); i++)
	{
		file << all_output_float[i] << "  ";
	}
	file.close();
}
#pragma endregion assignemnt 1 

#pragma region assingment part 2

// this version returns the matrix 
//eventually all of them must be 

static std::vector<int> compute_iso_curves(TrilateralMesh& m, int point_index_1)
{
	std::vector<int> predecessors = Geodesic_dijkstra_predecessors(m, point_index_1);
	std::vector<int> distances; // distances in terms of edge no not edge length 
	//init 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		distances.push_back(0);
	}
	predecessors[point_index_1] = point_index_1;
	//calculate geodesic distance  
	for (size_t i = 0; i < predecessors.size(); i++)
	{
		int distance = 1;
		int vertex = predecessors[i];
		while (vertex != point_index_1)
		{
			distance++;
			vertex = predecessors[vertex];
		}
		distances[i] = distance;
	}
	distances[point_index_1] = 0;
	// now we will draw the isocurves
	std::vector<float> temp_array;
	for (size_t i = 0; i < m.triangles.size(); i = i + 3)
	{
		int vertex_index_1 = m.triangles[i];
		int vertex_index_2 = m.triangles[i + 1];
		int vertex_index_3 = m.triangles[i + 2];

		// check the distance steps of each of the vertex distances
		int distance_step_1 = distances[vertex_index_1];
		int distance_step_2 = distances[vertex_index_2];
		int distance_step_3 = distances[vertex_index_3];

		// push the points accordingly 
		// if two vertex have same distance steps color it to blue

		// v1 and v2 
		temp_array.push_back(m.vertices[vertex_index_1].x);
		temp_array.push_back(m.vertices[vertex_index_1].y);
		temp_array.push_back(m.vertices[vertex_index_1].z);
		if (vertex_index_1 == point_index_1)
		{
			//red 
			temp_array.push_back(1.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		else
		{
			if (distance_step_1 == distance_step_2)
			{
				// if equal push blue
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(1.0f);
			}
			else
			{
				// if not equal push black
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}

		temp_array.push_back(m.vertices[vertex_index_2].x);
		temp_array.push_back(m.vertices[vertex_index_2].y);
		temp_array.push_back(m.vertices[vertex_index_2].z);
		if (vertex_index_2 == point_index_1)
		{
			//red 
			temp_array.push_back(1.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		else
		{
			if (distance_step_1 == distance_step_2)
			{
				// if equal push blue
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(1.0f);
			}
			else
			{
				// if not equal push black
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}


		// vertex v1 and v3 
		temp_array.push_back(m.vertices[vertex_index_1].x);
		temp_array.push_back(m.vertices[vertex_index_1].y);
		temp_array.push_back(m.vertices[vertex_index_1].z);
		if (vertex_index_1 == point_index_1)
		{
			//red 
			temp_array.push_back(1.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		else
		{
			if (distance_step_1 == distance_step_3)
			{
				// if equal push blue
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(1.0f);
			}
			else
			{
				// if not equal push black
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}
		temp_array.push_back(m.vertices[vertex_index_3].x);
		temp_array.push_back(m.vertices[vertex_index_3].y);
		temp_array.push_back(m.vertices[vertex_index_3].z);
		if (vertex_index_3 == point_index_1)
		{
			//red 
			temp_array.push_back(1.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		else
		{
			if (distance_step_1 == distance_step_3)
			{
				// if equal push blue
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(1.0f);
			}
			else
			{
				// if not equal push black
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}



		// vertex v2 and v3 
		temp_array.push_back(m.vertices[vertex_index_2].x);
		temp_array.push_back(m.vertices[vertex_index_2].y);
		temp_array.push_back(m.vertices[vertex_index_2].z);
		if (vertex_index_2 == point_index_1)
		{
			//red 
			temp_array.push_back(1.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		else
		{
			if (distance_step_2 == distance_step_3)
			{
				// if equal push blue
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(1.0f);
			}
			else
			{
				// if not equal push black
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}

		temp_array.push_back(m.vertices[vertex_index_3].x);
		temp_array.push_back(m.vertices[vertex_index_3].y);
		temp_array.push_back(m.vertices[vertex_index_3].z);

		if (vertex_index_3 == point_index_1)
		{
			//red 
			temp_array.push_back(1.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		else
		{
			if (distance_step_2 == distance_step_3)
			{
				// if equal push blue
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(1.0f);
			}
			else
			{
				// if not equal push black
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}



	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);
	;
	return distances;
}
static std::vector<float> compute_curve_distances(TrilateralMesh& m, std::vector<int> distances, int point_index1)
{
	std::vector<float> geodesic_distances = Geodesic_dijkstra(m, point_index1);
	std::vector<float> total_distances;

	// get the max from distances
	int max = 0;
	for (size_t i = 0; i < distances.size(); i++)
	{
		if (distances[i] > max)
		{
			max = distances[i];
		}
	}
	// init total_distances
	for (size_t i = 0; i < max + 1; i++)
	{
		total_distances.push_back(0);
	}

	// for each r i 
	// search every triangle that has two components with the same distance ri 
	// there must be two of these triangles, therefore the third index should be more far than v1 and v2 
	glm::vec3 origin_point = m.vertices[point_index1];

	for (size_t i = 0; i < m.triangles.size(); i = i + 3)
	{
		int index1 = m.triangles[i];
		int index2 = m.triangles[i + 1];
		int index3 = m.triangles[i + 2];
		if (distances[index1] == distances[index2])
		{
			// calculate the distances
			float dist_point_known = glm::distance(m.vertices[index1], m.vertices[point_index1]);
			float dist_point_unknown = glm::distance(m.vertices[index3], m.vertices[point_index1]);
			if (dist_point_known > dist_point_unknown)
			{
				continue;
			}
			// from now on we have found the correct triangle
			// index 3 is v0 
			int v0 = index3;
			int v1 = index2;
			int v2 = index1;
			int radius = v1;
			float g0 = geodesic_distances[v0];
			float g1 = geodesic_distances[v1];
			float g2 = geodesic_distances[v2];
			float alpha_1 = glm::abs(radius - g0) / glm::abs(g1 - g0);
			float alpha_2 = glm::abs(radius - g0) / glm::abs(g1 - g0);
			glm::vec3 p1 = (1 - alpha_1) * m.vertices[v0] + alpha_1 * m.vertices[v1];
			glm::vec3 p2 = (1 - alpha_2) * m.vertices[v0] + alpha_2 * m.vertices[v2];

			total_distances[distances[v1]] += glm::distance(p1, p2);

		}
		else if (distances[index1] == distances[index3])
		{
			float dist_point_known = glm::distance(m.vertices[index1], m.vertices[point_index1]);
			float dist_point_unknown = glm::distance(m.vertices[index2], m.vertices[point_index1]);
			if (dist_point_known > dist_point_unknown)
			{
				continue;
			}
			// from now on we have found the correct triangle
			// index 2 is v0 
			int v0 = index2;
			int v1 = index3;
			int v2 = index1;
			int radius = v1;
			float g0 = geodesic_distances[v0];
			float g1 = geodesic_distances[v1];
			float g2 = geodesic_distances[v2];
			float alpha_1 = glm::abs(radius - g0) / glm::abs(g1 - g0);
			float alpha_2 = glm::abs(radius - g0) / glm::abs(g1 - g0);
			glm::vec3 p1 = (1 - alpha_1) * m.vertices[v0] + alpha_1 * m.vertices[v1];
			glm::vec3 p2 = (1 - alpha_2) * m.vertices[v0] + alpha_2 * m.vertices[v2];

			total_distances[distances[v1]] += glm::distance(p1, p2);
		}
		else if (distances[index2] == distances[index3])
		{
			float dist_point_known = glm::distance(m.vertices[index2], m.vertices[point_index1]);
			float dist_point_unknown = glm::distance(m.vertices[index1], m.vertices[point_index1]);
			if (dist_point_known > dist_point_unknown)
			{
				continue;
			}
			// from now on we have found the correct triangle
			// index 2 is v0 
			int v0 = index1;
			int v1 = index3;
			int v2 = index2;
			int radius = v1;
			float g0 = geodesic_distances[v0];
			float g1 = geodesic_distances[v1];
			float g2 = geodesic_distances[v2];
			float alpha_1 = glm::abs(radius - g0) / glm::abs(g1 - g0);
			float alpha_2 = glm::abs(radius - g0) / glm::abs(g1 - g0);
			glm::vec3 p1 = (1 - alpha_1) * m.vertices[v0] + alpha_1 * m.vertices[v1];
			glm::vec3 p2 = (1 - alpha_2) * m.vertices[v0] + alpha_2 * m.vertices[v2];

			total_distances[distances[v1]] += glm::distance(p1, p2);
		}

	}
	float min_distance = INFINITY;
	// we will have 2 of those triangles , divide all of the vector data to half
	for (size_t i = 0; i < total_distances.size(); i++)
	{
		total_distances[i] /= 2;
		if (min_distance > total_distances[i] && total_distances[i] != 0.0f)
		{
			min_distance = total_distances[i];
		}

	}
	//normalize distances
	for (size_t i = 0; i < total_distances.size(); i++)
	{
		total_distances[i] /= min_distance;
		std::cout << (int)total_distances[i] << std::endl;
	}
	//returns the normalized distances
	return total_distances;
}
static int find_point_from_histogram(TrilateralMesh& m, std::vector<float> histogram)
{
	float min_distance = INFINITY;
	int best_index;

	std::vector<int> point_list;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		bool is_equal = true;
		std::vector<int> distances = compute_iso_curves(m, i);
		std::vector<float> new_histogram = compute_curve_distances(m, distances, i);
		// compare histogram

		// get the minimum of them
		int min_size = glm::min(new_histogram.size(), histogram.size());
		float total_difference = 0.0f;
		for (size_t j = 0; j < min_size; j++)
		{
			total_difference += (new_histogram[j] - histogram[j]);
		}
		if (min_distance > total_difference)
		{
			min_distance = total_difference;
			best_index = i;
		}
	}
	return best_index;

}
/*static void compute_iso_curves(TrilateralMesh &m , int point_index_1 , int k , float offset   ) //k is the number of times
{
	//number of
	std::vector<float> distances = Geodesic_dijkstra(m, point_index_1);
	// 1- get the max of it so we now how to divide them
	float max = 0;
	for (size_t i = 0; i < distances.size(); i++)
	{
		if (max < distances[i])
		{
			max = distances[i];
		}

	}
	// divide by k
	float distance_step = max / (float)k;

	std::vector<float> temp_array;
	for (size_t i = 0; i < m.triangles.size(); i = i + 3 )
	{
		int vertex_index_1 = m.triangles[i];
		int vertex_index_2 = m.triangles[i + 1 ];
		int vertex_index_3 = m.triangles[i + 2 ];

		// check the distance steps of each of the vertex distances
		int distance_step_1 = distances[vertex_index_1] / distance_step;
		int distance_step_2 = distances[vertex_index_2] / distance_step;
		int distance_step_3 = distances[vertex_index_3] / distance_step;

		// push the points accordingly
		// if two vertex have same distance steps color it to blue

		// v1 and v2
		temp_array.push_back(m.vertices[vertex_index_1].x);
		temp_array.push_back(m.vertices[vertex_index_1].y);
		temp_array.push_back(m.vertices[vertex_index_1].z);
		if (distance_step_1 == distance_step_2)
		{
			// if equal push blue
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(1.0f);
		}
		else
		{
			// if not equal push black
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		temp_array.push_back(m.vertices[vertex_index_2].x);
		temp_array.push_back(m.vertices[vertex_index_2].y);
		temp_array.push_back(m.vertices[vertex_index_2].z);
		if (distance_step_1 == distance_step_2)
		{
			// if equal push blue
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(1.0f);
		}
		else
		{
			// if not equal push black
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}

		// vertex v1 and v3
		temp_array.push_back(m.vertices[vertex_index_1].x);
		temp_array.push_back(m.vertices[vertex_index_1].y);
		temp_array.push_back(m.vertices[vertex_index_1].z);
		if (distance_step_1 == distance_step_3)
		{
			// if equal push blue
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(1.0f);
		}
		else
		{
			// if not equal push black
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		temp_array.push_back(m.vertices[vertex_index_3].x);
		temp_array.push_back(m.vertices[vertex_index_3].y);
		temp_array.push_back(m.vertices[vertex_index_3].z);
		if (distance_step_1 == distance_step_3)
		{
			// if equal push blue
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(1.0f);
		}
		else
		{
			// if not equal push black
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}


		// vertex v2 and v3
		temp_array.push_back(m.vertices[vertex_index_2].x);
		temp_array.push_back(m.vertices[vertex_index_2].y);
		temp_array.push_back(m.vertices[vertex_index_2].z);
		if (distance_step_2 == distance_step_3)
		{
			// if equal push blue
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(1.0f);
		}
		else
		{
			// if not equal push black
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}
		temp_array.push_back(m.vertices[vertex_index_3].x);
		temp_array.push_back(m.vertices[vertex_index_3].y);
		temp_array.push_back(m.vertices[vertex_index_3].z);
		if (distance_step_2 == distance_step_3)
		{
			// if equal push blue
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(1.0f);
		}
		else
		{
			// if not equal push black
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}

	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);

	//dispose matrix
}*/
/*void draw_edges_for_isocurve(TrilateralMesh &m )
{
	glDrawArrays(GL_LINES, 0, m.triangles.size() * 2);
} */
#pragma endregion assignment part 2
/*
#pragma region assignemnt part3 
static void buffer_vertices_with_triangles(TrilateralMesh& m)
{
	std::vector<float> temp_array;
	bool is_color_changed = false;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_color_changed = false;
		temp_array.push_back(m.vertices[i].x);
		temp_array.push_back(m.vertices[i].y);
		temp_array.push_back(m.vertices[i].z);
		temp_array.push_back(0.0f);
		temp_array.push_back(0.0f);
		temp_array.push_back(0.0f);


	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);
}
static void buffer_vertices_with_triangles(TrilateralMesh& m, std::vector<int> path)
{
	std::vector<float> temp_array;
	bool is_color_changed = false;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_color_changed = false;
		temp_array.push_back(m.vertices[i].x);
		temp_array.push_back(m.vertices[i].y);
		temp_array.push_back(m.vertices[i].z);
		for (size_t j = 0; j < path.size(); j++)
		{
			if (path[j] == i)
			{
				is_color_changed = true;
				temp_array.push_back(1.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}
		//default zero
		if (!is_color_changed)
		{
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}

	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);
}
static void buffer_vertices_with_triangles_region_of_interest(TrilateralMesh& m, std::vector<int> path)
{
	std::vector<float> temp_array;
	bool is_color_changed = false;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_color_changed = false;
		temp_array.push_back(m.vertices[i].x);
		temp_array.push_back(m.vertices[i].y);
		temp_array.push_back(m.vertices[i].z);
		for (size_t j = 0; j < path.size(); j++)
		{
			if (path[j] == i)
			{
				is_color_changed = true;
				temp_array.push_back(1.0f);
				temp_array.push_back(0.0f);
				temp_array.push_back(0.0f);
			}
		}
		//default zero
		if (!is_color_changed)
		{
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
			temp_array.push_back(0.0f);
		}

	}
	glBufferData(GL_ARRAY_BUFFER, temp_array.size() * sizeof(float), &temp_array[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);
}
static void draw_vertices_with_triangles(TrilateralMesh& m)
{
	glDrawElements(GL_TRIANGLES, m.triangles.size(), GL_UNSIGNED_INT, 0);
}
static float compute_triangle_area(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 cross1 = glm::cross(p1 - p2, p1 - p3);
	float length = glm::length(cross1) / 2;
	return length;
}
static void compute_bilateral_map(TrilateralMesh& m, int point_index1, int point_index2, float tau, int division_no) // tau is the closeness division_no is the no of how much you want to separate
{
	//1 - extract the path from point1 to point2 
	std::vector<int> path = Geodesic_between_two_points(m, point_index1, point_index2);
	//2 - we should calculate the geodesic distances for all of the path and than redraw
	//declare variable 
	std::vector<float> histogram;
	//init variable
	for (size_t i = 0; i < division_no; i++)
	{
		histogram.push_back(0);
	}

	int* is_close = new int[m.vertices.size()];
	int* is_in_path = new int[m.vertices.size()];
	float* minimum_distance = new float[m.vertices.size()];
	float* tau_chart = new float[division_no];
	//init variable
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_close[i] = 0; //false for all 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_in_path[i] = 0; //false for all 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		minimum_distance[i] = (float)INFINITY; //false for all 
	}

	float step = 0;
	float histogram_no = tau / division_no;
	for (size_t i = 0; i < division_no; i++)
	{
		tau_chart[i] = step;
		step += tau / division_no;
	}

	for (size_t i = 0; i < path.size(); i++)
	{
		is_in_path[path[i]] = true; //get the path ones 
	}
	//calculate the distance from every vertex within path
	for (size_t i = 0; i < path.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(m, path[i]);
		// if a point is close change the flag 
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (distances[j] <= tau)
			{
				is_close[j] = 1;
				if (distances[j] < minimum_distance[j]) //get the minimum distances to p 
				{
					minimum_distance[j] = distances[j];
				}
			}
		}
	}
	//3 - now rebuffer the data 
	// create a vector and push the correct vertices
	std::vector<int> close_indices;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		if (is_close[i] == 1) // the vertex is closer than tau 
		{
			close_indices.push_back(i);
		}
	}
	//rebuffer here
	std::vector<float> new_buffer;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		//now, if it is in the path continue in red
		// if not path adn tau blue
		// else black
		new_buffer.push_back(m.vertices[i].x);
		new_buffer.push_back(m.vertices[i].y);
		new_buffer.push_back(m.vertices[i].z);
		if (is_close[i] == 1)
		{
			if (is_in_path[i] == 1)
			{
				//red
				new_buffer.push_back(1.0f);
				new_buffer.push_back(0.0f);
				new_buffer.push_back(0.0f);

			}
			else
			{
				int div_no = minimum_distance[i] / division_no;
				if (div_no % 2 == 0)
				{
					//blue 
					new_buffer.push_back(0.0f);
					new_buffer.push_back(0.0f);
					new_buffer.push_back(1.0f);
				}
				else
				{
					//green 
					new_buffer.push_back(0.0f);
					new_buffer.push_back(1.0f);
					new_buffer.push_back(0.0f);
				}

			}
		}
		else
		{
			//black 
			new_buffer.push_back(0.0f);
			new_buffer.push_back(0.0f);
			new_buffer.push_back(0.0f);
		}

	}
	glBufferData(GL_ARRAY_BUFFER, new_buffer.size() * sizeof(float), &new_buffer[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);

	for (size_t i = 0; i < m.triangles.size() / 3; i = i + 3)
	{
		float dist1 = minimum_distance[m.triangles[i]];
		float dist2 = minimum_distance[m.triangles[i + 1]];
		float dist3 = minimum_distance[m.triangles[i + 2]];

		int no_1 = dist1 / histogram_no;
		int no_2 = dist2 / histogram_no;
		int no_3 = dist3 / histogram_no;

		int min = INFINITY;
		if (no_1 < min)
		{
			min = no_1;
		}
		if (no_2 < min)
		{
			min = no_2;
		}
		if (no_3 < min)
		{
			min = no_3;
		}
		//put the traingles into min index of histogram
		float area = compute_triangle_area(m.vertices[m.triangles[i]], m.vertices[m.triangles[i + 1]], m.vertices[m.triangles[i + 2]]);
		if (min < division_no)
		{
			histogram[min] += area;
		}
	}
	for (size_t i = 0; i < histogram.size(); i++)
	{
		std::cout << (int)histogram[i] << std::endl;
	}
	// finally delete all
	delete[] minimum_distance;
	delete[] is_in_path;
	delete[] tau_chart;
	delete[] is_close;
}

//compute the bilateral map accoridng to a point 
static void compute_bilateral_map_according_to_point(TrilateralMesh& m, int point1_index, float tau, int division_no) //tau is the distance d here 
{
	//declare variable 
	std::vector<float> histogram;
	//init variable
	for (size_t i = 0; i < division_no; i++)
	{
		histogram.push_back(0);
	}

	float histogram_no = tau / division_no;
	std::vector<float> distances = Geodesic_dijkstra(m, point1_index);

	int* is_close = new int[m.vertices.size()];
	int* is_in_path = new int[m.vertices.size()];
	float* minimum_distance = new float[m.vertices.size()];
	float* tau_chart = new float[division_no];
	//init variable
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_close[i] = 0; //false for all 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_close[i] = 0; //false for all 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_in_path[i] = 0; //false for all 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		minimum_distance[i] = (float)INFINITY; //false for all 
	}
	for (size_t j = 0; j < distances.size(); j++)
	{
		if (distances[j] <= tau)
		{
			is_close[j] = 1;
		}
	}
	for (size_t i = 0; i < m.triangles.size() / 3; i = i + 3)
	{
		float dist1 = distances[m.triangles[i]];
		float dist2 = distances[m.triangles[i + 1]];
		float dist3 = distances[m.triangles[i + 2]];

		int no_1 = dist1 / histogram_no;
		int no_2 = dist2 / histogram_no;
		int no_3 = dist3 / histogram_no;

		int min = INFINITY;
		if (no_1 < min)
		{
			min = no_1;
		}
		if (no_2 < min)
		{
			min = no_2;
		}
		if (no_3 < min)
		{
			min = no_3;
		}
		//put the traingles into min index of histogram
		float area = compute_triangle_area(m.vertices[m.triangles[i]], m.vertices[m.triangles[i + 1]], m.vertices[m.triangles[i + 2]]);
		if (min < division_no)
		{
			histogram[min] += area;
		}
	}
	//buffer 
	std::vector<float> new_buffer;
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		//now, if it is in the path continue in red
		// if not path adn tau blue
		// else black
		new_buffer.push_back(m.vertices[i].x);
		new_buffer.push_back(m.vertices[i].y);
		new_buffer.push_back(m.vertices[i].z);
		if (is_close[i] == 1)
		{
			if (i == point1_index)
			{
				//red
				new_buffer.push_back(1.0f);
				new_buffer.push_back(0.0f);
				new_buffer.push_back(0.0f);

			}
			else
			{
				int div_no = distances[i] / division_no;
				if (div_no % 2 == 0)
				{
					//blue 
					new_buffer.push_back(0.0f);
					new_buffer.push_back(0.0f);
					new_buffer.push_back(1.0f);
				}
				else
				{
					//green 
					new_buffer.push_back(0.0f);
					new_buffer.push_back(1.0f);
					new_buffer.push_back(0.0f);
				}

			}
		}
		else
		{
			//black 
			new_buffer.push_back(0.0f);
			new_buffer.push_back(0.0f);
			new_buffer.push_back(0.0f);
		}

	}
	for (size_t i = 0; i < histogram.size(); i++)
	{
		std::cout << (int)histogram[i] << std::endl;
	}
	glBufferData(GL_ARRAY_BUFFER, new_buffer.size() * sizeof(float), &new_buffer[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);
	// now we have the minimum distances 

}
#pragma endregion assignment part 3 */
