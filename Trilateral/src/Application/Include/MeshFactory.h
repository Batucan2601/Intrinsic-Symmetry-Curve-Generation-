#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <vector>
#include <glm/ext/vector_float3.hpp>
#include <glm/geometric.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "happly.h"
#include "../Include/Mesh.h"
#include "../Include/FibonacciHeap.h"

#define INFINITE 10000000


#pragma region functions
static float compute_triangle_area(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 vec1 = p2 - p1;
	glm::vec3 vec2 = p3 - p1;

	glm::vec3 cross1 = glm::cross( vec1 ,  vec2 );
	float length = glm::length(cross1) / 2;
	return length;
}
static std::vector<float> compute_geodesic_distances_fibonacci_heap_distances(Mesh& m, int point_index)
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
			//infinite
			matrix[i] = (float)INFINITE;
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
static std::vector<int> compute_geodesic_distances_fibonacci_heap(Mesh& m, int point_index)
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
			//infinite
			matrix[i] = (float)INFINITE;
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
static std::vector<int> draw_with_fib_heap_implementation(Mesh& m, int p1_index, int p2_index)
{
	std::vector<int> predecessors = compute_geodesic_distances_fibonacci_heap(m, p1_index);
	std::vector<int> consec_indices;
	int pred = p2_index;
	while (pred != p1_index )
	{
		consec_indices.push_back(pred);
		pred = predecessors[pred];
	}
	return consec_indices;
}
typedef struct
{
	std::vector<float> point_pairs;
	glm::mat4 model_mat;
}MeshPointPairs;
#pragma endregion functions
class MeshFactory
{
public:
	MeshFactory();
	~MeshFactory();
	void add_mesh(Mesh &m);
	void buffer_meshes();
	void draw_meshes();
	void draw_mesh(int mesh_no);
	void draw_lines();
	void create_lines(int mesh_no1 , int mesh_no2 , int  no_of_lines);
	void remove_mesh(int mesh_no);
	void remove_all();
	void add_all();
	void add_lines_between_meshes(int mesh_index1, int mesh_index2, int p1, int p2);
	void get_camera_and_projection(glm::mat4 view_, glm::mat4  projection_);
	void remove_lines();
	std::vector<Mesh> mesh_vec;
	std::vector<MeshPointPairs> mesh_point_pairs; //point paris for each mesh
	std::vector<glm::vec3> points;
	std::vector<glm::vec3> colors;
	std::vector<int> point_indices;
	std::vector<std::pair<int, int>> correspondence_lines;
	//skeleton part 

	
	
	// for skeleton detection from pose  video-to-pose
	std::vector<float> mesh_skeleton_vec;
	unsigned int skeleton_VAO; 
	// 
	glm::mat4 view;
	glm::mat4 projection;
	// line 
private:
	
};

