#include "../Include/MeshFactory.h"
#include <GL/glew.h>
#include <utility>
#include <algorithm>
#include <queue>
#include <stack>
#include <algorithm> 
#include <eigen/Eigen/Dense>
#include "Sampling.h"
#include "CoreTypeDefs.h"

#pragma once
static std::vector<float> compute_geodesic_distances_min_heap_distances(Mesh& m, int point_index)
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

static void trilateral_map(MeshFactory &mesh_fac ,  int &selected_index  ,int p1, int p2, int p3)
{
	Mesh m = mesh_fac.mesh_vec[selected_index];
	//extract two paths
	//1 - extract the path from p1  to p2  and p1 to p3  and p2 to p3 
	std::vector<int> path_1_2 = draw_with_fib_heap_implementation(m, p1 , p2 );
	std::vector<int> path_1_3 = draw_with_fib_heap_implementation(m, p1, p3  );
	std::vector<int> path_2_3 = draw_with_fib_heap_implementation(m, p2, p3  );

#pragma region caluclation of distances between points 2 ,3 and 1  

	float dist_1_2 = 0.0f ;
	float dist_1_3 = 0.0f;
	float dist_2_3 = 0.0f;

	for (size_t i = 0; i < path_1_2.size() - 1 ; i++)
	{
		dist_1_2 += glm::distance(m.vertices[path_1_2[i]], m.vertices[path_1_2[i + 1]]);
	}
	for (size_t i = 0; i < path_1_3.size() - 1; i++)
	{
		dist_1_3 += glm::distance(m.vertices[path_1_3[i]], m.vertices[path_1_3[i + 1]]);
	}
	for (size_t i = 0; i < path_2_3.size() - 1; i++)
	{
		dist_2_3 += glm::distance(m.vertices[path_2_3[i]], m.vertices[path_2_3[i + 1]]);
	}
	std::cout << "dist1 " << dist_1_2 << std::endl;
	std::cout << "dist 2 " << dist_1_3 << std::endl;
#pragma endregion 

#pragma region  angle
	// for p1
	float angle_p1 = 0.0f; // in degrees ( from angle betwen 2 , 1 ,3 )
	glm::vec3 p1_point = glm::vec3(m.vertices[p1]);
	glm::vec3 p2_point = glm::vec3(m.vertices[p2]);
	glm::vec3 p3_point = glm::vec3(m.vertices[p3]);
	glm::vec3 vec_2_1 = glm::vec3(p2_point - p1_point);
	glm::vec3 vec_3_1 = glm::vec3(p3_point - p1_point);
	double angle_radian_p1 =  glm::acos( glm::dot(vec_2_1 , vec_3_1 ) / (glm::length(vec_2_1) * glm::length(vec_3_1)));
	angle_p1 = glm::degrees(angle_radian_p1); //not invariant 

	std::cout << "angle for point p1 == " << angle_p1 << std::endl;

	// for p2 
	float angle_p2 = 0.0f; // in degrees ( from angle betwen 2 , 1 ,3 )
	glm::vec3 vec_1_2 = glm::vec3(p1_point - p2_point);
	glm::vec3 vec_3_2 = glm::vec3(p3_point - p2_point);

	double angle_radian_p2 = glm::acos(glm::dot(vec_1_2, vec_3_2) / (glm::length(vec_1_2) * glm::length(vec_3_2)));
	angle_p2 = glm::degrees(angle_radian_p2); //not invariant 

	std::cout << "angle for point p2 == " << angle_p2 << std::endl;

	// for p3
	float angle_p3 = 0.0f; // in degrees ( from angle betwen 2 , 1 ,3 )
	glm::vec3 vec_1_3 = glm::vec3(p1_point - p3_point);
	glm::vec3 vec_2_3 = glm::vec3(p2_point - p3_point);

	double angle_radian_p3 = glm::acos(glm::dot(vec_1_3, vec_2_3) / (glm::length(vec_1_3) * glm::length(vec_2_3)));
	angle_p3 = glm::degrees(angle_radian_p3); //not invariant 

	std::cout << "angle for point p3 == " << angle_p3 << std::endl;

#pragma endregion 
	
	std::vector<float> draw_buffer; 
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		glm::vec3 point = m.vertices[i];
		bool is_in_way = false; 
		for (size_t j = 0; j < path_1_2.size(); j++)
		{
			glm::vec3 point_between_1_2 = m.vertices[path_1_2[j]];

			if (point.x == point_between_1_2.x && point.y == point_between_1_2.y && point.z == point_between_1_2.z)
			{
				is_in_way = true; 
			}
		}
		for (size_t j = 0; j < path_1_3.size(); j++)
		{
			glm::vec3 point_between_1_3 = m.vertices[path_1_3[j]];

			if (point.x == point_between_1_3.x && point.y == point_between_1_3.y && point.z == point_between_1_3.z)
			{
				is_in_way = true;
			}
		}
		for (size_t j = 0; j < path_2_3.size(); j++)
		{
			glm::vec3 point_between_2_3 = m.vertices[path_2_3[j]];

			if (point.x == point_between_2_3.x && point.y == point_between_2_3.y && point.z == point_between_2_3.z)
			{
				is_in_way = true;
			}
		}
		// rebuffer data
		draw_buffer.push_back(point.x );
		draw_buffer.push_back(point.y);
		draw_buffer.push_back(point.z );

		//red 
		if (is_in_way)
		{
			if (i == p1 || i == p2 || i == p3) //point itself  (will draw red therefore skip )
			{
				draw_buffer.push_back(0.0f);
			}
			else
			{
				draw_buffer.push_back(1.0f); // on the way 

			}

			
		}
		else
		{
			draw_buffer.push_back(0.0f);
		}

		//green 
		draw_buffer.push_back(0.0f);


		//blue
		if (i == p1 || i == p2 || i == p3) //point itself  
		{
			draw_buffer.push_back(1.0f);
		}
		else
		{
			draw_buffer.push_back(0.0f);
		}



	}
	int point_size = 0;
	for (size_t i = 0; i < selected_index; i++)
	{
		point_size += mesh_fac.mesh_vec[i].vertices.size() * 6;
	}
	glBufferSubData(GL_ARRAY_BUFFER, point_size * sizeof(float), draw_buffer.size() * sizeof(float), &draw_buffer[0]);
	/*glBufferData(GL_ARRAY_BUFFER, draw_buffer.size() * sizeof(float), &draw_buffer[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);*/

}


/*static float* trialteral_ROI(MeshFactory& mesh_fac  , int& selected_index, int point_index1, int point_index2, int point_index3, int division_no) // tau is the closeness division_no is the no of how much you want to separate
{
	Mesh m = mesh_fac.mesh_vec[selected_index];
	int* is_close = new int[m.vertices.size()];
	int* is_in_path = new int[m.vertices.size()];
	float* minimum_distance = new float[m.vertices.size()];
	float* tau_chart = new float[division_no];
	//0 - find tau
	float maximum_length_of_minimums = 0.0f;
	// find tau 
	glm::vec3 vec1_ = glm::vec3(m.vertices[point_index2] - m.vertices[point_index1]);
	glm::vec3 vec2_ = glm::vec3(m.vertices[point_index3] - m.vertices[point_index1]);
	glm::vec3 vec3_ = glm::vec3(m.vertices[point_index3] - m.vertices[point_index2]);
	//get max edge  and use it as perimeter 
	if (glm::length(vec1_) >= glm::length(vec2_) && glm::length(vec1_) >= glm::length(vec3_))
	{
		maximum_length_of_minimums = glm::length(vec1_);
	}
	else if (glm::length(vec2_) >= glm::length(vec1_) && glm::length(vec2_) >= glm::length(vec3_))
	{
		maximum_length_of_minimums = glm::length(vec2_);
	}
	else
	{
		maximum_length_of_minimums = glm::length(vec3_);
	}
	float total = 0.0f;
	float tau = maximum_length_of_minimums;
	for (size_t i = 0; i < division_no; i++)
	{
		tau_chart[i] = total;
		float step = maximum_length_of_minimums / division_no;
		total += step;
	}
	//1 - extract the path from point1 to point2 
	std::vector<int> path_1_2 = draw_with_fib_heap_implementation(m, point_index1, point_index2);
	std::vector<int> path_1_3 = draw_with_fib_heap_implementation(m, point_index1, point_index3);
	std::vector<int> path_2_3 = draw_with_fib_heap_implementation(m, point_index2, point_index3);

	std::vector<float> path_1 = compute_geodesic_distances_min_heap_distances(m, point_index1);
	std::vector<float> path_2 = compute_geodesic_distances_min_heap_distances(m, point_index2);
	std::vector<float> path_3 = compute_geodesic_distances_min_heap_distances(m, point_index3);

	//2 - we should calculate the geodesic distances for all of the path and than redraw
	//declare variable 
	std::vector<float> histogram;
	//init variable
	for (size_t i = 0; i < division_no; i++)
	{
		histogram.push_back(0);
	}
	float total_area = compute_triangle_area(m.vertices[point_index1], m.vertices[point_index2], m.vertices[point_index3]);


	//init variable
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_close[i] = 0; //false for all 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		is_in_path[i] = 0; //false for all 
	}
	


	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		is_in_path[path_1_2[i]] = true; //get the path ones 
	}
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		is_in_path[path_1_3[i]] = true; //get the path ones 
	}
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		is_in_path[path_2_3[i]] = true; //get the path ones 
	}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		minimum_distance[i] = INFINITE;

	}
	//calculate the distance from every vertex within path
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(m, path_1_2[i]);
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
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(m, path_2_3[i]);
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
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(m, path_1_3[i]);
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
				//check if it is inside
				
				
				// first we need to project the point in the plane of  our triangle 
				// 1 - generate plane
				// two vector 
				glm::vec3 vec1  = glm::vec3(m.vertices[point_index2] - m.vertices[point_index1]); 
				glm::vec3 vec2  = glm::vec3(m.vertices[point_index3] - m.vertices[point_index1]);
				//cross those vectors
				glm::vec3 normal_vec =glm::normalize(glm::cross(vec1, vec2));

				glm::vec3 vec_to_point_from_plane = glm::vec3(m.vertices[i] -m.vertices[point_index1]);
				float dot_prod = glm::dot(vec_to_point_from_plane , normal_vec);
				glm::vec3 projected_point = m.vertices[i] - dot_prod * normal_vec;

				float inside_area_1 = compute_triangle_area(projected_point, m.vertices[point_index1], m.vertices[point_index2]);
				float inside_area_2 = compute_triangle_area(projected_point, m.vertices[point_index1], m.vertices[point_index3]);
				float inside_area_3 = compute_triangle_area(projected_point, m.vertices[point_index2], m.vertices[point_index3]);
				// 2 - now we have a normal and 3 points 

				//before calculating the area we have to do another mesh search for disposing the meshes from other side : ( 
				// via ray casting

				

				bool projection_check = glm::abs((inside_area_1 + inside_area_2 + inside_area_3) - total_area) <= 1e-5; //area 
				bool length_check = path_1[i] + path_2[i] + path_3[i] <= glm::length(vec1_) + glm::length(vec2_) + glm::length(vec3_); 
				if(length_check && projection_check)
				//if (glm::abs((inside_area_1 + inside_area_2 + inside_area_3) - total_area) <= 1e-5 )
				{
					
					for (size_t j = 0; j < division_no - 1 ; j++)
					{
						if (tau_chart[j] <= minimum_distance[i] && tau_chart[j + 1] > minimum_distance[i])
						{
							if (j % 2 == 0)
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
							break;
						}
					}
				}
				else //if not inside leave it black
				{
					
					new_buffer.push_back(0.0f);
					new_buffer.push_back(0.0f);
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
	int point_size = 0;
	for (size_t i = 0; i < selected_index; i++)
	{
		point_size += mesh_fac.mesh_vec[i].vertices.size() * 6;
	} 

	glBufferSubData(GL_ARRAY_BUFFER, point_size * sizeof(float), new_buffer.size() * sizeof(float) , &new_buffer[0]);
	/*glBufferData(GL_ARRAY_BUFFER, new_buffer.size() * sizeof(float), &new_buffer[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);

	for (size_t i = 0; i < m.triangles.size() / 3; i = i + 3)
	{
		if (is_close[m.triangles[i]] + is_close[m.triangles[i + 1]] + is_close[m.triangles[i + 2]] == 3)
		{
			float distances[3];
			distances[0] = minimum_distance[m.triangles[i]];
			distances[1] = minimum_distance[m.triangles[i + 1]];
			distances[2] = minimum_distance[m.triangles[i + 2]];

			float min_dist = INFINITY;
			for (size_t j = 0; j < 3; j++)
			{
				if (min_dist > distances[j])
				{
					min_dist = distances[j];
				}
			}

			float hist_index = min_dist / (tau / division_no);
			float area = compute_triangle_area(m.vertices[m.triangles[i]], m.vertices[m.triangles[i + 1]], m.vertices[m.triangles[i + 2]]);
			histogram[hist_index] += area;
		}
		
	}
	for (size_t i = 0; i < histogram.size(); i++)
	{
		std::cout << (float)histogram[i] << std::endl;
	}
	// finally delete all
	delete[] minimum_distance;
	delete[] is_in_path;
	delete[] is_close;
	return tau_chart; 
}*/



struct TrilateralDescriptor
{
	double area; // ROI
	double total_length;
	double lenght_1_2;//  geodesic length between 1 - 2
	double lenght_1_3;//  length between 1 - 3
	double lenght_2_3;//  length between 2 - 3
	double curvature_1_2; // (geodesic / euclidiean)
	double curvature_1_3;
	double curvature_2_3;
};


static int* trialteral_ROI(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no , bool & is_visited_interior )
{
	Mesh *m = &mesh_fac.mesh_vec[selected_index];
	std::vector<int> path_1_2 = draw_with_fib_heap_implementation(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = draw_with_fib_heap_implementation(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = draw_with_fib_heap_implementation(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index1);
	std::vector<float> distance_matrix_p2 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index2);
	std::vector<float> distance_matrix_p3 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index3);
	// do a breadth first search
	//	get a random number of a vertex from a mesh 
	int random_vertex_index = rand() % m->vertices.size(); // we know that this is inside the system

	
	
	//now check if this index is equal to point1 point2 or point3 
	while ( (random_vertex_index == point_index1 || random_vertex_index == point_index2 || random_vertex_index == point_index3) )
	{
		random_vertex_index += 1;
	}
	//if not equal continue 
	
	//create a stack for BFS ( breadth first search )
	std::stack<int> stack;  // a stack consisting of indices

	// get the adjacencies
	std::vector<std::vector<std::pair<int, float>>> mesh_adjacencies = m->adjacenies;

	//lastly get a int array with  size of vertices in order to check if the vertex has been visited ( -1 edge , 0 not visisted , 1 visited) 
	int* is_visited = new int[m->vertices.size()];
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		is_visited[i] = 0;
	}
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		is_visited[path_1_2[i]] = -1;
	}
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		is_visited[path_1_3[i]] = -1;

	}
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		is_visited[path_2_3[i]] = -1;
	}

	//push our point to stack
	stack.push(random_vertex_index);
	while (!stack.empty())
	{
		int index = stack.top();
		stack.pop(); //vertex index popped from stack
		if (is_visited[index] == 0 ) //not visited
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



	for (size_t i = 0; i < m->vertices.size(); i++)  //start checking with visited area , it is highly likelty that visited is outside
	{
		if (is_visited[i] == 1 )
		{
			for (size_t j = i; j < m->vertices.size(); j++)
			{
				if (i != j)
				{
					std::vector<int> path_i_j = draw_with_fib_heap_implementation(*m, i, j);

					for (size_t k = 0; k < path_i_j.size(); k++)
					{
						for (size_t t = 0; t < path_1_2.size(); t++)
						{
							if (path_i_j[k] == path_1_2[t])
							{
								is_visited_interior = true;
							}
						}
						for (size_t t = 0; t < path_1_3.size(); t++)
						{
							if (path_i_j[k] == path_1_3[t])
							{
								is_visited_interior = true;
							}
						}
						for (size_t t = 0; t < path_2_3.size(); t++)
						{
							if (path_i_j[k] == path_2_3[t])
							{
								is_visited_interior = true;
							}
						}
						if (is_visited_interior)
						{
							break; 
						}
					}

					if (is_visited_interior)
					{
						break;
					}
				}
			}
			if (is_visited_interior)
			{
				break;
			}
		}
	}
	 // is_visited_interior = true is_visited 0 is the correct interior 
	 // is_visited_interior = false is_visited 1 is the correct interior 
	// now recolor
	std::vector<glm::vec3> new_color_buffer;
	if (is_visited_interior)
	{
		for (size_t i = 0; i < m->colors.size(); i++)
		{

			if (is_visited[i] == -1) //edge 
			{
				new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 1) 
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 0)
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
			}
		}
	}
	else
	{
		for (size_t i = 0; i < m->colors.size(); i++)
		{

			if (is_visited[i] == -1) //edge 
			{
				new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 0) //not visited 
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 1)
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
			}
		}
	}
	m->colors = new_color_buffer;
	return  is_visited;
}
std::vector<float>  histogramROi(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no, int* is_visited, bool is_visited_interior)
{
	//histogram to be returned 
	std::vector<float> histogram; 
	// fill it with division no 
	for (size_t i = 0; i < division_no+1; i++)
	{
		histogram.push_back(0);
	}
	Mesh* m = &mesh_fac.mesh_vec[selected_index];
	std::vector<int> path_1_2 = draw_with_fib_heap_implementation(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = draw_with_fib_heap_implementation(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = draw_with_fib_heap_implementation(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index1);
	std::vector<float> distance_matrix_p2 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index2);
	std::vector<float> distance_matrix_p3 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index3);
	//find the maximum distance from is_visited and paths
	float max = -999;
	float min = 100000;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (is_visited[i] == 1) // if vertex is visited 
		{
			for (size_t j = 0; j < path_1_2.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_1_2[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_1_2[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_1_2[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_1_2[j]]);
				}
			}
			for (size_t j = 0; j < path_1_3.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_1_3[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_1_3[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_1_3[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_1_3[j]]);
				}
			}
			for (size_t j = 0; j < path_2_3.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_2_3[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_2_3[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_2_3[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_2_3[j]]);
				}
			}
		}
	}
	float step = max / division_no;

	//now recolor
	std::vector<glm::vec3> new_color_buffer;
	if (is_visited_interior)
	{
		for (size_t i = 0; i < m->colors.size(); i++)
		{

			if (is_visited[i] == -1) //edge 
			{
				new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 1) //not visited 
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 0) // get the max distance
			{
				//get distance
				float max_dist_from_index_i = std::max(distance_matrix_p1[i], distance_matrix_p2[i]);
				max_dist_from_index_i = std::max(max_dist_from_index_i, distance_matrix_p3[i]);
				int hist_index = max_dist_from_index_i / step;
				if (hist_index % 2 == 0)
				{
					new_color_buffer.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

				}
				else if (hist_index % 2 == 1)
				{
					new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
				}
			}
		}
	}
	else
	{
		for (size_t i = 0; i < m->colors.size(); i++)
		{
			if (is_visited[i] == -1) //edge 
			{
				new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 0) //not visited 
			{
				new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			}
			else if (is_visited[i] == 1) // get the max distance
			{
				//get distance
				float max_dist_from_index_i = std::max(distance_matrix_p1[i], distance_matrix_p2[i]);
				max_dist_from_index_i = std::max(max_dist_from_index_i, distance_matrix_p3[i]);
				int hist_index = max_dist_from_index_i / step;
				if (hist_index % 2 == 0)
				{
					new_color_buffer.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

				}
				else if (hist_index % 2 == 1)
				{
					new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
				}
			}
		}
	}
	m->colors = new_color_buffer;

	//getting the numbers
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		if (is_visited[m->triangles[i]] != 0 || is_visited[m->triangles[i + 1 ] ] != 0  || is_visited[m->triangles[i + 2 ] ] != 0) //if any vertex is visited
		{
			glm::vec3 p1 = m->vertices[m->triangles[i]];
			glm::vec3 p2 = m->vertices[m->triangles[i + 1]];
			glm::vec3 p3 = m->vertices[m->triangles[i + 2]];

			float min_dist_from_p1 = std::min(distance_matrix_p1[m->triangles[i]], distance_matrix_p2[m->triangles[i]]);
			min_dist_from_p1 = std::min(min_dist_from_p1, distance_matrix_p3[m->triangles[i]]);

			float min_dist_from_p2 = std::min(distance_matrix_p1[m->triangles[i+1]], distance_matrix_p2[m->triangles[i + 1]]);
			min_dist_from_p2 = std::min(min_dist_from_p2, distance_matrix_p3[m->triangles[i + 1]]);

			float min_dist_from_p3 = std::min(distance_matrix_p1[m->triangles[i+2]], distance_matrix_p2[m->triangles[i+2]]);
			min_dist_from_p3 = std::min(min_dist_from_p3, distance_matrix_p3[m->triangles[i+2]]);

			p1 *= 1e5;
			p2 *= 1e5;
			p3 *= 1e5;
			float area = compute_triangle_area(p1 , p2, p3);

			int hist_no_p1 = min_dist_from_p1 / step; 
			int hist_no_p2 = min_dist_from_p2 / step; 
			int hist_no_p3 = min_dist_from_p3 / step;

			if (hist_no_p1 > division_no)
			{
				hist_no_p1 = division_no;
			}
			if (hist_no_p2 > division_no)
			{
				hist_no_p2 = division_no;
			}
			if (hist_no_p3 > division_no)
			{
				hist_no_p3 = division_no;
			}
			area /= (3  * 1e5);
			histogram[hist_no_p1] += area;
			histogram[hist_no_p2] += area;
			histogram[hist_no_p3] += area;
		}
	}
	return histogram; 
}

// the ultimate trilateral descriptor generator 
static TrilateralDescriptor  generate_trilateral_descriptor(  MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, bool is_simplified)
{
	Mesh* m = &mesh_fac.mesh_vec[selected_index];

	TrilateralDescriptor trilateral_descriptor;//trialteral descriptor 
	//init descriptor
	trilateral_descriptor.area = 0;
	trilateral_descriptor.lenght_1_2 = 0;
	trilateral_descriptor.lenght_1_3 = 0;
	trilateral_descriptor.lenght_2_3 = 0;
		

	std::vector<int> path_1_2 = draw_with_fib_heap_implementation(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = draw_with_fib_heap_implementation(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = draw_with_fib_heap_implementation(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index1);
	std::vector<float> distance_matrix_p2 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index2);
	std::vector<float> distance_matrix_p3 = compute_geodesic_distances_fibonacci_heap_distances(*m, point_index3);

	// get distances
	trilateral_descriptor.lenght_1_2 = distance_matrix_p1[point_index2];
	trilateral_descriptor.lenght_1_3 = distance_matrix_p1[point_index3];
	trilateral_descriptor.lenght_2_3 = distance_matrix_p2[point_index3];

	trilateral_descriptor.curvature_1_2 = trilateral_descriptor.lenght_1_2 / glm::distance(m->vertices[point_index1] , m->vertices[point_index2]);
	trilateral_descriptor.curvature_1_3 = trilateral_descriptor.lenght_1_3 /  glm::distance(m->vertices[point_index1] , m->vertices[point_index3]);
	trilateral_descriptor.curvature_2_3 = trilateral_descriptor.lenght_2_3 / glm::distance(m->vertices[point_index2] , m->vertices[point_index3]);
	
	//for only brute force research 
	if (is_simplified)
	{
		return trilateral_descriptor;
	}
	int random_vertex_index = rand() % m->vertices.size(); // we know that this is inside the system

	bool is_visited_interior = false;

	//now check if this index is equal to point1 point2 or point3 
	while ((random_vertex_index == point_index1 || random_vertex_index == point_index2 || random_vertex_index == point_index3))
	{
		random_vertex_index += 1;
	}
	//if not equal continue 
	//create a stack for BFS ( breadth first search )
	std::stack<int> stack;  // a stack consisting of indices
	// get the adjacencies
	std::vector<std::vector<std::pair<int, float>>> mesh_adjacencies = m->adjacenies;
	//lastly get a int array with  size of vertices in order to check if the vertex has been visited ( -1 edge , 0 not visisted , 1 visited) 
	int* is_visited = new int[m->vertices.size()];
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		is_visited[i] = 0;
	}
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		is_visited[path_1_2[i]] = -1;
	}
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		is_visited[path_1_3[i]] = -1;
	}
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		is_visited[path_2_3[i]] = -1;
	}

	//push our point to stack
	stack.push(random_vertex_index);
	while (!stack.empty())
	{
		int index = stack.top();
		stack.pop(); //vertex index popped from stack
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



	for (size_t i = 0; i < m->vertices.size(); i++)  //start checking with visited area , it is highly likelty that visited is outside
	{
		if (is_visited[i] == 1)
		{
			for (size_t j = i; j < m->vertices.size(); j++)
			{
				if (i != j)
				{
					std::vector<int> path_i_j = draw_with_fib_heap_implementation(*m, i, j);

					for (size_t k = 0; k < path_i_j.size(); k++)
					{
						for (size_t t = 0; t < path_1_2.size(); t++)
						{
							if (path_i_j[k] == path_1_2[t])
							{
								is_visited_interior = true;
							}
						}
						for (size_t t = 0; t < path_1_3.size(); t++)
						{
							if (path_i_j[k] == path_1_3[t])
							{
								is_visited_interior = true;
							}
						}
						for (size_t t = 0; t < path_2_3.size(); t++)
						{
							if (path_i_j[k] == path_2_3[t])
							{
								is_visited_interior = true;
							}
						}
						if (is_visited_interior)
						{
							break;
						}
					}

					if (is_visited_interior)
					{
						break;
					}
				}
			}
			if (is_visited_interior)
			{
				break;
			}
		}
	}


	//area 
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		if (is_visited[m->triangles[i]] != 0 || is_visited[m->triangles[i + 1]] != 0 || is_visited[m->triangles[i + 2]] != 0) //if any vertex is visited
		{
			glm::vec3 p1 = m->vertices[m->triangles[i]];
			glm::vec3 p2 = m->vertices[m->triangles[i + 1]];
			glm::vec3 p3 = m->vertices[m->triangles[i + 2]];

			float min_dist_from_p1 = std::min(distance_matrix_p1[m->triangles[i]], distance_matrix_p2[m->triangles[i]]);
			min_dist_from_p1 = std::min(min_dist_from_p1, distance_matrix_p3[m->triangles[i]]);

			float min_dist_from_p2 = std::min(distance_matrix_p1[m->triangles[i + 1]], distance_matrix_p2[m->triangles[i + 1]]);
			min_dist_from_p2 = std::min(min_dist_from_p2, distance_matrix_p3[m->triangles[i + 1]]);

			float min_dist_from_p3 = std::min(distance_matrix_p1[m->triangles[i + 2]], distance_matrix_p2[m->triangles[i + 2]]);
			min_dist_from_p3 = std::min(min_dist_from_p3, distance_matrix_p3[m->triangles[i + 2]]);

			p1 *= 1e5;
			p2 *= 1e5;
			p3 *= 1e5;
			float area = compute_triangle_area(p1, p2, p3);



			area /= (3 * 1e5);
			trilateral_descriptor.area += area; // get area
		}
	}

	return trilateral_descriptor;
}
static void point_matching_with_dominant_symmetry_plane(MeshFactory& mesh_fac, int& selected_index, Plane* plane  , int sampling_no  )
{
	Mesh* mesh = &mesh_fac.mesh_vec[selected_index];
	
	std::vector<int> vertex_indices = furthest_point_sampling(mesh,sampling_no);
	// divide the points in two 
	std::vector<int> vertex_indices_right;
	std::vector<int> vertex_indices_left;

	float* vertex_indices_left_right = new float[vertex_indices.size() ];
	for (size_t i = 0; i < vertex_indices.size(); i++)
	{
		//check left or right
		float point_loc = get_point_status_from_plane(plane, &mesh->vertices[vertex_indices[i]]);

		vertex_indices_left_right[i] = point_loc;

	}
	bool* vertex_indices_location_array = new bool[mesh->vertices.size()];
	for (size_t i = 0; i < mesh->vertices.size(); i++)
	{
		bool is_in_indices = false; 
		for (size_t j = 0; j < vertex_indices.size(); j++)
		{
			if (vertex_indices[j] == i)
			{
				is_in_indices = true; 
				break; 
			}
		}
		vertex_indices_location_array[i] = is_in_indices;
	}
	std::vector<std::vector<int>> closest_3_pairs; // n x 3 matrix
	for (size_t i = 0; i < vertex_indices.size(); i++)
	{
		std::vector<float> closest_matrix = compute_geodesic_distances_fibonacci_heap_distances(*mesh, vertex_indices[i]);
		//get closest 3 
		int closest_index_1 = -1;
		int closest_index_2 = -1;
		float closest_distance1 = 1e5; 
		float closest_distance2 = 1e5;
		for (size_t j = 0; j < vertex_indices.size() ; j++)
		{
				if (closest_distance1 > closest_matrix[vertex_indices[j]])
				{
					if (j != i) //because it is 0 an himself
					{
						if (((vertex_indices_left_right[i] >= 0) && (vertex_indices_left_right[j] >= 0)) || ((vertex_indices_left_right[vertex_indices[i]] < 0) && (vertex_indices_left_right[j] < 0))) //check same sign 
						{
							closest_index_1 = j;
							closest_distance1 = closest_matrix[j];
						}

					}
				}
			
		}
		for (size_t j = 0; j < vertex_indices.size(); j++)
		{
				if (closest_distance2 > closest_matrix[vertex_indices[j]])
				{
					if (j != i && j != closest_index_1 ) //because it is 0 an himself
					{
						if (((vertex_indices_left_right[i] >= 0) && (vertex_indices_left_right[j] >= 0)) || ((vertex_indices_left_right[i] < 0) && (vertex_indices_left_right[j] < 0))) //check same sign 
						{
							closest_index_2 = j;
							closest_distance2 = closest_matrix[j];
						}
					}
				}
			
		}
		std::vector<int> temp; 
		closest_3_pairs.push_back(temp);
		closest_3_pairs[i].push_back(vertex_indices[i]);
		closest_3_pairs[i].push_back(vertex_indices[closest_index_1]);
		closest_3_pairs[i].push_back(vertex_indices[closest_index_2]);
	}
	TrilateralDescriptor* descriptors = new TrilateralDescriptor[vertex_indices.size()];
	for (size_t i = 0; i < vertex_indices.size(); i++)
	{
		descriptors[i] = generate_trilateral_descriptor(mesh_fac , selected_index , closest_3_pairs[i][0] , closest_3_pairs[i][1] , closest_3_pairs[i][2] , true );
	}

	std::vector<std::vector<int>> vertex_index_pairs; 
	//last part match points.
	for (size_t i = 0; i < vertex_indices.size() ; i++)
	{
		float min_diff_score = INFINITY; 
		int min_diff_index = -1; 
		for (size_t j = 0; j < vertex_indices.size() ; j++)
		{
			if (i != j)
			{
				float area_dif = std::abs(descriptors[i].area - descriptors[j].area );
				float roi_len_dif = std::abs(descriptors[i].total_length - descriptors[j].total_length);
				// easy sum for now
				float curv_dif = std::abs(descriptors[i].curvature_1_2 + descriptors[i].curvature_1_3 + descriptors[i].curvature_2_3 - descriptors[j].curvature_1_2 - descriptors[j].curvature_1_3 - descriptors[j].curvature_2_3);
				float dist_dif = std::abs(descriptors[i].lenght_1_2+ descriptors[i].lenght_1_3+ descriptors[i].lenght_2_3 - descriptors[j].lenght_1_2- descriptors[j].lenght_1_3- descriptors[j].lenght_2_3);

				if (area_dif + roi_len_dif + curv_dif + dist_dif <  min_diff_score)
				{
					min_diff_score = area_dif + roi_len_dif + curv_dif + dist_dif;
					min_diff_index = j;
				}
			}
		}
		std::vector<int> temp;
		temp.push_back(i);
		temp.push_back(min_diff_index);
		vertex_index_pairs.push_back(temp);
	}
	delete[] vertex_indices_location_array;
	delete[] descriptors;
	delete[] vertex_indices_left_right;

	std::vector<float> pair_points;
	for (size_t i = 0; i < vertex_index_pairs.size(); i++)
	{
		pair_points.push_back(mesh->vertices[vertex_indices[i]].x);
		pair_points.push_back(mesh->vertices[vertex_indices[i]].y);
		pair_points.push_back(mesh->vertices[vertex_indices[i]].z);

		pair_points.push_back(0.0f);
		pair_points.push_back(1.0f);
		pair_points.push_back(1.0f);
	}
	MeshPointPairs pair;
	pair.point_pairs = pair_points; 
	pair.model_mat = mesh->model_mat;
	mesh_fac.mesh_point_pairs.push_back(pair);
}

// sampling
void simple_sample(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int sample_size, int division_no) // sample size must be bigger than 3 
{
	Mesh m1 = mesh_fac.mesh_vec[mesh_index1];
	Mesh m2 = mesh_fac.mesh_vec[mesh_index2];

	std::vector<int> sample_indices_m1;
	std::vector<int> sample_indices_m2;
	srand(time(NULL));
	for (size_t i = 0; i < sample_size; i++)
	{
		int point_m1 = rand() % m1.vertices.size();
		int point_m2 = rand() % m2.vertices.size();

		sample_indices_m1.push_back(point_m1);
		sample_indices_m2.push_back(point_m2);
	}

	std::vector<float> m1_histograms;
	std::vector<float> m2_histograms;
	// C(n,3) iteration
	// compare every vertex from m1 to every vertex from m2

	//m1
	for (size_t i = 0; i < sample_indices_m1.size(); i++)
	{
		for (size_t j = 0; j < sample_indices_m1.size(); j++)
		{
			if (i != j)
			{
				for (size_t k = 0; k < sample_indices_m1.size(); k++)
				{
					if (i != k && k != j)
					{
						float* histogram_arr = new float[division_no];
						//histogram_arr = trialteral_ROI(mesh_fac, mesh_index1, i, k, j, division_no); // tau is the closeness division_no is the no of how much you want to separate
						delete[] histogram_arr;
					}
				}
			}

		}
	}
}

//void match_points_from2_mesh(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int sample_size, int division_no) // sample size must be bigger than 3 
//{
//
//	Mesh* m1 = &mesh_fac.mesh_vec[mesh_index1];
//	Mesh* m2 = &mesh_fac.mesh_vec[mesh_index2];
//	srand(time(NULL));
//	std::vector<int> sample_indices_m1;
//	std::vector<int> sample_indices_m2;
//	srand(time(NULL));
//	for (size_t i = 0; i < sample_size; i++)
//	{
//		int point_m1 = rand() % m1->vertices.size();
//		int point_m2 = rand() % m2->vertices.size();
//
//		sample_indices_m1.push_back(point_m1);
//		sample_indices_m2.push_back(point_m2);
//	}
//	std::vector<std::pair<int , std::vector<std::vector<float>>>> m1_histograms;
//	
//	//creations 
//	std::vector<std::pair<int , std::vector<std::vector<float>>>> m2_histograms;
//	for (size_t i = 0; i < sample_size; i++) // mesh 1 
//	{
//		std::vector<std::vector<float>> float_vec; 
//		for (size_t j = 0; j < sample_size - 1 ; j++)
//		{
//			std::vector<float> histogram_vec;
//			float_vec.push_back(histogram_vec);
//		}
//		auto p1 = std::make_pair(sample_indices_m1[i], float_vec);
//		m1_histograms.push_back(p1);
//	}
//	for (size_t i = 0; i < sample_size; i++) // mesh 1 
//	{
//		std::vector<std::vector<float>> float_vec;
//		for (size_t j = 0; j < sample_size - 1; j++)
//		{
//			std::vector<float> histogram_vec;
//			float_vec.push_back(histogram_vec);
//		}
//		auto p1 = std::make_pair(sample_indices_m2[i], float_vec);
//		m2_histograms.push_back(p1);
//	}
//
//	//generate histograms
//	for (size_t i = 0; i < sample_size; i++) // for every point i 
//	{
//		for (size_t j = 0; j < sample_size ; j++) // for every point  j 
//		{
//			for (size_t k = 0; k < sample_size; k++) //for every point k 
//			{
//				if (i != j && i != k && j != k) //if points are different
//				{
//					bool is_visited_interior = false;
//					int* is_visited = trialteral_ROI(mesh_fac, mesh_index1, m1_histograms[i].first, m1_histograms[j].first, m1_histograms[k].first, division_no, is_visited_interior);
//					std::vector<float> histogram_vec =  histogramROi(mesh_fac, mesh_index1, m1_histograms[i].first, m1_histograms[j].first, m1_histograms[k].first, division_no, is_visited, is_visited_interior);
//					m1_histograms[i].second.push_back(histogram_vec);
//				}
//			}
//		}
//	}
//	for (size_t i = 0; i < sample_size; i++) // for every point i 
//	{
//		for (size_t j = 0; j < sample_size; j++) // for every point  j 
//		{
//			for (size_t k = 0; k < sample_size; k++) //for every point k 
//			{
//				if (i != j && i != k && j != k) //if points are different
//				{
//					bool is_visited_interior = false;
//					int* is_visited = trialteral_ROI(mesh_fac, mesh_index2, m2_histograms[i].first, m2_histograms[j].first, m2_histograms[k].first, division_no, is_visited_interior);
//					std::vector<float> histogram_vec = histogramROi(mesh_fac, mesh_index2, m2_histograms[i].first, m2_histograms[j].first, m2_histograms[k].first, division_no, is_visited, is_visited_interior);
//					m2_histograms[i].second.push_back(histogram_vec);
//				}
//			}
//		}
//	}
//	std::vector<std::vector<float>> m1; 
//	for (size_t i = 0; i < division_no; i++) // init 
//	{
//		std::vector<float> v1;
//		for (size_t j = 0; j < division_no; j++)
//		{
//			v1.push_back(0);
//		}
//		closeness_of_points.push_back(v1);
//	}
//
//
//	//comparison of each histogram
//	for(size_t i = 0; i < m1_histograms.size(); i++) // for every point
//	{
//		std::vector<float> total_histogram_m1;
//		float total_size_m1 = 0;
//		for (size_t j = 0; j < m1_histograms[i].second.size(); j++) //fill it with 0's 
//		{
//			total_histogram_m1.push_back(0);
//		}
//		for (size_t j = 0; j < m1_histograms[i].second.size(); j++) //sum up everyone		
//		{
//			std::vector<float> histogram  = m1_histograms[i].second[j];
//			for (size_t k = 0; k < division_no; k++)
//			{
//				total_histogram_m1[k] += histogram[k];
//				total_size_m1 += histogram[k];
//			}
//		}
//		for (size_t j = 0; j < m1_histograms[i].second.size(); j++) //normalize 		
//		{
//			total_histogram_m1[j] /= total_size_m1; 
//		}
//
//		//repeat for  histogram m2 
//		std::vector<float> total_histogram_m2;
//		float total_size_m2 = 0;
//		for (size_t j = 0; j < m2_histograms[i].second.size(); j++) //fill it with 0's 
//		{
//			total_histogram_m2.push_back(0);
//		}
//		for (size_t j = 0; j < m2_histograms[i].second.size(); j++) //sum up everyone		
//		{
//			std::vector<float> histogram = m2_histograms[i].second[j];
//			for (size_t k = 0; k < division_no; k++)
//			{
//				total_histogram_m2[k] += histogram[k];
//				total_size_m2 += histogram[k];
//			}
//		}
//		for (size_t j = 0; j < m2_histograms[i].second.size(); j++) //normalize 		
//		{
//			total_histogram_m2[j] /= total_size_m2;
//		}
//
//		closeness_of_points[i][j] = 
//	}
//}
void match_points_from2_mesh(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2,  int division_no) // sample size must be bigger than 3 
{
	//get the histogram of the traingle with the most near 2 points
	Mesh *m1 = &mesh_fac.mesh_vec[mesh_index1];
	Mesh *m2 = &mesh_fac.mesh_vec[mesh_index2];

	std::vector<int> sample_indices_m1;
	std::vector<int> sample_indices_m2;
	srand(time(NULL));
	int set[10] = {0 ,100 , 200 , 1500 , 2500,3200 , 5000 , 4000 , 400 ,4200 };
	for (size_t i = 0; i < division_no; i++)
	{
		int point_m1 = rand() % m1->vertices.size();
		int point_m2 = rand() % m2->vertices.size();

		sample_indices_m1.push_back(set[i]);
		sample_indices_m2.push_back(set[i]);
	}

	// for first mesh 
	std::vector<std::vector< float>> histograms_m1;
	for (size_t i = 0; i < sample_indices_m1.size(); i++) //get 2 nearest points 
	{
		std::vector<float> distances; // first is point index second is the distance 
		for (size_t j = 0; j < sample_indices_m1.size(); j++)//init and fill
		{
			distances.push_back(0);
		}
		int p1_index = sample_indices_m1[i];
		std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m1, p1_index);
		//traverse through distance matrix
		for (size_t j = 0; j < sample_indices_m1.size(); j++)
		{
			float distance = distance_matrix_p1[sample_indices_m1[j]];
			distances[j] = distance;
		}
			
		int first_min = -1; 
		int second_min = -1 ;

		float min = 999999;
		for (size_t j = 0; j < distances.size() ; j++)
		{
			if (j != i)
			{
				if (distances[j] < min)
				{
					min = distances[j];
					first_min = j;
				}
			}
		}

		min = 999999;
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (j != i && j != first_min)
			{
				if (distances[j] < min)
				{
					min = distances[j];
					second_min = j;
				}
			}
		}
		bool is_visited_interior = false;
		int* is_visited = trialteral_ROI(mesh_fac, mesh_index2, sample_indices_m1[i], sample_indices_m1[first_min], sample_indices_m1[second_min], division_no, is_visited_interior);
		std::vector<float> histogram_of_nearest_2_points = histogramROi(mesh_fac, mesh_index1, sample_indices_m1[i], sample_indices_m1[first_min], sample_indices_m1[second_min], division_no , is_visited , is_visited_interior);
		
		//normalize
		float total = 0;
		for (size_t j = 0; j < histogram_of_nearest_2_points.size(); j++)
		{
			total += histogram_of_nearest_2_points[j];
		}
		for (size_t j = 0; j < histogram_of_nearest_2_points.size(); j++)
		{
			histogram_of_nearest_2_points[j] /= total;
		}

		histograms_m1.push_back(histogram_of_nearest_2_points);
	}

	// for second mesh 
	std::vector<std::vector< float>> histograms_m2;
	for (size_t i = 0; i < sample_indices_m2.size(); i++) //get 2 nearest points 
	{
		std::vector<float> distances; // first is point index second is the distance 
		for (size_t j = 0; j < sample_indices_m2.size(); j++)//init and fill
		{
			distances.push_back(0);
		}
		int p1_index = sample_indices_m2[i];
		std::vector<float> distance_matrix_p1 = compute_geodesic_distances_fibonacci_heap_distances(*m1, p1_index);
		//traverse through distance matrix
		for (size_t j = 0; j < sample_indices_m1.size(); j++)
		{
			float distance = distance_matrix_p1[sample_indices_m2[j]];
			distances[j] = distance;
		}

		int first_min = -1;
		int second_min = -1;

		float min = 999999;
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (j != i)
			{
				if (distances[j] < min)
				{
					min = distances[j];
					first_min = j;
				}
			}
		}

		min = 999999;
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (j != i)
			{
				if (distances[j] < min)
				{
					min = distances[j];
					second_min = j;
				}
			}
		}
		bool is_visited_interior = false;
		int* is_visited = trialteral_ROI(mesh_fac, mesh_index2, sample_indices_m2[i], sample_indices_m2[first_min], sample_indices_m2[second_min], division_no, is_visited_interior);
		std::vector<float> histogram_of_nearest_2_points = histogramROi(mesh_fac, mesh_index1, sample_indices_m2[i], sample_indices_m2[first_min], sample_indices_m2[second_min], division_no, is_visited, is_visited_interior);

		//normalize
		float total = 0;
		for (size_t j = 0; j < histogram_of_nearest_2_points.size(); j++)
		{
			total += histogram_of_nearest_2_points[j];
		}
		for (size_t j = 0; j < histogram_of_nearest_2_points.size(); j++)
		{
			histogram_of_nearest_2_points[j] /= total;
		}

		histograms_m2.push_back(histogram_of_nearest_2_points);
	}
	std::vector<std::pair<int, int>> final_pairs; 
	
	for (size_t i = 0; i < histograms_m1.size(); i++)
	{
		float min = 99999;
		int min_index = -1; 
		std::vector<float> differences; 
		for (size_t j = 0; j < division_no; j++)
		{
			differences.push_back(0);
		}
		for (size_t j = 0; j < histograms_m2.size(); j++)
		{
			float total = 0;
			for (size_t k = 0; k < division_no; k++)
			{
				total += glm::abs(histograms_m1[i][k] - histograms_m2[j][k]);
			}
			differences.push_back(total);
		}

		//find minimum difference 
		for (size_t j = 0; j < differences.size(); j++)
		{
			if (differences[j] < min)
			{
				min = differences[j];
				min_index = j;
			}
		}

		final_pairs.push_back(std::make_pair(i , min_index));
	}
	glm::vec3 red = glm::vec3(1.0f , 0.0f , 0.0f );
	glm::vec3 green = glm::vec3(0.0f , 1.0f , 0.0f );
	glm::vec3 blue = glm::vec3(0.0f, 0.0f, 1.0f);
	glm::vec3 renk1 = glm::vec3(1.0f, 1.0f, 0.0f);
	glm::vec3 renk2 = glm::vec3(1.0f, 0.0f, 1.0f);
	glm::vec3 renk3 = glm::vec3(0.0f, 1.0f, 1.0f);
	glm::vec3 beyaz = glm::vec3(1.0f, 1.0f, 1.0f);
	glm::vec3 renk4 = glm::vec3(0.0f, 0.25f, 1.0f);
	glm::vec3 renk5 = glm::vec3(0.0f, 1.0f, 0.25f);
	glm::vec3 renk6 = glm::vec3(0.25f, 1.0f, 0.25f);


	std::vector<glm::vec3> m1_colors = m1->colors;
	std::vector<glm::vec3> m2_colors = m2->colors;

	//red 
	int indexm1 = sample_indices_m1[0];
	int indexm2 = sample_indices_m2[0];
	m1_colors[indexm1] = red;
	m2_colors[indexm2] = red;

	indexm1 = sample_indices_m1[1];
	indexm2 = sample_indices_m2[1];
	m1_colors[indexm1] = green;
	m2_colors[indexm2] = green;

	indexm1 = sample_indices_m1[2];
	indexm2 = sample_indices_m2[2];
	m1_colors[indexm1] = blue;
	m2_colors[indexm2] = blue;

	indexm1 = sample_indices_m1[9];
	indexm2 = sample_indices_m2[9];
	m1_colors[indexm1] = renk1;
	m2_colors[indexm2] = renk1;

	indexm1 = sample_indices_m1[3];
	indexm2 = sample_indices_m2[3];
	m1_colors[indexm1] = renk2;
	m2_colors[indexm2] = renk2;

	indexm1 = sample_indices_m1[4];
	indexm2 = sample_indices_m2[4];
	m1_colors[indexm1] = renk3;
	m2_colors[indexm2] = renk3;

	indexm1 = sample_indices_m1[5];
	indexm2 = sample_indices_m2[5];
	m1_colors[indexm1] = renk4;
	m2_colors[indexm2] = renk4;

	indexm1 = sample_indices_m1[6];
	indexm2 = sample_indices_m2[6];
	m1_colors[indexm1] = beyaz;
	m2_colors[indexm2] = beyaz;

	indexm1 = sample_indices_m1[7];
	indexm2 = sample_indices_m2[7];
	m1_colors[indexm1] = renk5;
	m2_colors[indexm2] = renk5;

	indexm1 = sample_indices_m1[8];
	indexm2 = sample_indices_m2[8];
	m1_colors[indexm1] = renk6;
	m2_colors[indexm2] = renk6;

	m1->colors = m1_colors;
	m2->colors = m2_colors;
}

std::vector<float> match_points_from2_mesh_mock(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int division_no )
{
	//get the histogram of the traingle with the most near 2 points
	Mesh* m1 = &mesh_fac.mesh_vec[mesh_index1];
	Mesh* m2 = &mesh_fac.mesh_vec[mesh_index2];

	std::vector<int> sample_indices_m1;
	std::vector<int> sample_indices_m2;
	srand(time(NULL));

	std::vector<float> lines_between_meshes;

	for (size_t i = 0; i < division_no; i++)
	{
		int point_m1 = rand() % m1->vertices.size();
		int point_m2 = rand() % m2->vertices.size();

		sample_indices_m1.push_back(point_m1);
		sample_indices_m2.push_back(point_m2);
	}

	std::vector<float> line_points; 

	//sample from meshes
	for (size_t i = 0; i < division_no; i++)
	{
		//get p1 from m1 
		glm::vec3 p1_v3 = m1->vertices[sample_indices_m1[i]];
		// get p2 from m2
		glm::vec3 p2_v3 = m2->vertices[sample_indices_m2[i]];



		line_points.push_back(p1_v3[0]);
		line_points.push_back(p1_v3[1]);
		line_points.push_back(p1_v3[2]);
		line_points.push_back(-1.0f); //indicator for mesh1 

		line_points.push_back(p2_v3[0]);
		line_points.push_back(p2_v3[1]);
		line_points.push_back(p2_v3[2]);
		line_points.push_back(2.0f); //indicator for mesh2 

	}
	return line_points;
}


void brute_force_symmetry_extraction(MeshFactory& mesh_fac, int& selected_index)
{
	Mesh *m = &mesh_fac.mesh_vec[selected_index];

	float epsilon_length = 1;
	float epsilon_curvature = 1e-2;

	std::vector<int> symmetry_vertex_indices;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		std::cout << "% " <<(float)i / m->vertices.size() << std::endl;
		for (size_t j = 0; j < m->vertices.size(); j++)
		{
			for (size_t k = 0; k < m->vertices.size(); k++)
			{
				if (i == j || i == k || j == k)
				{
					continue;
				}
				//else
				TrilateralDescriptor t1 = generate_trilateral_descriptor(mesh_fac , selected_index , i , j, k , true );
				
				//find the middle vertex to add to symmetry
				if (glm::abs(t1.lenght_1_2 - t1.lenght_1_3) <= epsilon_length)
				{
					//check the curvature index
					if (glm::abs(t1.curvature_1_2 - t1.curvature_1_3) <= epsilon_curvature)
					{
						symmetry_vertex_indices.push_back(i);
					}
				}
				else if (glm::abs(t1.lenght_2_3 - t1.lenght_1_3) <= epsilon_length)
				{
					//check the curvature index
					if (glm::abs(t1.curvature_2_3 - t1.curvature_1_3) <= epsilon_curvature)
					{
						symmetry_vertex_indices.push_back(k);
					}
				}
				else if (glm::abs(t1.lenght_1_2 - t1.lenght_2_3) <= epsilon_length)
				{
					//check the curvature index
					if (glm::abs(t1.curvature_1_2 - t1.curvature_2_3) <= epsilon_curvature)
					{
						symmetry_vertex_indices.push_back(j);
					}
				}
			}
		}
	}
	std::vector<glm::vec3> new_color_buffer; 
	for (size_t i = 0; i < m->colors.size(); i++)
	{
		bool is_index_present = false;
		for (int j = 0; j < symmetry_vertex_indices.size(); j++)
		{
			if (symmetry_vertex_indices[j] == i)
			{
				//color it red
				new_color_buffer.push_back(glm::vec3(1.0f , 0.0f ,0.0f) );
				is_index_present = true; 
				break; 
			}
		}

		if (is_index_present)
		{
			new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
		}
	}
	m->colors = new_color_buffer;
}


//// from the paper Dominant Symmetry Plane Detection for Point-Based 3D Models
//
//Eigen::Matrix<double ,  3 ,3 > symmetry_finding_with_centroid(MeshFactory& mesh_fac, int& selected_index)
//{
//
//	Mesh mesh = mesh_fac.mesh_vec[selected_index]; //mesh itself
//	
//	int mesh_size = mesh.vertices.size();
//	//initialise the wights for the mesh
//	std::vector<double> mesh_weights;
//	for (int i = 0; i < mesh_size; i++)
//	{
//		mesh_weights.push_back(1.0); //assume every vertex is the same 
//	}
//	//double mesh_weights[mesh.vertices.size()]; //initialise weights according to paper 
//	// steps
//	// 1 - get centroid
//	glm::vec3 centroid(0.0f,0.f,0.0f); //initialisiation
//
//	for (size_t i = 0; i < mesh_size; i++)
//	{
//		centroid += mesh.vertices[i];
//	}
//	centroid /= mesh_size; //now we got the centroid 
//	// now convert the centrid glm vec to matrix
//	Eigen::Matrix<double, 3, 1 > centroid_mat;
//	centroid_mat(0) = centroid[0];
//	centroid_mat(1) = centroid[1];
//	centroid_mat(2) = centroid[2];
//	// get the covariance matrix
//	Eigen::Matrix<double, 3, 3> covariance; //covariance mtrix
//	//fill covariance matrix
//	covariance(0, 0) = 0.0;
//	covariance(0, 1) = 0.0;
//	covariance(0, 2) = 0.0;
//	covariance(1, 0) = 0.0;
//	covariance(1, 1) = 0.0;
//	covariance(1, 2) = 0.0;
//	covariance(2, 0) = 0.0;
//	covariance(2, 1) = 0.0;
//	covariance(2, 2) = 0.0;
//
//	for (size_t i = 0; i < mesh_size; i++)
//	{
//		//get P_i and turn it to a matrix 
//		Eigen::Matrix<double, 3, 1 > p_i;
//		p_i(0) = mesh.vertices[i][0];
//		p_i(1) = mesh.vertices[i][1];
//		p_i(2) = mesh.vertices[i][2];
//
//		Eigen::Matrix<double, 3, 1 > p_i_minus_centroid = (p_i - centroid_mat);
//
//		Eigen::Matrix<double, 3, 3 > covariance_i =  mesh_weights[i] * p_i_minus_centroid * p_i_minus_centroid.transpose();
//		
//		covariance += covariance_i;
//	}
//	double s = 0;
//	for (int i = 0; i < mesh_size; i++)
//	{
//		s += mesh_weights[i];
//	}
//	covariance /= s; //covariance is ready for the first part 
//
//	return covariance;
//}
//
//void plane_calculations_from_covariance(Eigen::Matrix<double, 3, 3 >& covariance)
//{
//	//define planes 
//	Eigen::Matrix<double, 4, 1 > p1; 
//	Eigen::Matrix<double, 4, 1 > p2; 
//	Eigen::Matrix<double, 4, 1 > p3;
//
//}