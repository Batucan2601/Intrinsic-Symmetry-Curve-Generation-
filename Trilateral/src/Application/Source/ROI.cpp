#include "../Include/ROI.h"
#include "../Include/Geodesic.h"

static std::vector<int> check_vertices_visited(TrilateralMesh* m, std::vector<int>& path_1_2, std::vector<int>& path_1_3, std::vector<int>& path_2_3);
static std::vector<int> breadth_first_search(TrilateralMesh* m, int point_index, std::vector<int> is_visited);

std::vector<int> ROI_trilateral(TrilateralMesh* m,TrilateralDescriptor& desc, int division_no, bool is_color)
{
	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, desc.p1, desc.p2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, desc.p1, desc.p3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, desc.p2, desc.p3);
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, desc.p1);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, desc.p2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, desc.p3);

	std::vector<int> is_visited = check_vertices_visited(m, path_1_2, path_1_3, path_2_3);

	if (is_color)
	{
		for (size_t i = 0; i < m->colors.size(); i++)
		{
			if (is_visited[i] == OUTSIDE) //edge 
			{
				continue;
			}
			else
			{
				if (is_visited[i] == EDGE) //edge 
				{
					//new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
					m->raylib_mesh.colors[i * 4] = 255;
					m->raylib_mesh.colors[i * 4 + 1] = 0;
					m->raylib_mesh.colors[i * 4 + 2] = 0;
					m->raylib_mesh.colors[i * 4 + 3] = 255;
				}
				else if (is_visited[i] == INSIDE)
				{
					m->raylib_mesh.colors[i * 4] = 0;
					m->raylib_mesh.colors[i * 4 + 1] = 255;
					m->raylib_mesh.colors[i * 4 + 2] = 0;
					m->raylib_mesh.colors[i * 4 + 3] = 255;
				}
				else if (is_visited[i] == MIDPOINT)
				{
					m->raylib_mesh.colors[i * 4] = 255;
					m->raylib_mesh.colors[i * 4 + 1] = 255;
					m->raylib_mesh.colors[i * 4 + 2] = 255;
					m->raylib_mesh.colors[i * 4 + 3] = 255;
				}
				desc.visited_indices.push_back(i);
			}
			
		}
	}

	return  is_visited;
}


static std::vector<int> check_vertices_visited(TrilateralMesh* m, std::vector<int>& path_1_2, std::vector<int>& path_1_3, std::vector<int>& path_2_3)
{
	std::vector<int> start_vertices;
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		start_vertices.push_back(path_1_2[i]);
	}
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		start_vertices.push_back(path_1_3[i]);
	}
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		start_vertices.push_back(path_2_3[i]);
	}

	//make the list unique
	/*std::sort(start_vertices.begin(), start_vertices.end());
	start_vertices.erase(std::unique(start_vertices.begin(), start_vertices.end()), start_vertices.end());
	start_vertices.push_back(path_1_2[0]);
	start_vertices.push_back(path_1_3[path_1_3.size() - 1]);
	start_vertices.push_back(path_2_3[0]);*/
	std::vector<int> is_visited;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		is_visited.push_back(OUTSIDE);
	}
	std::vector<int> start_vertices_neighbours;

	//set the vertices
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		is_visited[i] = OUTSIDE;
	}
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		is_visited[path_1_2[i]] = EDGE;
	}
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		is_visited[path_1_3[i]] = EDGE;

	}
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		is_visited[path_2_3[i]] = EDGE;
	}
	std::vector<std::vector<int>> visited_vertices_list;
	// get the neighbours
	for (size_t i = 0; i < start_vertices.size(); i++)
	{
		for (size_t j = 0; j < /*m->adjacenies[start_vertices[i]].size() */ 1; j++) //only check 1 neighbour this should emirically work ?
		{
			int point_index = m->adjacenies[start_vertices[i]][j].first;
			if (is_visited[point_index] != EDGE)
			{
				start_vertices_neighbours.push_back(point_index);
				std::vector<int> visited_points = breadth_first_search(m, point_index, is_visited);
				visited_vertices_list.push_back(visited_points);
			}
		}
	}
	int minimum_val = m->vertices.size(); // maximum an be mesh size 
	int minimum_index = -1;
	for (size_t i = 0; i < visited_vertices_list.size(); i++)
	{
		if (minimum_val > visited_vertices_list[i].size())
		{
			minimum_val = visited_vertices_list[i].size();
			minimum_index = i;
		}
	}
	for (size_t i = 0; i < visited_vertices_list[minimum_index].size(); i++)
	{
		int inside_index = visited_vertices_list[minimum_index][i];
		is_visited[inside_index] = INSIDE;
	}
	std::vector<int> return_visited_values;

	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		return_visited_values.push_back(is_visited[i]);
	}

	//debug ! 
	//there shouldo nly be 2 unique values
	std::vector<int> unique;
	for (size_t i = 0; i < visited_vertices_list.size(); i++)
	{
		bool is_exists = false;
		for (size_t j = 0; j < unique.size(); j++)
		{
			if (unique[j] == visited_vertices_list[i].size())
			{
				is_exists = true;
				break;
			}
		}
		if (!is_exists)
		{
			unique.push_back(visited_vertices_list[i].size());
		}
	}
	if (unique.size() > 2)
	{
		int debug = 1;
	}
	return return_visited_values;

}

//given edges and point do breadth first search on an area 
static std::vector<int> breadth_first_search(TrilateralMesh* m, int point_index, std::vector<int> is_visited)
{
	std::vector<int> visited_vertices;
	//push our point to stack
	std::stack<int> stack;
	stack.push(point_index);
	while (!stack.empty())
	{
		int index = stack.top();
		stack.pop(); //vertex index popped from stack
		if (is_visited[index] == OUTSIDE) //not visited
		{
			is_visited[index] = INSIDE; // now te vertex has been visited
			visited_vertices.push_back(index);
			// this region of loop assumes index is not edge, therefore add the adjacencies
			for (size_t i = 0; i < m->adjacenies[index].size(); i++) //process pairs 
			{
				stack.push(m->adjacenies[index][i].first);
			}
		}
		if (is_visited[index] == EDGE) //do nothing 
		{
			;
		}
	}
	return visited_vertices;
}