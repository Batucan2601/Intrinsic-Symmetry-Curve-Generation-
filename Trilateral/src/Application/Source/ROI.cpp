#include "../Include/ROI.h"
#include "../Include/Geodesic.h"

static std::vector<unsigned int> check_vertices_visited(TrilateralMesh* m, std::vector<int>& path_1_2, std::vector<int>& path_1_3, std::vector<int>& path_2_3);
static std::vector<unsigned int> breadth_first_search(TrilateralMesh* m, int point_index, std::vector<int> is_visited);

void ROI_trilateral(TrilateralMesh* m,TrilateralDescriptor& desc, int division_no, bool is_color)
{
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, desc.p1);
	if (distance_matrix_p1[desc.p3] < distance_matrix_p1[desc.p2])
	{
		int temp = desc.p3;
		desc.p3 = desc.p2;
		desc.p2 = temp;
	}
	
	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, desc.p1, desc.p2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, desc.p1, desc.p3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, desc.p2, desc.p3);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, desc.p2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, desc.p3);
	desc.path_1_2 = path_1_2;
	desc.path_1_3 = path_1_3;
	desc.path_2_3 = path_2_3;

	desc.geodesic_lenght_1_2 = distance_matrix_p1[desc.p2];
	desc.geodesic_lenght_1_3 = distance_matrix_p1[desc.p3];
	desc.geodesic_lenght_2_3 = distance_matrix_p2[desc.p3];
	std::vector<unsigned int> visited_indices = check_vertices_visited(m, path_1_2, path_1_3, path_2_3);
	desc.visited_indices = visited_indices;

	//generate histograms with the areas
	std::vector<float> path_1_2_area;
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		int index = path_1_2[i];
		path_1_2_area.push_back(m->areas[index]);
	}
	desc.hist_path_1_2.histogram = path_1_2_area;
	desc.hist_path_1_2.normalize(1);
	std::vector<float> path_1_3_area;
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		int index = path_1_3[i];
		path_1_3_area.push_back(m->areas[index]);
	}
	desc.hist_path_1_3.histogram = path_1_3_area;
	desc.hist_path_1_3.normalize(1);

	std::vector<float> path_2_3_area;
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		int index = path_2_3[i];
		path_2_3_area.push_back(m->areas[index]);
	}
	desc.hist_path_2_3.histogram = path_2_3_area;
	desc.hist_path_2_3.normalize(1);

}


static std::vector<unsigned int> check_vertices_visited(TrilateralMesh* m, std::vector<int>& path_1_2, std::vector<int>& path_1_3, std::vector<int>& path_2_3)
{
	std::vector<int> edge_vertices;
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		edge_vertices.push_back(path_1_2[i]);
	}
	for (size_t i = 0; i < path_1_3.size(); i++)
	{
		edge_vertices.push_back(path_1_3[i]);
	}
	for (size_t i = 0; i < path_2_3.size(); i++)
	{
		edge_vertices.push_back(path_2_3[i]);
	}

	std::vector<int> is_visited;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		is_visited.push_back(OUTSIDE);
	}
	std::vector<int> edge_vertices_neighbours;

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
	//check if colinear 2D
	std::vector<std::vector<unsigned int>> visited_vertices_list;
	// get the neighbours
	for (size_t i = 0; i < edge_vertices.size(); i++)
	{
		for (size_t j = 0; j < 1; j++) //only check 1 neighbour this should emirically work ?
		{
			int point_index = m->adjacenies[edge_vertices[i]][j].first;
			if (is_visited[point_index] != EDGE)
			{
				edge_vertices_neighbours.push_back(point_index);
				std::vector<unsigned int> visited_points = breadth_first_search(m, point_index, is_visited);
				visited_vertices_list.push_back(visited_points);
			}
		}
	}

	
	//the minimum sized batch is our inside 
	int minimum_size = INFINITY; // maximum an be mesh size 
	int minimum_index = -1;
	for (size_t i = 0; i < visited_vertices_list.size(); i++)
	{
		if (minimum_size > visited_vertices_list[i].size())
		{
			minimum_size = visited_vertices_list[i].size();
			minimum_index = i;
		}
	}

	// colinear no Inside 
	if (minimum_size > 95.0 / 100.0 * m->vertices.size())
	{
		std::vector<unsigned int> empty_list;
		return empty_list;
	}


	// if every index has same length they are colinear
	return visited_vertices_list[minimum_index];
}

//given edges and point do breadth first search on an area 
static std::vector<unsigned int> breadth_first_search(TrilateralMesh* m, int point_index, std::vector<int> is_visited)
{
	std::vector<unsigned int> visited_vertices;
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