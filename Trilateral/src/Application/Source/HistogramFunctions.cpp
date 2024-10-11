#include "../Include/HistogramFunctions.h"
#include "../Include/Geodesic.h"
#define _USE_MATH_DEFINES
#include "math.h"
Histogram histogram_roi_area_detailed(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited)
{
	//histogram to be returned 
	Histogram histogram(division_no);
	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, point_index1);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, point_index2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, point_index3);
	//find the maximum distance from is_visited and paths
	float max = -999;
	float min = 100000;
	std::vector<unsigned int> triangles_inside;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (is_visited[i] == INSIDE) // if vertex is visited 
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
	//now recolor
	std::vector<glm::vec3> new_color_buffer;
	for (size_t i = 0; i < m->colors.size(); i++)
	{

		if (is_visited[i] == EDGE) //edge 
		{
			global_is_visited[i] = EDGE;
		}
		else if (is_visited[i] == OUTSIDE) //not visited 
		{
		}
		else if (is_visited[i] == INSIDE) // get the max distance
		{
			if (global_is_visited[i] != EDGE)
			{
				global_is_visited[i] = INSIDE;
			}
		}
	}
	float max_dist_inside = -1;
	float min_dist_inside = -1;
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		if (is_visited[index1] == INSIDE || is_visited[index2] == INSIDE || is_visited[index3] == INSIDE) //if any vertex is visited
		{
			triangles_inside.push_back(index1);
			triangles_inside.push_back(index2);
			triangles_inside.push_back(index3);
			if (distance_matrix_p1[index1] > max_dist_inside)
			{
				max_dist_inside = distance_matrix_p1[index1];
			}
			if (distance_matrix_p1[index2] > max_dist_inside)
			{
				max_dist_inside = distance_matrix_p1[index2];
			}
			if (distance_matrix_p1[index3] > max_dist_inside)
			{
				max_dist_inside = distance_matrix_p1[index3];
			}

			if (distance_matrix_p1[index1] < min_dist_inside)
			{
				min_dist_inside = distance_matrix_p1[index1];
			}
			if (distance_matrix_p1[index2] < min_dist_inside)
			{
				min_dist_inside = distance_matrix_p1[index2];
			}
			if (distance_matrix_p1[index3] < min_dist_inside)
			{
				min_dist_inside = distance_matrix_p1[index3];
			}
		}
	}
	float step = max_dist_inside / division_no;
	for (size_t i = 0; i < triangles_inside.size(); i += 3)
	{
		int i1 = triangles_inside[i];
		int i2 = triangles_inside[i + 1];
		int i3 = triangles_inside[i + 2];
		glm::vec3 p1 = m->vertices[i1];
		glm::vec3 p2 = m->vertices[i2];
		glm::vec3 p3 = m->vertices[i3];
		float i1_dist = distance_matrix_p1[i1];
		float i2_dist = distance_matrix_p1[i2];
		float i3_dist = distance_matrix_p1[i3];

		int step_no_i1 = i1_dist / step; // floor 
		int step_no_i2 = i2_dist / step; // floor 
		int step_no_i3 = i3_dist / step; // floor
		if (step_no_i1 == division_no)
		{
			step_no_i1--;
		}
		if (step_no_i2 == division_no)
		{
			step_no_i2--;
		}
		if (step_no_i3 == division_no)
		{
			step_no_i3--;
		}
		float triangle_area = compute_triangle_area(p1, p2, p3);
		if (step_no_i1 == step_no_i2 && step_no_i1 == step_no_i3)
		{
			histogram[step_no_i1] += triangle_area;
		}
		//one of them is in other step 
		else /*if (step_no_i1 == step_no_i2 && step_no_i1 != step_no_i3 &&
			 step_no_i1 == step_no_i3 && step_no_i1 != step_no_i2 &&
			 step_no_i2 == step_no_i3 && step_no_i2 != step_no_i3)*/
		{
			std::vector<std::pair<int, int>> steps;
			steps.push_back(std::pair<int, int>(step_no_i1, i1));
			steps.push_back(std::pair<int, int>(step_no_i2, i2));
			steps.push_back(std::pair<int, int>(step_no_i3, i3));
			CoreType_sort_by_value(steps);
			//they are now in ascending order
			// 1.1.1 small arc
			float dist_step_0 = distance_matrix_p1[steps[0].second];
			//find the distance to other hist
			float small_r = (step * (steps[0].first + 1)) - dist_step_0;
			glm::vec3 edge1 = m->vertices[steps[1].second] - m->vertices[steps[0].second];
			glm::vec3 edge2 = m->vertices[steps[2].second] - m->vertices[steps[0].second];
			float cosine = glm::dot(edge1, edge2) / (glm::length(edge1) * glm::length(edge2));
			float small_radian = acos(cosine); //radian

			float area_arc_small = M_PI * small_r * small_r * small_radian / (2 * M_PI);

			histogram[steps[0].first] += area_arc_small;


			// 1.1.2 closest arc
			//a miscalculation but think it like also an arc although it is a reverse arc
			float dist_step_2 = distance_matrix_p1[steps[2].second];
			float big_r = dist_step_2 - (steps[2].first * step);
			edge1 = m->vertices[steps[1].second] - m->vertices[steps[2].second];
			edge2 = m->vertices[steps[0].second] - m->vertices[steps[2].second];
			cosine = glm::dot(edge1, edge2) / (glm::length(edge1) * glm::length(edge2));
			small_radian = acos(cosine); //radian
			float area_arc_big = M_PI * big_r * big_r * small_radian / (2 * M_PI);
			histogram[steps[2].first] += area_arc_big;


			float area_left = triangle_area - (area_arc_small + area_arc_big);

			int steps_left = steps[2].first - (steps[0].first + 1);
			if (steps_left == 0)
			{
				histogram[steps[0].first] += area_left;
				histogram[steps[1].first] += area_left;
			}
			else
			{
				int total_fraction = steps_left * steps_left;
				for (size_t step_no = steps[0].first + 1; step_no < steps[2].first; step_no++)
				{
					int fraction_no = step_no - (steps[0].first);
					float fraction_percentage = (pow(fraction_no, 2) - pow(fraction_no - 1, 2)) / fraction_no;

					histogram[step_no] += fraction_percentage * area_left;
				}
			}



		}

	}

	histogram.normalize(1);

	return histogram;

}


Histogram  Histogram_triangle_area(TrilateralMesh* m, TrilateralDescriptor& desc, int division_no, std::vector<int>& global_is_visited)
{
	Histogram histogram(division_no);
	std::vector<int> path_1_2 = desc.path_1_2;
	std::vector<int> path_1_3 = desc.path_1_3;
	std::vector<int> path_2_3 = desc.path_2_3;
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, desc.p1);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, desc.p2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, desc.p3);
	//find the maximum distance from is_visited and paths
	float max = -999;
	float min = 100000;
	std::vector<unsigned int> is_visited = desc.visited_indices;
	if (is_visited.size() == 0)
	{
		return Histogram(division_no);
	}

	for (size_t i = 0; i < is_visited.size(); i++)
	{
		float len = distance_matrix_p1[is_visited[i]];
		if (len > max)
		{
			max = len; 
		}
	}

	float step = max / division_no;

	
	//m->colors = new_color_buffer;
	std::vector<unsigned int>all_vertex_status(m->vertices.size(), OUTSIDE);
	for (size_t i = 0; i < is_visited.size(); i++)
	{
		all_vertex_status[is_visited[i]] = INSIDE; 
	}
	//getting the areas
	for (size_t i = 0; i < m->triangles.size(); i+=3)
	{
		if (all_vertex_status[m->triangles[i]] != OUTSIDE || all_vertex_status[m->triangles[i + 1]] != OUTSIDE || all_vertex_status[m->triangles[i + 2]] != OUTSIDE) //if any vertex is visited
		{
			glm::vec3 p1 = m->vertices[m->triangles[i]];
			glm::vec3 p2 = m->vertices[m->triangles[i + 1]];
			glm::vec3 p3 = m->vertices[m->triangles[i + 2]];

			float distp1 = distance_matrix_p1[m->triangles[i]];
			float distp2 = distance_matrix_p1[m->triangles[i + 1]];
			float distp3 = distance_matrix_p1[m->triangles[i + 2]];

			float area = compute_triangle_area(p1, p2, p3);

			int hist_no_p1 = distp1 / step;
			int hist_no_p2 = distp2 / step;
			int hist_no_p3 = distp3 / step;

			if (hist_no_p1 >= division_no)
			{
				hist_no_p1 = division_no - 1;
			}
			if (hist_no_p2 >= division_no)
			{
				hist_no_p2 = division_no - 1;
			}
			if (hist_no_p3 >= division_no)
			{
				hist_no_p3 = division_no - 1;
			}
			histogram[hist_no_p1] += area;
			histogram[hist_no_p2] += area;
			histogram[hist_no_p3] += area;
		}
	}

	//normalize histogram.
	float histogram_sum = 0;
	for (size_t i = 0; i < histogram.size(); i++)
	{
		histogram_sum += histogram[i];
	}
	for (size_t i = 0; i < histogram.size(); i++)
	{
		histogram[i] /= histogram_sum;
	}
	return histogram;
}

//use m_inner 
Histogram Histogram_triangle_area_w_res(TrilateralMesh* m, TrilateralDescriptor& desc, int division_no , int resolution)
{
	Histogram histogram(division_no);
	std::vector<int> path_1_2 = desc.path_1_2;
	std::vector<int> path_1_3 = desc.path_1_3;
	std::vector<int> path_2_3 = desc.path_2_3;

	int p1_index = desc.p1;
	glm::vec3 p1_point = m->vertices[p1_index];
	//generate m_inner

	int p1_new_index = -1;
	TrilateralDescriptor_generate_mesh_with_resolution(m ,desc , resolution);
	TrilateralMesh* m_inside = &desc.m_inside;
	for (size_t i = 0; i < desc.m_inside.vertices.size(); i++)
	{
		if (glm::length(desc.m_inside.vertices[i] - p1_point) < 1e-7)
		{
			p1_new_index = i;
		}
	}
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(desc.m_inside, p1_new_index);

	float max = -INFINITY; 
	for (size_t i = 0; i < m_inside->vertices.size(); i++)
	{
		float len = distance_matrix_p1[i];
		if (len > max && i != p1_new_index)
		{
			max = len;
		}
	}

	float step = max / division_no;


	
	for (size_t i = 0; i < m_inside->triangles.size(); i += 3)
	{
		glm::vec3 p1 = m_inside->vertices[m_inside->triangles[i]];
		glm::vec3 p2 = m_inside->vertices[m_inside->triangles[i + 1]];
		glm::vec3 p3 = m_inside->vertices[m_inside->triangles[i + 2]];

		float distp1 = distance_matrix_p1[m_inside->triangles[i]];
		float distp2 = distance_matrix_p1[m_inside->triangles[i + 1]];
		float distp3 = distance_matrix_p1[m_inside->triangles[i + 2]];

		float area = compute_triangle_area(p1, p2, p3);

		int hist_no_p1 = distp1 / step;
		int hist_no_p2 = distp2 / step;
		int hist_no_p3 = distp3 / step;

		if (hist_no_p1 >= division_no)
		{
			hist_no_p1 = division_no - 1;
		}
		if (hist_no_p2 >= division_no)
		{
			hist_no_p2 = division_no - 1;
		}
		if (hist_no_p3 >= division_no)
		{
			hist_no_p3 = division_no - 1;
		}
		histogram[hist_no_p1] += area;
		histogram[hist_no_p2] += area;
		histogram[hist_no_p3] += area;
	}

	//normalize histogram.
	float histogram_sum = 0;
	for (size_t i = 0; i < histogram.size(); i++)
	{
		histogram_sum += histogram[i];
	}
	for (size_t i = 0; i < histogram.size(); i++)
	{
		histogram[i] /= histogram_sum;
	}
	return histogram;
}