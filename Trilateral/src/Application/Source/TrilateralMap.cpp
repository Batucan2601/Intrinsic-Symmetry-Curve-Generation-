#include "../Include/TrilateralMap.h"
#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/MetricCalculations.h"
#include "../Include/Geodesic.h"
#include "../Include/Skeleton.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/ROI.h"
#include "../Include/Histogram.h"
#include "../Include/HistogramFunctions.h"
#include "../Include/SpinImage.h"
#include "../Include/VarianceMinimizingTransportPlan.h"
#include "../Include/NLateralDescriptor.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "../Include/HeatKernelSignature.h"


static std::vector<float> computePrincipalCurvatures(MeshFactory& meshFac, int selectedIndex, std::vector<int>& is_visited,
	std::vector<double>& principalCurvatures1, std::vector<double>& principalCurvatures2, int partition_no);
static std::vector<int> getPointsInside(MeshFactory& meshFac, int selectedIndex, std::vector<int>& is_visited);
static Eigen::Matrix3d computeCurvatureTensor(const TrilateralMesh& mesh, int vertexIndex);
static std::vector<std::pair<unsigned int, unsigned int>>  trilateral_unique_pairing_histogram(std::vector<TrilateralDescriptor>& trilateral_desc_left,
	std::vector<TrilateralDescriptor>& trilateral_desc_right);
static std::vector<std::pair<unsigned int, unsigned int>>  trilateral_unique_pairing_histogram_w_HKS(TrilateralMesh* m , std::vector<TrilateralDescriptor>& trilateral_desc_left,
	std::vector<TrilateralDescriptor>& trilateral_desc_right);
static std::vector<std::pair<unsigned int, unsigned int>>  trilateral_area_histogram(TrilateralMesh* m, std::vector<TrilateralDescriptor>& trilateral_desc_left,
	std::vector<TrilateralDescriptor>& trilateral_desc_right);
std::vector<float> compute_geodesic_distances_min_heap_distances(TrilateralMesh& m, int point_index)
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

void trilateral_map_drawing_using_three_points(MeshFactory& mesh_fac, int& selected_index, int p1, int p2, int p3)
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];
	//extract two paths
	//1 - extract the path from p1  to p2  and p1 to p3  and p2 to p3 
	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, p1, p2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, p1, p3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, p2, p3);

#pragma region calculation of distances between points 2 ,3 and 1  

	float dist_1_2 = 0.0f;
	float dist_1_3 = 0.0f;
	float dist_2_3 = 0.0f;

	for (size_t i = 0; i < path_1_2.size() - 1; i++)
	{
		dist_1_2 += glm::distance(m->vertices[path_1_2[i]], m->vertices[path_1_2[i + 1]]);
	}
	for (size_t i = 0; i < path_1_3.size() - 1; i++)
	{
		dist_1_3 += glm::distance(m->vertices[path_1_3[i]], m->vertices[path_1_3[i + 1]]);
	}
	for (size_t i = 0; i < path_2_3.size() - 1; i++)
	{
		dist_2_3 += glm::distance(m->vertices[path_2_3[i]], m->vertices[path_2_3[i + 1]]);
	}
	std::cout << "dist1 " << dist_1_2 << std::endl;
	std::cout << "dist 2 " << dist_1_3 << std::endl;
#pragma endregion 

#pragma region  angle
	// for p1
	float angle_p1 = 0.0f; // in degrees ( from angle betwen 2 , 1 ,3 )
	glm::vec3 p1_point = glm::vec3(m->vertices[p1]);
	glm::vec3 p2_point = glm::vec3(m->vertices[p2]);
	glm::vec3 p3_point = glm::vec3(m->vertices[p3]);
	glm::vec3 vec_2_1 = glm::vec3(p2_point - p1_point);
	glm::vec3 vec_3_1 = glm::vec3(p3_point - p1_point);
	double angle_radian_p1 = glm::acos(glm::dot(vec_2_1, vec_3_1) / (glm::length(vec_2_1) * glm::length(vec_3_1)));
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

	std::vector<glm::vec3> draw_buffer;
	std::vector<glm::vec3> draw_buffer_color;
	
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		glm::vec3 point = m->vertices[i];
		bool is_in_way = false;
		for (size_t j = 0; j < path_1_2.size(); j++)
		{
			glm::vec3 point_between_1_2 = m->vertices[path_1_2[j]];

			if (point.x == point_between_1_2.x && point.y == point_between_1_2.y && point.z == point_between_1_2.z)
			{
				is_in_way = true;
			}
		}
		for (size_t j = 0; j < path_1_3.size(); j++)
		{
			glm::vec3 point_between_1_3 = m->vertices[path_1_3[j]];

			if (point.x == point_between_1_3.x && point.y == point_between_1_3.y && point.z == point_between_1_3.z)
			{
				is_in_way = true;
			}
		}
		for (size_t j = 0; j < path_2_3.size(); j++)
		{
			glm::vec3 point_between_2_3 = m->vertices[path_2_3[j]];

			if (point.x == point_between_2_3.x && point.y == point_between_2_3.y && point.z == point_between_2_3.z)
			{
				is_in_way = true;
			}
		}
		// rebuffer data
		//draw_buffer.push_back(point.x);
		//draw_buffer.push_back(point.y);
		//draw_buffer.push_back(point.z);

		glm::vec3 color;
		//red 
		if (is_in_way)
		{
			if (i == p1 || i == p2 || i == p3) //point itself  (will draw red therefore skip )
			{
				//draw_buffer.push_back(0.0f);
				color.r = 0.0f; 
			}
			else
			{
				//draw_buffer.push_back(1.0f); // on the way 
				color.r = 1.0f;

			}


		}
		else
		{
			//draw_buffer.push_back(0.0f);
			color.r = 0.0f;
		}

		//green 
		//draw_buffer.push_back(0.0f);
		color.g = 0.0f;


		//blue
		if (i == p1 || i == p2 || i == p3) //point itself  
		{
			//draw_buffer.push_back(1.0f);
			color.b = 1.0f;
		}
		else
		{
			//draw_buffer.push_back(0.0f);
			color.b = 0.0f;
		}

		m->colors[i] = color; 

	}
	int point_size = 0;
	for (size_t i = 0; i < selected_index; i++)
	{
		point_size += mesh_fac.mesh_vec[i].vertices.size() * 6;
	}
	
	//glBufferSubData(GL_ARRAY_BUFFER, point_size * sizeof(float), draw_buffer.size() * sizeof(float), &draw_buffer[0]);
	/*glBufferData(GL_ARRAY_BUFFER, draw_buffer.size() * sizeof(float), &draw_buffer[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.triangles.size() * sizeof(int), &m.triangles[0], GL_STATIC_DRAW);*/

}


/*static float* trilateral_ROI(MeshFactory& mesh_fac  , int& selected_index, int point_index1, int point_index2, int point_index3, int division_no) // tau is the closeness division_no is the no of how much you want to separate
{
	TrilateralMesh m = mesh_fac.mesh_vec[selected_index];
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
	std::vector<int> path_1_2 = Geodesic_between_two_points(m, point_index1, point_index2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(m, point_index1, point_index3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(m, point_index2, point_index3);

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
		minimum_distance[i] = INFINITY;

	}
	//calculate the distance from every vertex within path
	for (size_t i = 0; i < path_1_2.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(m, path_1_2[i]);
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
		std::vector<float> distances = Geodesic_dijkstra(m, path_2_3[i]);
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
		std::vector<float> distances = Geodesic_dijkstra(m, path_1_3[i]);
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
 std::vector<TrilateralDescriptor> get_trilateral_points_using_closest_pairs(MeshFactory& mesh_fac, const int& selected_index, std::vector<unsigned int>& indices)
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];

	std::vector<TrilateralDescriptor> trilateralDesc;
	for (size_t i = 0; i < indices.size(); i++)
	{
		TrilateralDescriptor desc;
		//get two of the closed indexed points
		std::vector<float> geodesic_distances = Geodesic_dijkstra(*m, indices[i]);
		std::vector<std::pair<float, unsigned int >> distances;
		for (size_t j = 0; j < geodesic_distances.size(); j++)
		{
			bool is_in_indices = false;
			for (size_t k = 0; k < indices.size(); k++)
			{
				{
					is_in_indices = true;
					break;
				}
			}

			if (is_in_indices)
			{
				float dist = geodesic_distances[j];
				unsigned int index = j;
				std::pair<float, unsigned int > pair;
				pair.first = dist;
				pair.second = j;
				distances.push_back(pair);
			}
		}
		//get first and second closest
		float minVal = INFINITY;
		float minIndexFirst = -1;
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (minVal > distances[j].first)
			{
				minIndexFirst = distances[j].second;
				minVal = distances[j].first;
			}
		}

		minVal = INFINITY;
		float minIndexSecond = -1;
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (minVal > distances[j].first && minIndexFirst != distances[j].second && distances[j].second != indices[i])
			{
				minIndexSecond = distances[j].second;
				minVal = distances[j].first;
			}
		}
		desc.p1 = indices[i];
		desc.p2 = minIndexFirst;
		desc.p3 = minIndexSecond;

		int mesh_no = selected_index; //casting problems 

		desc = generate_trilateral_descriptor(mesh_fac, mesh_no, desc.p1, desc.p2, desc.p3, true); // do not compute area for now
		trilateralDesc.push_back(desc);
	}

	return trilateralDesc;
}
 std::vector<TrilateralDescriptor> get_trilateral_points_using_closest_pairs(TrilateralMesh* m, std::vector<unsigned int>& indices)
 {
	 std::vector<TrilateralDescriptor> trilateralDesc;
	 for (size_t i = 0; i < indices.size(); i++)
	 {
		 TrilateralDescriptor desc;
		 //get two of the closed indexed points
		 std::vector<float> geodesic_distances = Geodesic_dijkstra(*m, indices[i]);
		 std::vector<std::pair<float, unsigned int >> distances;
		 desc.p1 = indices[i];
		 float smallestVal = INFINITY; 
		 int smallestIndex = -1;
		 for (size_t j = 0; j < indices.size(); j++)
		 {
			 if (j == i)
			 {
				 continue; 
			 }
			 if (geodesic_distances[indices[j]] < smallestVal)
			 {
				 smallestVal = geodesic_distances[indices[j]];
				 smallestIndex = indices[j];
			 }
		 }
		 desc.p2 = smallestIndex;
		 geodesic_distances[smallestIndex] = INFINITY; 
		 smallestVal = INFINITY; 
		 for (size_t j = 0; j < indices.size(); j++)
		 {
			 if (j == i)
			 {
				 continue;
			 }
			 if (geodesic_distances[indices[j]] < smallestVal)
			 {
				 smallestVal = geodesic_distances[indices[j]];
				 smallestIndex = indices[j];
			 }
		 }
		 desc.p3 = smallestIndex;
		 desc = TrilateralDescriptor_create(m, desc.p1, desc.p2, desc.p3, true); // do not compute area for now
		 trilateralDesc.push_back(desc);
	 }
	 return trilateralDesc;
 }
 std::vector<TrilateralDescriptor> get_trilateral_points_using_furthest_pairs(TrilateralMesh* m, std::vector<unsigned int>& indices)
 {
	 std::vector<TrilateralDescriptor> trilateralDesc;
	 for (size_t i = 0; i < indices.size(); i++)
	 {
		 TrilateralDescriptor desc;
		 //get two of the closed indexed points
		 std::vector<float> geodesic_distances = Geodesic_dijkstra(*m, indices[i]);
		 std::vector<std::pair<float, unsigned int >> distances;
		 for (size_t j = 0; j < geodesic_distances.size(); j++)
		 {
			 bool is_in_indices = false;
			 for (size_t k = 0; k < indices.size(); k++)
			 {
				 if (j == indices[k] && j != indices[i])
				 {
					 is_in_indices = true;
					 break;
				 }
			 }

			 if (is_in_indices)
			 {
				 float dist = geodesic_distances[j];
				 unsigned int index = j;
				 std::pair<float, unsigned int > pair;
				 pair.first = dist;
				 pair.second = j;
				 distances.push_back(pair);
			 }
		 }
		 //get first and second farthest
		 float maxVal = -INFINITY;
		 float maxIndexFirst = -1;
		 for (size_t j = 0; j < distances.size(); j++)
		 {
			 if (maxVal < distances[j].first)
			 {
				 maxIndexFirst = distances[j].second;
				 maxVal = distances[j].first;
			 }
		 }

		 maxVal = -INFINITY;
		 float maxIndexSecond = -1;
		 for (size_t j = 0; j < distances.size(); j++)
		 {
			 if (maxVal < distances[j].first && maxIndexFirst != distances[j].second && distances[j].second != indices[i])
			 {
				 maxIndexSecond = distances[j].second;
				 maxVal = distances[j].first;
			 }
		 }
		 desc.p1 = indices[i];
		 desc.p2 = maxIndexFirst;
		 desc.p3 = maxIndexSecond;


		 desc = TrilateralDescriptor_create(m, desc.p1, desc.p2, desc.p3, true); // do not compute area for now
		 trilateralDesc.push_back(desc);
	 }
	 return trilateralDesc;
 }
 

 static int getClosestMeshIndex(TrilateralMesh* m,int point_index1, int point_index2, int point_index3)
 {
	 glm::vec3 median = (m->vertices[point_index1] + m->vertices[point_index2] + m->vertices[point_index3]) / 3.0f;
	 float closest_dist = INFINITY;
	 int closest_index = -1; 
	 for (size_t i = 0; i < m->vertices.size(); i++)
	 {
		 float len = glm::distance(m->vertices[i], median);
		 if (len < closest_dist)
		 {
			 closest_dist = len;
			 closest_index = i;
		 }
	 }
	 return closest_index;
 }





std::vector<float> histogramROi_w_HKS(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited)
{
	//histogram to be returned 
	std::vector<float> histogram;
	// fill it with division no 
	for (size_t i = 0; i < division_no; i++)
	{
		histogram.push_back(0);
	}
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];
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
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i+1];
		int index3 = m->triangles[i+2];
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
		float triangle_area= compute_triangle_area(p1, p2, p3);
		if (step_no_i1 == step_no_i2 && step_no_i1 == step_no_i3)
		{
			histogram[step_no_i1] += triangle_area;
		}
		//one of them is in other step 
		else /*if (step_no_i1 == step_no_i2 && step_no_i1 != step_no_i3 &&
			 step_no_i1 == step_no_i3 && step_no_i1 != step_no_i2 &&
			 step_no_i2 == step_no_i3 && step_no_i2 != step_no_i3)*/
		{
			histogram[step_no_i1] += triangle_area/3.0f;
			histogram[step_no_i2] += triangle_area/3.0f;
			histogram[step_no_i3] += triangle_area/3.0f;

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

// the ultimate trilateral descriptor generator 
 TrilateralDescriptor  generate_trilateral_descriptor(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, bool is_simplified)
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];

	TrilateralDescriptor trilateral_descriptor;//trialteral descriptor 
	//init descriptor
	trilateral_descriptor.area = 0;
	trilateral_descriptor.geodesic_lenght_1_2 = 0;
	trilateral_descriptor.geodesic_lenght_1_3 = 0;
	trilateral_descriptor.geodesic_lenght_2_3 = 0;


	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, point_index1);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, point_index2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, point_index3);

	// get distances
	trilateral_descriptor.geodesic_lenght_1_2 = distance_matrix_p1[point_index2];
	trilateral_descriptor.geodesic_lenght_1_3 = distance_matrix_p1[point_index3];
	trilateral_descriptor.geodesic_lenght_2_3 = distance_matrix_p2[point_index3];

	trilateral_descriptor.curvature_1_2 = trilateral_descriptor.geodesic_lenght_1_2 / glm::distance(m->vertices[point_index1], m->vertices[point_index2]);
	trilateral_descriptor.curvature_1_3 = trilateral_descriptor.geodesic_lenght_1_3 / glm::distance(m->vertices[point_index1], m->vertices[point_index3]);
	trilateral_descriptor.curvature_2_3 = trilateral_descriptor.geodesic_lenght_2_3 / glm::distance(m->vertices[point_index2], m->vertices[point_index3]);

	trilateral_descriptor.euclidian_lenght_1_2 = glm::distance(m->vertices[point_index1], m->vertices[point_index2]);
	trilateral_descriptor.euclidian_lenght_1_3 = glm::distance(m->vertices[point_index1], m->vertices[point_index3]);
	trilateral_descriptor.euclidian_lenght_2_3 = glm::distance(m->vertices[point_index2], m->vertices[point_index3]);

	trilateral_descriptor.curvature_1_2 = trilateral_descriptor.geodesic_lenght_1_2 / trilateral_descriptor.euclidian_lenght_1_2;
	trilateral_descriptor.curvature_1_3 = trilateral_descriptor.geodesic_lenght_1_3 / trilateral_descriptor.euclidian_lenght_1_3;
	trilateral_descriptor.curvature_2_3 = trilateral_descriptor.geodesic_lenght_2_3 / trilateral_descriptor.euclidian_lenght_2_3;

	trilateral_descriptor.p1 = point_index1;
	trilateral_descriptor.p2 = point_index2;
	trilateral_descriptor.p3 = point_index3;


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
					std::vector<int> path_i_j = Geodesic_between_two_points(*m, i, j);

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

			float area = compute_triangle_area(p1, p2, p3);

			trilateral_descriptor.area += area; // get area


		}
	}


	return trilateral_descriptor;
}

 std::vector<std::pair<unsigned int, unsigned int>>  point_match_trilateral_weights(TrilateralMesh*mesh, std::vector<TrilateralDescriptor>& trilateralDescVec, const float& curvWeight,
	const float& geodesicWeight, const float& areaWeight)
{
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	std::vector< TrilateralError>  errorVector;
	for (size_t i = 0; i < trilateralDescVec.size(); i++)
	{
		TrilateralDescriptor desc_i = trilateralDescVec[i];
		float least_error = INFINITY;
		int least_error_index = -1;
		float least_among_three = 0;
		for (size_t j = 0; j < trilateralDescVec.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			TrilateralDescriptor desc_j = trilateralDescVec[j];
			float areaError = abs(desc_j.area - desc_i.area) / std::max(desc_j.area, desc_i.area);
			areaError = 0;

			//generate 4 vectors ( 2 x 2 )
			// size here is 1 - area 2 - curv_error 3 - curv_error  4 - euclidian_error 5 - euclidian error 6 - curvError 7 - curvError
			Eigen::VectorXf i_1(10);
			Eigen::VectorXf i_2(10);

			Eigen::VectorXf j_1(10);
			Eigen::VectorXf j_2(10);


			//fil those
			i_1(0) = desc_i.area;
			i_2(0) = desc_i.area;

			i_1(1) = desc_i.geodesic_lenght_1_2;
			i_1(2) = desc_i.geodesic_lenght_1_3;

			i_2(1) = desc_i.geodesic_lenght_1_3;
			i_2(2) = desc_i.geodesic_lenght_1_2;

			i_1(3) = 0;//desc_i.euclidian_lenght_1_2;
			i_1(4) = 0;//desc_i.euclidian_lenght_1_3;

			i_2(3) = 0;//desc_i.euclidian_lenght_1_3;
			i_2(4) = 0;//desc_i.euclidian_lenght_1_2;

			i_1(5) = 0;//desc_i.curvature_1_2;
			i_1(6) = 0;//desc_i.curvature_1_3;

			i_2(5) = 0;//desc_i.curvature_1_3;
			i_2(6) = 0;//desc_i.curvature_1_2;

			i_1(7) = desc_i.geodesic_lenght_2_3;
			i_1(8) = 0;//desc_i.euclidian_lenght_2_3;
			i_1(9) = 0;//desc_i.curvature_2_3;

			i_2(7) = desc_i.geodesic_lenght_2_3;
			i_2(8) = 0;//desc_i.euclidian_lenght_2_3;
			i_2(9) = 0;//desc_i.curvature_2_3;

			//fill j 
			j_1(0) = desc_j.area;
			j_2(0) = desc_j.area;

			j_1(1) = desc_j.geodesic_lenght_1_2;
			j_1(2) = desc_j.geodesic_lenght_1_3;

			j_2(2) = desc_j.geodesic_lenght_1_2;
			j_2(1) = desc_j.geodesic_lenght_1_3;

			j_1(3) = 0;//desc_j.euclidian_lenght_1_2;
			j_1(4) = 0;//desc_j.euclidian_lenght_1_3;

			j_2(3) = 0;//desc_j.euclidian_lenght_1_3;
			j_2(4) = 0;//desc_j.euclidian_lenght_1_2;

			j_1(5) = 0;//desc_j.curvature_1_2;
			j_1(6) = 0;//desc_j.curvature_1_3;

			j_2(5) = 0;//desc_j.curvature_1_3;
			j_2(6) = 0;//desc_j.curvature_1_2;

			j_1(7) = desc_j.geodesic_lenght_2_3;
			j_1(8) = 0;//desc_j.euclidian_lenght_2_3;
			j_1(9) = 0;//desc_j.curvature_2_3;

			j_2(7) = desc_j.geodesic_lenght_2_3;
			j_2(8) = 0;//desc_j.euclidian_lenght_2_3;
			j_2(9) = 0;//desc_j.curvature_2_3;

			//normalize 4 vectors
			i_1 = i_1.normalized();
			i_2 = i_2.normalized();
			j_1 = j_1.normalized();
			j_2 = j_2.normalized();

			float dist_i_j_1_1 = sqrt(pow((i_1(0) - j_1(0)), 2) + pow((i_1(1) - j_1(1)), 2) + pow((i_1(2) - j_1(2)), 2) + pow((i_1(3) - j_1(3)), 2) + pow((i_1(4) - j_1(4)), 2)
				+ pow((i_1(5) - j_1(5)), 2) + pow((i_1(6) - j_1(6)), 2) + pow((i_1(7) - j_1(7)), 2) + pow((i_1(8) - j_1(8)), 2) + pow((i_1(9) - j_1(9)), 2));
			float dist_i_j_1_2 = sqrt(pow((i_1(0) - j_2(0)), 2) + pow((i_1(1) - j_2(1)), 2) + pow((i_1(2) - j_2(2)), 2) + pow((i_1(3) - j_2(3)), 2) + pow((i_1(4) - j_2(4)), 2)
				+ pow((i_1(5) - j_2(5)), 2) + pow((i_1(6) - j_2(6)), 2) + pow((i_1(7) - j_2(7)), 2) + pow((i_1(8) - j_2(8)), 2) + pow((i_1(9) - j_2(9)), 2));
			float dist_i_j_2_1 = sqrt(pow((i_2(0) - j_1(0)), 2) + pow((i_2(1) - j_1(1)), 2) + pow((i_2(2) - j_1(2)), 2) + pow((i_2(3) - j_1(3)), 2) + pow((i_2(4) - j_1(4)), 2)
				+ pow((i_2(5) - j_1(5)), 2) + pow((i_2(6) - j_1(6)), 2) + pow((i_2(7) - j_1(7)), 2) + pow((i_2(8) - j_1(8)), 2) + pow((i_2(9) - j_1(9)), 2));
			float dist_i_j_2_2 = sqrt(pow((i_2(0) - j_2(0)), 2) + pow((i_2(1) - j_2(1)), 2) + pow((i_2(2) - j_2(2)), 2) + pow((i_2(3) - j_2(3)), 2) + pow((i_2(4) - j_2(4)), 2)
				+ pow((i_2(5) - j_2(5)), 2) + pow((i_2(6) - j_2(6)), 2) + pow((i_2(7) - j_2(7)), 2) + pow((i_2(8) - j_2(8)), 2) + pow((i_2(9) - j_2(9)), 2));

			if (dist_i_j_1_1 <= dist_i_j_1_2 && dist_i_j_1_1 <= dist_i_j_2_1 && dist_i_j_1_1 <= dist_i_j_2_2)
			{
				least_among_three = dist_i_j_1_1;
			}
			else if (dist_i_j_1_2 <= dist_i_j_1_1 && dist_i_j_1_2 <= dist_i_j_2_1 && dist_i_j_1_2 <= dist_i_j_2_2)
			{
				least_among_three = dist_i_j_1_2;
			}
			if (dist_i_j_2_1 <= dist_i_j_1_2 && dist_i_j_2_1 <= dist_i_j_1_1 && dist_i_j_1_1 <= dist_i_j_2_2)
			{
				least_among_three = dist_i_j_2_1;
			}
			if (dist_i_j_2_2 <= dist_i_j_1_2 && dist_i_j_1_1 <= dist_i_j_2_1 && dist_i_j_2_2 <= dist_i_j_1_1)
			{
				least_among_three = dist_i_j_2_2;
			}
			if (least_among_three < least_error)
			{
				least_error = least_among_three;
				least_error_index = j;
			}
			//check for each point
			//p1
			/*float curvErrorP1P2 = abs(desc_j.curvature_1_2- desc_i.curvature_1_2) / std::max(desc_j.curvature_1_2 , desc_i.curvature_1_2);
			float euclideanErrorP1P2 = abs(desc_j.euclidian_lenght_1_2- desc_i.euclidian_lenght_1_2) / std::max(desc_j.euclidian_lenght_1_2, desc_i.euclidian_lenght_1_2);
			float geodesicErrorP1P2 = abs(desc_j.geodesic_lenght_1_2- desc_i.geodesic_lenght_1_2) / std::max(desc_j.geodesic_lenght_1_2, desc_i.geodesic_lenght_1_2);

			float curvErrorP1P3 = abs(desc_j.curvature_1_3 - desc_i.curvature_1_3) / std::max(desc_j.curvature_1_3, desc_i.curvature_1_3);
			float euclideanErrorP1P3 = abs(desc_j.euclidian_lenght_1_3 - desc_i.euclidian_lenght_1_3) / std::max(desc_j.euclidian_lenght_1_3, desc_i.euclidian_lenght_1_3);
			float geodesicErrorP1P3 = abs(desc_j.geodesic_lenght_1_3 - desc_i.geodesic_lenght_1_3) / std::max(desc_j.geodesic_lenght_1_3, desc_i.geodesic_lenght_1_3);

			//p2
			float curvErrorP2P1 = abs(desc_j.curvature_1_2 - desc_i.curvature_1_2) / std::max(desc_j.curvature_1_2, desc_i.curvature_1_2);
			float euclideanErrorP2P1 = abs(desc_j.euclidian_lenght_1_2 - desc_i.euclidian_lenght_1_2) / std::max(desc_j.euclidian_lenght_1_2, desc_i.euclidian_lenght_1_2);
			float geodesicErrorP2P1 = abs(desc_j.geodesic_lenght_1_2 - desc_i.geodesic_lenght_1_2) / std::max(desc_j.geodesic_lenght_1_2, desc_i.geodesic_lenght_1_2);

			float curvErrorP2P3 = abs(desc_j.curvature_2_3 - desc_i.curvature_2_3) / std::max(desc_j.curvature_2_3, desc_i.curvature_2_3);
			float euclideanErrorP2P3 = abs(desc_j.euclidian_lenght_2_3- desc_i.euclidian_lenght_2_3) / std::max(desc_j.euclidian_lenght_2_3, desc_i.euclidian_lenght_2_3);
			float geodesicErrorP2P3 = abs(desc_j.geodesic_lenght_2_3 - desc_i.geodesic_lenght_2_3) / std::max(desc_j.geodesic_lenght_2_3, desc_i.geodesic_lenght_2_3);

			//p3
			float curvErrorP3P1 = curvErrorP1P3;
			float euclideanErrorP3P1 = euclideanErrorP1P3;
			float geodesicErrorP3P1 = geodesicErrorP1P3;

			float curvErrorP3P2 = curvErrorP2P3;
			float euclideanErrorP3P2 = euclideanErrorP2P3;
			float geodesicErrorP3P2 = geodesicErrorP2P3;



			//check with each
			// p1
			if (least_error > (areaError + curvErrorP1P2 + curvErrorP1P3 + geodesicErrorP1P2 + geodesicErrorP1P3 + euclideanErrorP1P2 + euclideanErrorP1P3))
			{
				if ( desc_j.p1 != desc_i.p2 && desc_j.p2 != desc_i.p1)
				{
					least_error = areaError + curvErrorP1P2 + curvErrorP1P3 + geodesicErrorP1P2 + geodesicErrorP1P3 + euclideanErrorP1P2 + euclideanErrorP1P3;
					least_error_index = desc_j.p1;
				}
			}
			//p2
			if (least_error > (areaError + curvErrorP2P1 + curvErrorP2P3 + geodesicErrorP2P1 + geodesicErrorP2P3 + euclideanErrorP2P1 + euclideanErrorP2P3))
			{
				if (desc_j.p2 != desc_i.p3 && desc_j.p3 != desc_i.p2)
				{
					least_error = areaError + curvErrorP2P1 + curvErrorP2P3 + geodesicErrorP2P1 + geodesicErrorP2P3 + euclideanErrorP2P1 + euclideanErrorP2P3;
					least_error_index = desc_j.p2;
				}
			}
			//p3
			if (least_error > (areaError + curvErrorP3P1 + curvErrorP3P2 + geodesicErrorP3P1 + geodesicErrorP3P2 + euclideanErrorP3P1 + euclideanErrorP3P2))
			{
				if (desc_j.p3 != desc_i.p1 && desc_j.p1 != desc_i.p3)
				{
					least_error = areaError + curvErrorP3P1 + curvErrorP3P2 + geodesicErrorP3P1 + geodesicErrorP3P2 + euclideanErrorP3P1 + euclideanErrorP3P2;
					least_error_index = desc_j.p3;
				}
			}*/

		}
		resemblance_pairs.push_back(std::pair<unsigned int, unsigned int >(trilateralDescVec[i].p1, trilateralDescVec[least_error_index].p1));

	}
	return resemblance_pairs;
}
std::vector<std::pair<unsigned int, unsigned int>>  point_match_trilateral_weights(TrilateralMesh* m, std::vector<TrilateralDescriptor>& trilateralDescVecLeft, std::vector<TrilateralDescriptor>& trilateralDescVecRight, const float& curvWeight,
	const float& geodesicWeight, const float& areaWeight)
{
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	std::vector< TrilateralError>  errorVector;
	for (size_t i = 0; i < trilateralDescVecLeft.size(); i++)
	{
		TrilateralDescriptor desc_i = trilateralDescVecLeft[i];
		float least_error = INFINITY;
		int least_error_index = -1;
		float least_among_three = 0;
		for (size_t j = 0; j < trilateralDescVecRight.size(); j++)
		{

			TrilateralDescriptor desc_j = trilateralDescVecRight[j];
			float areaError = abs(desc_j.area - desc_i.area) / std::max(desc_j.area, desc_i.area);
			areaError = 0;

			//generate 4 vectors ( 2 x 2 )
			// size here is 1 - area 2 - curv_error 3 - curv_error  4 - euclidian_error 5 - euclidian error 6 - curvError 7 - curvError
			Eigen::VectorXf i_1(13);
			Eigen::VectorXf i_2(13);

			Eigen::VectorXf j_1(13);
			Eigen::VectorXf j_2(13);

			//calculate  one ring area ( here for now )

			//fil those
			i_1(0) = desc_i.area;
			i_2(0) = desc_i.area;

			i_1(1) = desc_i.geodesic_lenght_1_2;
			i_1(2) = desc_i.geodesic_lenght_1_3;

			i_2(1) = desc_i.geodesic_lenght_1_3;
			i_2(2) = desc_i.geodesic_lenght_1_2;

			i_1(3) = desc_i.euclidian_lenght_1_2;
			i_1(4) = desc_i.euclidian_lenght_1_3;

			i_2(3) = desc_i.euclidian_lenght_1_3;
			i_2(4) = desc_i.euclidian_lenght_1_2;

			i_1(5) = desc_i.curvature_1_2;
			i_1(6) = desc_i.curvature_1_3;

			i_2(5) = desc_i.curvature_1_3;
			i_2(6) = desc_i.curvature_1_2;

			i_1(7) = desc_i.geodesic_lenght_2_3 ;
			i_1(8) = desc_i.euclidian_lenght_2_3;
			i_1(9) = desc_i.curvature_2_3;

			i_2(7) = desc_i.geodesic_lenght_2_3;
			i_2(8) = desc_i.euclidian_lenght_2_3;
			i_2(9) = desc_i.curvature_2_3;


			i_1(10) = 0;//desc_i.n_ring_area_p1;
			i_1(11) = 0;//desc_i.n_ring_area_p2;
			i_1(12) = 0;//desc_i.n_ring_area_p3;
			i_2(10) = 0;//desc_i.n_ring_area_p1;
			i_2(11) = 0;//desc_i.n_ring_area_p3;
			i_2(12) = 0;//desc_i.n_ring_area_p2;
			//fill j 
			j_1(0) = desc_j.area;
			j_2(0) = desc_j.area;

			j_1(1) = desc_j.geodesic_lenght_1_2;
			j_1(2) = desc_j.geodesic_lenght_1_3;

			j_2(2) = desc_j.geodesic_lenght_1_2;
			j_2(1) = desc_j.geodesic_lenght_1_3;

			j_1(3) = desc_j.euclidian_lenght_1_2;
			j_1(4) = desc_j.euclidian_lenght_1_3;

			j_2(3) = desc_j.euclidian_lenght_1_3;
			j_2(4) = desc_j.euclidian_lenght_1_2;

			j_1(5) = desc_j.curvature_1_2;
			j_1(6) = desc_j.curvature_1_3;

			j_2(5) = desc_j.curvature_1_3;
			j_2(6) = desc_j.curvature_1_2;

			j_1(7) = desc_j.geodesic_lenght_2_3;
			j_1(8) = desc_j.euclidian_lenght_2_3;
			j_1(9) = desc_j.curvature_2_3;

			j_2(7) = desc_j.geodesic_lenght_2_3;
			j_2(8) = desc_j.euclidian_lenght_2_3;
			j_2(9) = desc_j.curvature_2_3;

			j_1(10) = 0;//desc_j.n_ring_area_p1;
			j_1(11) = 0;//desc_j.n_ring_area_p2;
			j_1(12) = 0;//desc_j.n_ring_area_p3;
			j_2(10) = 0;//desc_j.n_ring_area_p1;
			j_2(11) = 0;//desc_j.n_ring_area_p3;
			j_2(12) = 0;//desc_j.n_ring_area_p2;
			//normalize 4 vectors
			i_1 = i_1.normalized();
			i_2 = i_2.normalized();
			j_1 = j_1.normalized();
			j_2 = j_2.normalized();

			float dist_i_j_1_1 = sqrt(pow((i_1(0) - j_1(0)), 2) + pow((i_1(1) - j_1(1)), 2) + pow((i_1(2) - j_1(2)), 2) + pow((i_1(3) - j_1(3)), 2) + pow((i_1(4) - j_1(4)), 2)
				+ pow((i_1(5) - j_1(5)), 2) + pow((i_1(6) - j_1(6)), 2) + pow((i_1(7) - j_1(7)), 2) + pow((i_1(8) - j_1(8)), 2) + pow((i_1(9) - j_1(9)), 2)
				+ pow((i_1(10) - j_1(10)), 2) + +pow((i_1(10) - j_1(10)), 2) + pow((i_1(10) - j_1(10)), 2));
			float dist_i_j_1_2 = sqrt(pow((i_1(0) - j_2(0)), 2) + pow((i_1(1) - j_2(1)), 2) + pow((i_1(2) - j_2(2)), 2) + pow((i_1(3) - j_2(3)), 2) + pow((i_1(4) - j_2(4)), 2)
				+ pow((i_1(5) - j_2(5)), 2) + pow((i_1(6) - j_2(6)), 2) + pow((i_1(7) - j_2(7)), 2) + pow((i_1(8) - j_2(8)), 2) + pow((i_1(9) - j_2(9)), 2)
				+ pow((i_1(10) - j_1(10)), 2) + +pow((i_1(10) - j_1(10)), 2) + pow((i_1(10) - j_1(10)), 2));
			float dist_i_j_2_1 = sqrt(pow((i_2(0) - j_1(0)), 2) + pow((i_2(1) - j_1(1)), 2) + pow((i_2(2) - j_1(2)), 2) + pow((i_2(3) - j_1(3)), 2) + pow((i_2(4) - j_1(4)), 2)
				+ pow((i_2(5) - j_1(5)), 2) + pow((i_2(6) - j_1(6)), 2) + pow((i_2(7) - j_1(7)), 2) + pow((i_2(8) - j_1(8)), 2) + pow((i_2(9) - j_1(9)), 2)
				+ pow((i_1(10) - j_1(10)), 2) + +pow((i_1(10) - j_1(10)), 2) + pow((i_1(10) - j_1(10)), 2));
			float dist_i_j_2_2 = sqrt(pow((i_2(0) - j_2(0)), 2) + pow((i_2(1) - j_2(1)), 2) + pow((i_2(2) - j_2(2)), 2) + pow((i_2(3) - j_2(3)), 2) + pow((i_2(4) - j_2(4)), 2)
				+ pow((i_2(5) - j_2(5)), 2) + pow((i_2(6) - j_2(6)), 2) + pow((i_2(7) - j_2(7)), 2) + pow((i_2(8) - j_2(8)), 2) + pow((i_2(9) - j_2(9)), 2)
				+ pow((i_1(10) - j_1(10)), 2) + +pow((i_1(10) - j_1(10)), 2) + pow((i_1(10) - j_1(10)), 2));
			if (dist_i_j_1_1 <= dist_i_j_1_2 && dist_i_j_1_1 <= dist_i_j_2_1 && dist_i_j_1_1 <= dist_i_j_2_2)
			{
				least_among_three = dist_i_j_1_1;
			}
			else if (dist_i_j_1_2 <= dist_i_j_1_1 && dist_i_j_1_2 <= dist_i_j_2_1 && dist_i_j_1_2 <= dist_i_j_2_2)
			{
				least_among_three = dist_i_j_1_2;
			}
			if (dist_i_j_2_1 <= dist_i_j_1_2 && dist_i_j_2_1 <= dist_i_j_1_1 && dist_i_j_1_1 <= dist_i_j_2_2)
			{
				least_among_three = dist_i_j_2_1;
			}
			if (dist_i_j_2_2 <= dist_i_j_1_2 && dist_i_j_1_1 <= dist_i_j_2_1 && dist_i_j_2_2 <= dist_i_j_1_1)
			{
				least_among_three = dist_i_j_2_2;
			}
			if (least_among_three < least_error)
			{
				least_error = least_among_three;
				least_error_index = j;
			}
			

		}
		resemblance_pairs.push_back(std::pair<unsigned int, unsigned int >(trilateralDescVecLeft[i].p1, trilateralDescVecRight[least_error_index].p1));

	}
	return resemblance_pairs;
}
 void display_accuracy(TrilateralMesh* m, std::vector<std::pair<unsigned int, unsigned int>>& calculated_symmetry_pairs)
{
	float correct_sample = 0;
	float total_sample = calculated_symmetry_pairs.size();

	std::vector<std::pair<unsigned int, unsigned int>> relevant_symmetry_pairs;
	std::vector<bool> is_sample_correct(calculated_symmetry_pairs.size(), false);
	// 1 - first fect the relevant symmetry pairs 
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		for (size_t j = 0; j < calculated_symmetry_pairs.size(); j++)
		{
			if ((m->symmetry_pairs[i].first == calculated_symmetry_pairs[j].first
				&& m->symmetry_pairs[i].second == calculated_symmetry_pairs[j].second)
				|| (m->symmetry_pairs[i].second == calculated_symmetry_pairs[j].first
					&& m->symmetry_pairs[i].first == calculated_symmetry_pairs[j].second))
			{
				correct_sample += 1;
				is_sample_correct[j] = true;
				break;
			}
		}
	}

	std::vector<float> pairdata;
	glBindBuffer(GL_ARRAY_BUFFER, 3); // VBO will lok into that later 
	for (int i = 0; i < calculated_symmetry_pairs.size(); i++)
	{
		glm::vec3 p1 = m->vertices[calculated_symmetry_pairs[i].first];
		glm::vec3 p2 = m->vertices[calculated_symmetry_pairs[i].second];
		pairdata.push_back(p1.x);
		pairdata.push_back(p1.y);
		pairdata.push_back(p1.z);

		if (is_sample_correct[i])
		{
			pairdata.push_back(255);
			pairdata.push_back(255);
			pairdata.push_back(255);
		}
		else
		{
			pairdata.push_back(255);
			pairdata.push_back(0);
			pairdata.push_back(0);
		}



		pairdata.push_back(p2.x);
		pairdata.push_back(p2.y);
		pairdata.push_back(p2.z);
		if (is_sample_correct[i])
		{
			pairdata.push_back(255);
			pairdata.push_back(255);
			pairdata.push_back(255);
		}
		else
		{
			pairdata.push_back(255);
			pairdata.push_back(0);
			pairdata.push_back(0);
		}

	}
	glBufferData(GL_ARRAY_BUFFER, pairdata.size() * sizeof(float), &pairdata[0], GL_STATIC_DRAW);
	MeshPointPairs p;
	p.point_pairs = pairdata;

	std::cout << " total accuracy " << correct_sample / total_sample << std::endl;
}

 void point_matching_with_dominant_symmetry_plane(MeshFactory& mesh_fac, int& selected_index, Plane* plane, int sampling_no)
{
	TrilateralMesh* mesh = &mesh_fac.mesh_vec[selected_index];

	std::vector<unsigned int> vertex_indices = furthest_point_sampling(mesh, sampling_no);
	// divide the points in two 
	std::vector<int> vertex_indices_right;
	std::vector<int> vertex_indices_left;

	float* vertex_indices_left_right = new float[vertex_indices.size()];
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
		std::vector<float> closest_matrix = Geodesic_dijkstra(*mesh, vertex_indices[i]);
		//get closest 3 
		int closest_index_1 = -1;
		int closest_index_2 = -1;
		float closest_distance1 = 1e5;
		float closest_distance2 = 1e5;
		for (size_t j = 0; j < vertex_indices.size(); j++)
		{
			if (closest_distance1 > closest_matrix[vertex_indices[j]])
			{
				if (j != i) //because it is 0 an himself
				{
					if (((vertex_indices_left_right[i] >= 0) && (vertex_indices_left_right[j] >= 0)) || ((vertex_indices_left_right[i] < 0) && (vertex_indices_left_right[j] < 0))) //check same sign 
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
				if (j != i && j != closest_index_1) //because it is 0 an himself
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
		descriptors[i] = generate_trilateral_descriptor(mesh_fac, selected_index, closest_3_pairs[i][0], closest_3_pairs[i][1], closest_3_pairs[i][2], true);
	}

	std::vector<std::vector<int>> vertex_index_pairs;
	//last part match points.
	for (size_t i = 0; i < vertex_indices.size(); i++)
	{
		float min_diff_score = INFINITY;
		int min_diff_index = -1;
		for (size_t j = 0; j < vertex_indices.size(); j++)
		{
			if (i != j)
			{
				float area_dif = std::abs(descriptors[i].area - descriptors[j].area);
				float roi_len_dif = std::abs(descriptors[i].total_length - descriptors[j].total_length);
				// easy sum for now
				float curv_dif = std::abs(descriptors[i].curvature_1_2 + descriptors[i].curvature_1_3 + descriptors[i].curvature_2_3 - descriptors[j].curvature_1_2 - descriptors[j].curvature_1_3 - descriptors[j].curvature_2_3);
				float dist_dif = std::abs(descriptors[i].geodesic_lenght_1_2 + descriptors[i].geodesic_lenght_1_3 + descriptors[i].geodesic_lenght_2_3 - descriptors[j].geodesic_lenght_1_2 - descriptors[j].geodesic_lenght_1_3 - descriptors[j].geodesic_lenght_2_3);

				if (area_dif + roi_len_dif + curv_dif + dist_dif < min_diff_score)
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
	TrilateralMesh m1 = mesh_fac.mesh_vec[mesh_index1];
	TrilateralMesh m2 = mesh_fac.mesh_vec[mesh_index2];

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
						//histogram_arr = trilateral_ROI(mesh_fac, mesh_index1, i, k, j, division_no); // tau is the closeness division_no is the no of how much you want to separate
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
//	TrilateralMesh* m1 = &mesh_fac.mesh_vec[mesh_index1];
//	TrilateralMesh* m2 = &mesh_fac.mesh_vec[mesh_index2];
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
//					int* is_visited = trilateral_ROI(mesh_fac, mesh_index1, m1_histograms[i].first, m1_histograms[j].first, m1_histograms[k].first, division_no, is_visited_interior);
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
//					int* is_visited = trilateral_ROI(mesh_fac, mesh_index2, m2_histograms[i].first, m2_histograms[j].first, m2_histograms[k].first, division_no, is_visited_interior);
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
/*void match_points_from2_mesh(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int division_no) // sample size must be bigger than 3 
{
	//get the histogram of the traingle with the most near 2 points
	TrilateralMesh* m1 = &mesh_fac.mesh_vec[mesh_index1];
	TrilateralMesh* m2 = &mesh_fac.mesh_vec[mesh_index2];

	std::vector<int> sample_indices_m1;
	std::vector<int> sample_indices_m2;
	srand(time(NULL));
	int set[10] = { 0 ,100 , 200 , 1500 , 2500,3200 , 5000 , 4000 , 400 ,4200 };
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
		std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m1, p1_index);
		//traverse through distance matrix
		for (size_t j = 0; j < sample_indices_m1.size(); j++)
		{
			float distance = distance_matrix_p1[sample_indices_m1[j]];
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
		int* is_visited = trilateral_ROI(mesh_fac, mesh_index2, sample_indices_m1[i], sample_indices_m1[first_min], sample_indices_m1[second_min], division_no, is_visited_interior);
		std::vector<float> histogram_of_nearest_2_points = histogramROi(mesh_fac, mesh_index1, sample_indices_m1[i], sample_indices_m1[first_min], sample_indices_m1[second_min], division_no, is_visited, is_visited_interior);

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
		std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m1, p1_index);
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
		int* is_visited = trilateral_ROI(mesh_fac, mesh_index2, sample_indices_m2[i], sample_indices_m2[first_min], sample_indices_m2[second_min], division_no, is_visited_interior);
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

		final_pairs.push_back(std::make_pair(i, min_index));
	}
	glm::vec3 red = glm::vec3(1.0f, 0.0f, 0.0f);
	glm::vec3 green = glm::vec3(0.0f, 1.0f, 0.0f);
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
}*/

//from the paper Robust3DShapeCorrespondenceintheSpectralDomain 4.1 
 std::vector<glm::vec3> generate_spectral_embedding(MeshFactory& meshFac, int mesh_index, std::vector<unsigned int> landmark_vertices)
{
	int k = 3; // 2 , 3 and  4th eignvectors 
	TrilateralMesh* m = &meshFac.mesh_vec[mesh_index];
	std::sort(landmark_vertices.begin(), landmark_vertices.end());
	Eigen::MatrixXf affinity_A(landmark_vertices.size(), landmark_vertices.size());
	affinity_A.setZero();
	// 1 - generate affinitiy matrix
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, landmark_vertices[i]);
		int landmark_vertex_index = 0;
		for (size_t j = 0; j < m->vertices.size(); j++)
		{
			if (j == landmark_vertices[landmark_vertex_index])
			{
				affinity_A(i, landmark_vertex_index) = distance_matrix_p1[j];
				if (landmark_vertex_index + 1 == landmark_vertices.size())
				{
					break;
				}
				landmark_vertex_index++;

			}
		}

	}

	//diagonal 0 because diagonal is distance to itself
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		affinity_A(i, i) = 0;
	}
	// R = AAt
	// Eigen::MatrixXf R(landmark_vertices.size(), landmark_vertices.size());
	// R = affinity_A * affinity_A.transpose();

	//solves the eignevector
	Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
	eigensolver.compute(affinity_A);
	Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
	Eigen::MatrixXf eigen_vectors = eigensolver.eigenvectors().real();

	/*std::cout << "eigen values before sorted  " << std::endl;
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		std::cout << eigen_values(i) << std::endl;
	}*/
	//sort them to descending order !!!!
	std::vector<unsigned int> sorted_eigen_vertices(landmark_vertices.size());
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		int biggest_index = -1;
		float biggest_value = -INFINITY;

		for (size_t j = i; j < landmark_vertices.size(); j++)
		{
			if (eigen_values(j) > biggest_value)
			{
				biggest_value = eigen_values(j);
				biggest_index = j;
			}
		}
		sorted_eigen_vertices[i] = biggest_index;
		float temp = eigen_values(i);
		eigen_values(i) = eigen_values(biggest_index);
		eigen_values(biggest_index) = temp;
	}
	/*std::cout << "eigen values sorted  " << std::endl;
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		std::cout << eigen_values(i) << std::endl;
	}*/
	// now sort eigen_vectors
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		Eigen::MatrixXf temp = eigen_vectors.col(i);
		eigen_vectors.col(i) = eigen_vectors.col(sorted_eigen_vertices[i]);
		eigen_vectors.col(sorted_eigen_vertices[i]) = temp;
	}


	//get get 2 3 and 4'th eigenvector 
	Eigen::MatrixXf first_k_eigen_vectors(landmark_vertices.size(), k);
	for (size_t i = 0; i < k; i++)
	{
		first_k_eigen_vectors.col(i) = eigen_vectors.col(i + 1); // +1 is because of 2 , 3 and 4 th eigenvectors 
	}
	//make eigenvalues in diagonal  matrix form
	Eigen::MatrixXf eigen_values_diag_matrix(k, k);
	eigen_values_diag_matrix.setZero();
	for (size_t i = 0; i < k; i++)
	{
		eigen_values_diag_matrix(i, i) = eigen_values(i + 1);// +1 is because of 2 , 3 and 4 th eigenvectors
	}

	// fill the landmark vector 
	std::vector<glm::vec3> embedded_points_vec(m->vertices.size());
	std::fill(embedded_points_vec.begin(), embedded_points_vec.end(), glm::vec3(0, 0, 0));

	Eigen::MatrixXf A_k(k, landmark_vertices.size());
	A_k = eigen_values_diag_matrix * first_k_eigen_vectors.transpose();

	// end of 4.1 


	std::cout << " eigen values diag matrix " << std::endl;
	/*for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			std::cout << i << " " << j << " " << eigen_values_diag_matrix(i, j) << std::endl;
		}
	}
	std::cout << " first k eigen vectors " << std::endl;
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			std::cout << i << " " << j << " " << first_k_eigen_vectors(i, j) << std::endl;
		}
	}
	std::cout << "A_k" << std::endl;
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < landmark_vertices.size(); j++)
		{
			std::cout << i << " " << j << " " << A_k(i, j) << std::endl;
		}
	}*/
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		embedded_points_vec[landmark_vertices[i]] = glm::vec3(A_k(0, i), A_k(1, i), A_k(2, i));
	}
	//now landmark points have been filled

	// continuing with Sparse multidimensional scaling
	// using landmark points Silva et al 2004 
	Eigen::MatrixXf sigma(landmark_vertices.size(), m->vertices.size()); //sigma
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		glm::vec3 p1 = m->vertices[i];
		for (size_t j = 0; j < landmark_vertices.size(); j++)
		{
			float dist = glm::distance(p1, embedded_points_vec[landmark_vertices[j]]);
			sigma(j, i) = dist;
		}
	}
	//calculate mean of sigma
	Eigen::VectorXf mean_sigma(landmark_vertices.size());
	mean_sigma.setZero();
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		mean_sigma += sigma.col(i);
	}
	mean_sigma /= m->vertices.size();

	/*std::cout << " mean sigma " << std::endl;
	for (size_t i = 0; i < landmark_vertices.size(); i++)
	{
		std::cout << mean_sigma(i) << std::endl;
	}*/


	Eigen::MatrixXf L_pseudo_inv(k, landmark_vertices.size());
	//fill l pseudo inv 
	for (size_t i = 0; i < k; i++)
	{
		L_pseudo_inv.row(i) = first_k_eigen_vectors.col(i) * sqrt(eigen_values_diag_matrix(i, i));

		/*for (size_t j = 0; j < landmark_vertices.size(); j++)
		{
			L_pseudo_inv(i, j) = first_k_eigen_vectors(j , i ) / sqrt(eigen_values(i));
			std::cout << first_k_eigen_vectors(j, i) << " " << sqrt(eigen_values(i)) << std::endl;;
		}*/
	}
	/*std::cout << " L_pseudor inv " << std::endl;
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < landmark_vertices.size(); j++)
		{
			std::cout << i << " " << j << " " << L_pseudo_inv(i, j) << std::endl;
		}
	}*/
	//for each point calculate the embedding vector
	int landmark_vertices_index = 0;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		Eigen::MatrixXf embed_vec = L_pseudo_inv * (sigma.col(i) - mean_sigma);
		embedded_points_vec[i] = glm::vec3(embed_vec(0, 0), embed_vec(1, 0), embed_vec(2, 0));
	}

	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		m->vertices[i] = glm::vec3(embedded_points_vec[i].x * m->vertices[i].x, embedded_points_vec[i].y * m->vertices[i].y, embedded_points_vec[i].z * m->vertices[i].z);
	}

	// try to create these vecs
	//m->vertices = embedded_points_vec;
	return embedded_points_vec;
}

 void reset_points(MeshFactory& mesh_fac, int meshIndex)
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[meshIndex];

	for (size_t i = 0; i < m->colors.size(); i++)
	{
		m->colors[i].r = 0;
		m->colors[i].g = 0;
		m->colors[i].b = 0;
	}
}

 float get_N_ring_area(TrilateralMesh* m, float point_index , int N )
 {
	 float area = 0;
	 std::vector<unsigned int> indices;
	 std::vector<bool> is_vertex_on_ring(m->vertices.size() , false);
	 if (N == 1)
	 {
		 return m->areas[point_index];
	 }
	 for (size_t i = 0; i < N+1; i++)
	 {
		 if (i == 0) // one ring 
		 {
			 indices.push_back(point_index);
			 is_vertex_on_ring[point_index] = true;
			 continue;
		 }
		 int static_indices_size = indices.size();
		 for (size_t j = 0; j < static_indices_size; j++)
		 {
			 for (size_t k = 0; k < m->adjacenies[indices[j]].size(); k++)
			 {
				 if (!is_vertex_on_ring[m->adjacenies[indices[j]][k].first] )
				 {
					 is_vertex_on_ring[m->adjacenies[indices[j]][k].first] = true; 
					 indices.push_back(m->adjacenies[indices[j]][k].first);
				 }
			 }
		 }
	 }
	 for (size_t i = 0; i < m->triangles.size(); i += 3)
	 {
		 if(is_vertex_on_ring[m->triangles[i]] || is_vertex_on_ring[m->triangles[i + 1]] || is_vertex_on_ring[m->triangles[i + 2]] )
		 {
			 area += compute_triangle_area(m->vertices[m->triangles[i]], m->vertices[m->triangles[i + 1]], m->vertices[m->triangles[i + 2]]);
		 }
	 }
	 return area;
 }

 void trilateral_self_matching_with_dominant_sym(MeshFactory& mesh_fac,int selected_index,  const int number_of_n_lateral_points , const std::string& method_name , int sym_iter_no, bool is_LMDS)
{
	 TrilateralMesh mesh = mesh_fac.mesh_vec[selected_index];
	 // 1 - part one is same for now, generate a symmetry plane 
#pragma region symmetry plane 
	 if (is_LMDS)
	 {
		 mesh = compute_landmark_MDS(&mesh, 3);
	 }
	 //calculate center of the plane 
	 glm::vec3 plane_center(0, 0, 0);
	 for (size_t i = 0; i < mesh.vertices.size(); i++)
	 {
		 plane_center += mesh.vertices[i];
	 }
	 plane_center /= mesh.vertices.size();
	 Plane plane = generate_dominant_symmetry_plane(&mesh , sym_iter_no);

#pragma endregion

#pragma region separation of points on the mesh from plane 
	 std::vector<unsigned int> fps_points = furthest_point_sampling(&mesh, number_of_n_lateral_points, true);
#pragma endregion

#pragma region do point sampling on partial regions 

	 std::vector<unsigned int > fps_positive;
	 std::vector<unsigned int > fps_negative;
	 // now do two distinct fps
	 
	 for (size_t i = 0; i < number_of_n_lateral_points; i++)
	 {
		 if (get_point_status_from_plane(&plane, &mesh.vertices[fps_points[i]]) > 0)
		 {
			 fps_positive.push_back(fps_points[i]);
		 }
		 else
		 {
			 fps_negative.push_back(fps_points[i]);
		 }

	 }
	 
#pragma endregion

#pragma region n_lateral_algorithm
		
	 // n_lateral computation
	 
	 std::vector<TrilateralDescriptor> positive_mesh_trilateral_descriptor = get_trilateral_points_using_closest_pairs(&mesh, fps_positive);
	 std::vector<TrilateralDescriptor> negative_mesh_trilateral_descriptor = get_trilateral_points_using_closest_pairs(&mesh, fps_negative);
#pragma endregion

	 for (size_t i = 0; i < fps_positive.size(); i++)
	 {
		 std::vector<int> global_is_visited(mesh.vertices.size(), OUTSIDE);
		 ROI_trilateral(&mesh, positive_mesh_trilateral_descriptor[i], 3, false);
		 Histogram histogram = Histogram_triangle_area(&mesh, positive_mesh_trilateral_descriptor[i],3);
		 positive_mesh_trilateral_descriptor[i].histogram = histogram;
	 }
	 for (size_t i = 0; i < fps_negative.size(); i++)
	 {
		 std::vector<int> global_is_visited(mesh.vertices.size(), OUTSIDE);
		 ROI_trilateral(&mesh, negative_mesh_trilateral_descriptor[i], 3, false);
		 Histogram histogram = Histogram_triangle_area(&mesh, negative_mesh_trilateral_descriptor[i], 3);
		 negative_mesh_trilateral_descriptor[i].histogram = histogram;
	 }
	 // write a function for comparing two descriptor
	 //irrelevant constants 
	 std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = trilateral_unique_pairing_histogram(positive_mesh_trilateral_descriptor,
		 negative_mesh_trilateral_descriptor);//point_match_trilateral_weights(&L_MDS_mesh, positive_mesh_trilateral_descriptor, negative_mesh_trilateral_descriptor, const1, const2, const3);

	 //forge it into two list
	 std::vector<unsigned int> left_correspondences;
	 std::vector<unsigned int> right_correspondences;
	 for (size_t i = 0; i < resemblance_pairs.size(); i++)
	 {
		 left_correspondences.push_back(resemblance_pairs[i].first);
		 right_correspondences.push_back(resemblance_pairs[i].second);
	 }
	 float total_error = Metric_get_geodesic_cost_with_list(&mesh, left_correspondences, right_correspondences);

	 // now use fps points to get maximum distance in order to compare to 
	 float maximum_geodesic_distance = 0;
	 for (size_t i = 0; i < fps_positive.size(); i++)
	 {
		 std::vector<float> distances = Geodesic_dijkstra(mesh, fps_positive[i]);
		 for (size_t j = 0; j < distances.size(); j++)
		 {
			 if (maximum_geodesic_distance < distances[j])
			 {
				 maximum_geodesic_distance = distances[j];
			 }
		 }
	 }

	 // color left red
	 std::vector<unsigned int> is_selected(mesh.vertices.size(), 0);
	 for (size_t i = 0; i < resemblance_pairs.size(); i++)
	 {
		 mesh.colors[resemblance_pairs[i].first].r = 255;
		 mesh.colors[resemblance_pairs[i].first].g = 0;
		 mesh.colors[resemblance_pairs[i].first].b = 0;

		 mesh.colors[resemblance_pairs[i].second].r = 0;
		 mesh.colors[resemblance_pairs[i].second].g = 0;
		 mesh.colors[resemblance_pairs[i].second].b = 255;
	 }

	 mesh.calculated_symmetry_pairs = resemblance_pairs;

	 Metric_write_to_file(&mesh, "../../Results/Trilateral_W_LMDS_AND_DOMINANTPLANE.txt");

	 TrilateralMesh plane_mesh = generate_mesh_from_plane(&plane, &plane.point);
	 mesh_fac.add_mesh(plane_mesh);

	 mesh_fac.mesh_vec[selected_index] = mesh ;

}


static std::vector<std::pair<unsigned int, unsigned int>>  trilateral_unique_pairing_histogram(std::vector<TrilateralDescriptor>& trilateral_desc,
int N,std::vector<Histogram>& trilateral_histograms)
{
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	//compare each
	std::vector<std::pair<float, std::pair<int, int>>> compareResults;

	for (size_t i = 0; i < N; i++)
	{
		Eigen::VectorXd histogram_i = stdVectorToEigenVectorXd(trilateral_histograms[i].histogram);
		int minimum_index = -1;
		float minimum_value = INFINITY;
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				continue;
			}

			Eigen::VectorXd histogram_j = stdVectorToEigenVectorXd(trilateral_histograms[j].histogram);
			//float resemblance = (histogram_i - histogram_j).norm(); // L2 norm
			float resemblance = histogram_i.dot(histogram_j) / (histogram_i.norm() * histogram_j.norm()); // cosine similarity 
			
			compareResults.push_back({ resemblance, {i,j} });

		}
	}

	std::vector<std::pair<int, int>> selectedPairs;
	float skeleton_resemblance_error = 0;
	std::vector<bool> used(N, false);
	std::sort(compareResults.begin(), compareResults.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ trilateral_desc[i].p1 , trilateral_desc[j].p1 });
			used[i] = used[j] = true;  // Mark these objects as used
		}
	}
	return resemblance_pairs;

}
std::vector<float>  histogramROi_hks(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited)
{
	//histogram to be returned 
	std::vector<float> histogram;
	// fill it with division no 
	for (size_t i = 0; i < division_no; i++)
	{
		histogram.push_back(0);
	}
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
		float heat_kernel = ( m->normalized_heat_kernel_signature[i1] + m->normalized_heat_kernel_signature[i2] + m->normalized_heat_kernel_signature[i3] ) / 3;
		float triangle_area = ( m->normalized_heat_kernel_signature[i1] + m->normalized_heat_kernel_signature[i2] + m->normalized_heat_kernel_signature[i3] ) / 3;
		if (step_no_i1 == step_no_i2 && step_no_i1 == step_no_i3)
		{
			histogram[step_no_i1] += heat_kernel  *  triangle_area;
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

			histogram[steps[0].first] += area_arc_small * m->normalized_heat_kernel_signature[steps[0].second];


			// 1.1.2 closest arc
			//a miscalculation but think it like also an arc although it is a reverse arc
			float dist_step_2 = distance_matrix_p1[steps[2].second];
			float big_r = dist_step_2 - (steps[2].first * step);
			edge1 = m->vertices[steps[1].second] - m->vertices[steps[2].second];
			edge2 = m->vertices[steps[0].second] - m->vertices[steps[2].second];
			cosine = glm::dot(edge1, edge2) / (glm::length(edge1) * glm::length(edge2));
			small_radian = acos(cosine); //radian
			float area_arc_big = M_PI * big_r * big_r * small_radian / (2 * M_PI);
			histogram[steps[2].first] += area_arc_big * m->normalized_heat_kernel_signature[steps[2].second];


			float area_left = triangle_area - (area_arc_small + area_arc_big);

			int steps_left = steps[2].first - (steps[0].first + 1);
			if (steps_left == 0)
			{
				histogram[steps[0].first] += area_left * m->normalized_heat_kernel_signature[steps[0].second];
				histogram[steps[1].first] += area_left * m->normalized_heat_kernel_signature[steps[1].second];
			}
			else
			{
				int total_fraction = steps_left * steps_left;
				for (size_t step_no = steps[0].first + 1; step_no < steps[2].first; step_no++)
				{
					int fraction_no = step_no - (steps[0].first);
					float fraction_percentage = (pow(fraction_no, 2) - pow(fraction_no - 1, 2)) / fraction_no;

					histogram[step_no] += fraction_percentage * area_left * heat_kernel;
				}
			}
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

static std::vector<std::vector<Pair_w>> trilateral_hks_scalar(TrilateralMesh* m , 
std::vector<TrilateralDescriptor>& trilateral_desc_left,
std::vector<TrilateralDescriptor>& trilateral_desc_right,
int N ) // return N of the most resembled for each point
{
	std::vector<std::pair< float, std::pair<unsigned int, unsigned int>>> resemblance_pairs_left_right;
	std::vector<std::pair<float, std::pair<int, int>>> compareResults_left;

	int N_left = trilateral_desc_left.size();
	int N_right = trilateral_desc_right.size();
	int histogram_size = trilateral_desc_left[0].histogram.size();
	std::vector<std::vector<Pair_w>> pair_vec( N_left + N_right );
	for (size_t i = 0; i < N_left; i++)
	{
		float hks_i = m->normalized_heat_kernel_signature[trilateral_desc_left[i].p1];
		for (size_t j = 0; j < N_right; j++)
		{
			float hks_j = m->normalized_heat_kernel_signature[trilateral_desc_right[j].p1];
			
			float resemblance = std::abs(hks_i - hks_j); // L2 norm
			//float resemblance_c = histogram_i.dot(histogram_j) / (histogram_i.norm() * histogram_j.norm()); // cosine similarity 
			compareResults_left.push_back({ resemblance, {i,j} });
		}
	}
	// Greedily select pairs for left 
	for (const auto& entry : compareResults_left) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (pair_vec[i].size() < N  ) {
			Pair_w w = { i , j , entry.first};
			pair_vec[i].push_back(w);
		}
	}
	
}
static std::vector<std::pair<unsigned int, unsigned int>> trilateral_hks_histogram(TrilateralMesh* m,
	std::vector<TrilateralDescriptor>& trilateral_desc_left,
	std::vector<TrilateralDescriptor>& trilateral_desc_right,
	int histogram_size,
	int N) // return N of the most resembled for each point
{
	std::vector<std::pair< float, std::pair<unsigned int, unsigned int>>> resemblance_pairs_left_right;
	std::vector<std::pair<float, std::pair<int, int>>> compareResults_left;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs; 
	int N_left = trilateral_desc_left.size();
	int N_right = trilateral_desc_right.size();
	int sample_no = N_left + N_right;

	std::vector<int> global_is_visited(m->vertices.size(), OUTSIDE);
	std::vector<Histogram> histograms_left;
	std::vector<Histogram> histograms_right;
	for (size_t i = 0; i <  N_left; i++)
	{
		Histogram histogram = Histogram_triangle_area_w_res(m, trilateral_desc_left[i], histogram_size, 4);
		histograms_left.push_back(histogram);
	}
	for (size_t i = 0; i < N_right; i++)
	{
		Histogram histogram = Histogram_triangle_area_w_res(m, trilateral_desc_right[i], histogram_size,4);
		histograms_right.push_back(histogram);
	}

	for (size_t i = 0; i < N_left; i++)
	{
		int minimum_index = -1;
		float minimum_value = INFINITY;
		Histogram histogram_vec_i = histograms_left[i];

		for (size_t j = 0; j < N_right; j++)
		{
			Histogram histogram_vec_j = histograms_right[j];
			float resemblance = Histogram_L2Norm(histogram_vec_i, histogram_vec_j);
			
			float resemblance_1_2 = Histogram_L2Norm_DifferentSize(trilateral_desc_left[i].hist_path_1_2, trilateral_desc_right[j].hist_path_1_2);
			float resemblance_1_3 = Histogram_L2Norm_DifferentSize(trilateral_desc_left[i].hist_path_1_3, trilateral_desc_right[j].hist_path_1_3);
			float resemblance_2_3 = Histogram_L2Norm_DifferentSize(trilateral_desc_left[i].hist_path_2_3, trilateral_desc_right[j].hist_path_2_3);

			if (resemblance_1_2 + resemblance_1_3 + resemblance_2_3 < 0.3)
			{
				compareResults_left.push_back({ resemblance, {i,j} });
			}

		}
	}

	std::vector<std::pair<int, int>> selectedPairs;
	std::vector<bool> used_left_l(N_left, false);
	std::vector<bool> used_right_l(N_right, false);
	std::sort(compareResults_left.begin(), compareResults_left.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults_left) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_left_l[i]) 
		{
			std::pair<float, std::pair<unsigned int, unsigned int>> pair_with_weight;
			pair_with_weight.first = entry.first;
			pair_with_weight.second = std::pair<unsigned int, unsigned int>(i, j);
			resemblance_pairs_left_right.push_back(pair_with_weight);
			used_left_l[i] = used_right_l[j] = true;  // Mark these objects as used
		}
	}
	std::vector<std::pair<float, std::pair<int, int>>> compareResults_right;
	//do it for right to let also

	for (size_t i = 0; i < N_right; i++)
	{
		int minimum_index = -1;
		float minimum_value = INFINITY;
		Histogram histogram_vec_i = histograms_right[i];
		for (size_t j = 0; j < N_left; j++)
		{
			Histogram histogram_vec_j = histograms_left[j];
			float resemblance = Histogram_L2Norm(histogram_vec_i, histogram_vec_j);
			;			compareResults_left.push_back({ resemblance, {i,j} });
		}
	}
	std::vector<std::pair<int, int>> selectedPairs_r;
	std::vector<bool> used_left_r(N_left, false);
	std::vector<bool> used_right_r(N_right, false);
	std::sort(compareResults_right.begin(), compareResults_right.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults_right) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_right_r[i] ) 
		{
			std::pair<float, std::pair<unsigned int, unsigned int >> pair_with_weight;
			pair_with_weight.first = entry.first;
			pair_with_weight.second = std::pair<unsigned int, unsigned int>(j, i);
			resemblance_pairs_left_right.push_back(pair_with_weight);
			used_right_r[i] = used_left_r[j] = true;  // Mark these objects as used
		}
	}

	std::sort(resemblance_pairs_left_right.begin(), resemblance_pairs_left_right.end());
	std::vector<bool> used;
	if (N_left > N_right)
	{
		used.resize(N_left, false);
	}
	else
	{
		used.resize(N_right, false);
	}
	std::vector<bool> used_right(N_right, false);
	std::vector<bool> used_left(N_left, false);
	for (const auto& entry : resemblance_pairs_left_right)
	{
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_left[i] )
		{
			used_left[i] = used_right[j] = true;
			unsigned int p_left = trilateral_desc_left[i].p1;;
			unsigned int p_right = trilateral_desc_right[j].p1;;
			resemblance_pairs.push_back(std::pair<unsigned int, unsigned int >(p_left, p_right));
		}
	}
	//lastly check best 
	return resemblance_pairs;

}
static std::vector<std::pair<unsigned int, unsigned int>>  trilateral_area_histogram(TrilateralMesh* m, std::vector<TrilateralDescriptor>& trilateral_desc_left,
	std::vector<TrilateralDescriptor>& trilateral_desc_right)
{
	//dvorak externsion
	std::vector<unsigned int> left_indices;
	std::vector<unsigned int> right_indices;
	for (size_t i = 0; i < trilateral_desc_left.size(); i++)
	{
		left_indices.push_back(trilateral_desc_left[i].p1);
	}
	for (size_t i = 0; i < trilateral_desc_right.size(); i++)
	{
		right_indices.push_back(trilateral_desc_right[i].p1);
	}
	std::vector<DvorakPairs> left_dvorak = dvorak_extraction_of_significant_points(m, left_indices);
	std::vector<DvorakPairs> right_dvorak = dvorak_extraction_of_significant_points(m, right_indices);
	std::vector<DvorakPairs> dvorak_all = dvorak_extraction_of_significant_points(m, right_indices);
	dvorak_all.insert(dvorak_all.begin(), left_dvorak.begin() , left_dvorak.end());
	dvorak_all.insert(dvorak_all.begin(), right_dvorak.begin() , right_dvorak.end());


	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	std::vector<std::pair< float, std::pair<unsigned int, unsigned int>>> resemblance_pairs_left_right;
	//compare each
	std::vector<std::pair<float, std::pair<int, int>>> compareResults_left;

	int N_left = trilateral_desc_left.size();
	int N_right = trilateral_desc_right.size();
	int histogram_size = trilateral_desc_left[0].histogram.size();
	for (size_t i = 0; i < N_left; i++)
	{
		int minimum_index = -1;
		float minimum_value = INFINITY;
		float hks_i = m->normalized_heat_kernel_signature[trilateral_desc_left[i].p1];
		std::vector<float> histogram_vec_i = trilateral_desc_left[i].histogram.histogram;
		/*for (size_t j = 0; j < histogram_size; j++)
		{
			histogram_vec_i.push_back(hks_i);
		}*/
		Eigen::VectorXd histogram_i = stdVectorToEigenVectorXd(trilateral_desc_left[i].histogram.histogram);
		for (size_t j = 0; j < N_right; j++)
		{
			float hks_j = m->normalized_heat_kernel_signature[trilateral_desc_right[j].p1];
			std::vector<float> histogram_vec_j = trilateral_desc_right[j].histogram.histogram;
			
			//check curvature
			//bool is_curvature  = dvorak_curvature_similarity_criterion(dvorak_all , 0.5, i , j + left_dvorak.size() );
			//bool is_angle  = dvorak_normal_angle_criterion(m ,dvorak_all , i , j + left_dvorak.size() , 0.985  );
			bool is_hks_dif = std::abs( (hks_i - hks_j) ) < 0.015;
			for (size_t k = 0; k < histogram_size; k++)
			{
				histogram_vec_j.push_back(hks_j);
			}
			Eigen::VectorXd histogram_j = stdVectorToEigenVectorXd(trilateral_desc_right[j].histogram.histogram);
			float resemblance = (histogram_i - histogram_j).norm() ; // L2 norm
			//float resemblance_c = histogram_i.dot(histogram_j) / (histogram_i.norm() * histogram_j.norm()); // cosine similarity 
			if (is_hks_dif)
			{
				compareResults_left.push_back({ resemblance, {i,j} });
			}
		}
	}

	std::vector<std::pair<int, int>> selectedPairs;
	float skeleton_resemblance_error = 0;
	std::vector<bool> used_left_l(N_left, false);
	std::vector<bool> used_right_l(N_right, false);
	std::sort(compareResults_left.begin(), compareResults_left.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults_left) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_left_l[i]) {
			std::pair<float, std::pair<unsigned int, unsigned int>> pair_with_weight;
			pair_with_weight.first = entry.first;
			pair_with_weight.second = std::pair<unsigned int, unsigned int>(i, j);
			resemblance_pairs_left_right.push_back(pair_with_weight);
			used_left_l[i] = used_right_l[j] = true;  // Mark these objects as used
		}
	}
	std::vector<std::pair<float, std::pair<int, int>>> compareResults_right;
	//do it for right to let also
	for (size_t i = 0; i < N_right; i++)
	{
		int minimum_index = -1;
		float minimum_value = INFINITY;
		float hks_i = m->normalized_heat_kernel_signature[trilateral_desc_right[i].p1];
		std::vector<float> histogram_vec_i = trilateral_desc_right[i].histogram.histogram;
		/*for (size_t j = 0; j < histogram_size; j++)
		{
			histogram_vec_i.push_back(hks_i);
		}*/
		Eigen::VectorXd histogram_i = stdVectorToEigenVectorXd(trilateral_desc_right[i].histogram.histogram);
		for (size_t j = 0; j < N_left; j++)
		{
			float hks_j = m->normalized_heat_kernel_signature[trilateral_desc_left[j].p1];
			std::vector<float> histogram_vec_j = trilateral_desc_left[j].histogram.histogram;
			//check curvature
			//bool is_curvature  = dvorak_curvature_similarity_criterion(dvorak_all , 0.5, i , j + left_dvorak.size() );
			//bool is_angle  = dvorak_normal_angle_criterion(m ,dvorak_all , i , j + left_dvorak.size() , 0.985  );
			bool is_hks_dif = std::abs((hks_i - hks_j)) < 0.03;
			/*for (size_t k = 0; k < histogram_size; k++)
			{
				histogram_vec_j.push_back(hks_j);
			}*/
			Eigen::VectorXd histogram_j = stdVectorToEigenVectorXd(trilateral_desc_left[j].histogram.histogram);
			float resemblance = (histogram_i - histogram_j).norm(); // L2 norm
			//float resemblance_c = histogram_i.dot(histogram_j) / (histogram_i.norm() * histogram_j.norm()); // cosine similarity 
			if (is_hks_dif)
			{
				compareResults_right.push_back({ resemblance, {i,j} });
			}
		}
	}

	std::vector<std::pair<int, int>> selectedPairs_r;
	std::vector<bool> used_left_r(N_left, false);
	std::vector<bool> used_right_r(N_right, false);
	std::sort(compareResults_right.begin(), compareResults_right.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults_right) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_right_r[i]  ) {
			std::pair<float, std::pair<unsigned int, unsigned int >> pair_with_weight;
			pair_with_weight.first = entry.first;
			pair_with_weight.second = std::pair<unsigned int, unsigned int>(j, i);
			resemblance_pairs_left_right.push_back(pair_with_weight);
			used_right_r[i] = used_left_r[j] = true;  // Mark these objects as used
		}
	}

	std::sort(resemblance_pairs_left_right.begin(), resemblance_pairs_left_right.end());
	std::vector<bool> used;
	if (N_left > N_right)
	{
		used.resize(N_left, false);
	}
	else
	{
		used.resize(N_right, false);
	}
	std::vector<bool> used_right(N_right, false);
	std::vector<bool> used_left(N_left, false);
	for (const auto& entry : resemblance_pairs_left_right)
	{
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_left[i] || !used_right[j])
		{
			used_left[i] = used_right[j] = true;
			unsigned int p_left = trilateral_desc_left[i].p1;;
			unsigned int p_right = trilateral_desc_right[j].p1;;
			resemblance_pairs.push_back(std::pair<unsigned int, unsigned int >(p_left, p_right));
		}
	}
	//lastly check best 
	return resemblance_pairs;
}
static std::vector<std::pair<unsigned int, unsigned int>>  trilateral_unique_pairing_histogram(std::vector<TrilateralDescriptor>& trilateral_desc_left,
	std::vector<TrilateralDescriptor>& trilateral_desc_right )
{
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	//compare each
	std::vector<std::pair<float, std::pair<int, int>>> compareResults;

	int N_left = trilateral_desc_left.size();
	int N_right = trilateral_desc_right.size();
	for (size_t i = 0; i < N_left; i++)
	{
		Eigen::VectorXd histogram_i = stdVectorToEigenVectorXd(trilateral_desc_left[i].histogram.histogram);
		int minimum_index = -1;
		float minimum_value = INFINITY;
		for (size_t j = 0; j < N_right; j++)
		{
			//check dvorak
			Eigen::VectorXd histogram_j = stdVectorToEigenVectorXd(trilateral_desc_right[j].histogram.histogram);
			float resemblance = (histogram_i - histogram_j).norm(); // L2 norm
			//float resemblance = histogram_i.dot(histogram_j) / (histogram_i.norm() * histogram_j.norm()); // cosine similarity 
			compareResults.push_back({ resemblance, {i,j} });
		}
	}

	std::vector<std::pair<int, int>> selectedPairs;
	float skeleton_resemblance_error = 0;
	std::vector<bool> used_left(N_left, false);
	std::vector<bool> used_right(N_right, false);
	std::sort(compareResults.begin(), compareResults.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used_left[i] && !used_right[j]) {
			resemblance_pairs.push_back({ trilateral_desc_left[i].p1 , trilateral_desc_right[j].p1 });
			used_left[i] = used_right[j] = true;  // Mark these objects as used
		}
	}
	return resemblance_pairs;

}
void trilateral_FPS_histogram_matching(MeshFactory& mesh_fac, const int& selected_index, int sample_no , int division_no , bool recordTxt)
{
	TrilateralMesh* mesh = &mesh_fac.mesh_vec[selected_index];
	int N = mesh->vertices.size();
	std::vector<Histogram> trilateral_histograms;
	std::vector<unsigned int> sampled_points = furthest_point_sampling(mesh , sample_no , true );
	
	// use dijkstra to get each beest neihbours
	std::vector<TrilateralDescriptor> trilateral_desc = get_trilateral_points_using_closest_pairs(mesh_fac, selected_index, sampled_points);
	std::vector<int> global_is_visited( N , OUTSIDE); 
	for (size_t i = 0; i < sample_no; i++)
	{
		ROI_trilateral(mesh, trilateral_desc[i], division_no , false );
		Histogram histogram = Histogram_triangle_area(mesh,trilateral_desc[i], division_no);
		trilateral_histograms.push_back(histogram);
	}
	for (size_t i = 0; i < mesh->colors.size(); i++)
	{
		if (global_is_visited[i] == EDGE)
		{
			mesh->colors[i] = glm::vec3(255.0f , 0.0f ,0.0f);
		}
		if (global_is_visited[i] == INSIDE)
		{
			mesh->colors[i] = glm::vec3(0.0f, 255.0f, 0.0f);

		}
	}
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	resemblance_pairs = trilateral_unique_pairing_histogram(trilateral_desc, sample_no, trilateral_histograms);
	//compare each
	/*for (size_t i = 0; i < sample_no; i++)
	{
		Eigen::VectorXd histogram_i = stdVectorToEigenVectorXd(trilateral_histograms[i]);
		int minimum_index = -1;
		float minimum_value = INFINITY;
		for (size_t j = 0; j < sample_no; j++)
		{
			if (i == j)
			{
				continue;
			}
			Eigen::VectorXd histogram_j = stdVectorToEigenVectorXd(trilateral_histograms[j]);
			if ( (histogram_i - histogram_j).norm() < minimum_value)
			{
				minimum_value = (histogram_i - histogram_j).norm();
				minimum_index = j;
			}

		}
		resemblance_pairs.push_back(std::pair<int,int>(trilateral_desc[i].p1 , trilateral_desc[minimum_index].p1));
	}*/

	//buffer
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		mesh->colors[resemblance_pairs[i].first].r = 255;
		mesh->colors[resemblance_pairs[i].first].g = 0;
		mesh->colors[resemblance_pairs[i].first].b = 0;

		mesh->colors[resemblance_pairs[i].second].r = 0;
		mesh->colors[resemblance_pairs[i].second].g = 0;
		mesh->colors[resemblance_pairs[i].second].b = 255;
	}

	mesh->calculated_symmetry_pairs = resemblance_pairs;

	if (recordTxt)
	{
		Metric_write_to_file(mesh, "../../Results/Base_Trilateral_FPS_histogram.txt");
	}

}


// exact same with function above except do the histogram in MDS 
void trilateral_FPS_histogram_matching_w_principal_comp(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no)
{
	TrilateralMesh* mesh = &mesh_fac.mesh_vec[selected_index];
	std::vector<std::vector<float>> trilateral_histograms;
	std::vector<unsigned int> sampled_points = furthest_point_sampling(mesh, sample_no, true);
	int N = mesh->vertices.size();

	// use dijkstra to get each beest neihbours
	std::vector<TrilateralDescriptor> trilateral_desc = get_trilateral_points_using_closest_pairs(mesh_fac, selected_index, sampled_points);
	std::vector<int> global_is_visited;
	std::vector<std::vector<int>> is_visited_list(mesh->vertices.size());
	std::vector<std::vector<int>> points_inside_list(mesh->vertices.size());
	for (size_t i = 0; i < mesh->vertices.size(); i++)
	{
		global_is_visited.push_back(OUTSIDE);
	}
	for (size_t i = 0; i < sample_no; i++)
	{
		if (i == 55)
		{
			int a = 1; 
		}
		std::vector<float> histogram  = computePrincipalCurvatures(mesh_fac, selected_index, is_visited_list[i],
		division_no);
		trilateral_histograms.push_back(histogram);

	}
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	//compare each
	for (size_t i = 0; i < sample_no; i++)
	{
		float smallest_dist = INFINITY;
		int smallest_index = -1;
		for (size_t j = 0; j < sample_no; j++)
		{
			if (j == i)
			{
				continue;
			}
			float dist = 0; //= fabs(chi_squre_distance(trilateral_histograms[i], trilateral_histograms[j]));
			for (size_t k = 0; k < division_no; k++)
			{
				dist += powf(trilateral_histograms[i][k] - trilateral_histograms[j][k], 2);
			}
			if (smallest_dist > dist)
			{
				smallest_dist = dist;
				smallest_index = j;
			}
		}
		resemblance_pairs.push_back(std::pair<int, int>(sampled_points[i], sampled_points[smallest_index]));
	}

	//buffer
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		mesh->colors[resemblance_pairs[i].first].r = 255;
		mesh->colors[resemblance_pairs[i].first].g = 0;
		mesh->colors[resemblance_pairs[i].first].b = 0;

		mesh->colors[resemblance_pairs[i].second].r = 0;
		mesh->colors[resemblance_pairs[i].second].g = 0;
		mesh->colors[resemblance_pairs[i].second].b = 255;
	}

	mesh->calculated_symmetry_pairs = resemblance_pairs;

}
static std::vector<int> getPointsInside(MeshFactory& meshFac, int selectedIndex, std::vector<int>& is_visited)
{
	TrilateralMesh* m = &meshFac.mesh_vec[selectedIndex];
	int numPoints = is_visited.size();
	std::vector<int> vertices_inside; 
	for (size_t i = 0; i < numPoints; i++)
	{
		if (is_visited[i] != OUTSIDE)
		{
			vertices_inside.push_back(i);
		}
	}
	return vertices_inside;
}
// A sampler of  Useful  page 128 
static std::vector<float> computePrincipalCurvatures(MeshFactory& meshFac, int selectedIndex, std::vector<int>& is_visited,
	int division_no) 
{
	
	TrilateralMesh* m = &meshFac.mesh_vec[selectedIndex];
	std::vector<float> principalCurvatures1;
	std::vector<float> principalCurvatures2;
	std::vector<float> shapeIndices;
	int numVertices = m->vertices.size();
	int numTriangles = m->triangles.size();
	principalCurvatures1.resize(numVertices, 0.0);
	principalCurvatures2.resize(numVertices, 0.0);
	std::vector<int> vertices_inside = getPointsInside(meshFac, selectedIndex, is_visited);

	//get numbe of points inside
	int numTrianglesInside = 0;
	for (int i = 0; i < numTriangles; i += 3) 
	{
		int p1_index = m->triangles[i];
		int p2_index = m->triangles[i + 1];
		int p3_index = m->triangles[i + 2];
		if( is_visited[p1_index] != OUTSIDE && is_visited[p2_index] != OUTSIDE
		&& is_visited[p3_index] != OUTSIDE)
		{
			numTrianglesInside++;
		}
	}
	// For each vertex, we need to estimate the curvature tensor
	for (int i = 0; i < vertices_inside.size(); ++i)
	{
		Eigen::Matrix3d curvatureTensor = Eigen::Matrix3d::Zero();

		/*for (int j = 0; j < numTriangles; j += 3)
		{
			if (m->triangles[j] == i || m->triangles[j + 1] == i || m->triangles[j + 2] == i) 
			{
				glm::vec3 v1 = m->vertices[m->triangles[j]];
				glm::vec3 v2 = m->vertices[m->triangles[j + 1]];
				glm::vec3 v3 = m->vertices[m->triangles[j + 2]];

				glm::vec3 cross_res = glm::normalize(glm::cross(v3 - v1, v2 - v1));
				Eigen::Vector3d normal(cross_res.x,cross_res.y,cross_res.z);
				vertexNormal += normal;

				// Compute the edge vectors
				Eigen::Vector3d e1 = Eigen::Vector3d(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z).normalized();
				Eigen::Vector3d e2 = Eigen::Vector3d(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z).normalized();

				// Compute the curvature tensor for the current triangle
				curvatureTensor += (normal * normal.transpose());
			}
		}*/
		curvatureTensor  = computeCurvatureTensor(meshFac.mesh_vec[selectedIndex], vertices_inside[i]);

		// Diagonalize the curvature tensor to get the principal curvatures
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(curvatureTensor);
		Eigen::Vector3d eigenvalues = solver.eigenvalues();

		principalCurvatures1.push_back(eigenvalues(0));
		principalCurvatures2.push_back(eigenvalues(2)); // Smallest and largest eigenvalues correspond to principal curvatures
		float shape_index = -2 * M_PI * atan( (eigenvalues(2) - eigenvalues(0)) / (eigenvalues(0) + eigenvalues(2)));
		if (_isnanf(shape_index))
		{
			int d = 1;
		}
		shapeIndices.push_back(shape_index);
	}
	//generate histogram
	// 1 - get average point
	glm::vec3 avg_point(0.0f, 0.0f, 0.0f);
	for (size_t i = 0; i < vertices_inside.size(); i++)
	{
		avg_point += m->vertices[vertices_inside[i]];
	}
	avg_point /= vertices_inside.size();

	//get the point nearest to avg_point
	int nearest_vertex = -1;
	float nearest_diff = INFINITY;
	for (size_t i = 0; i < vertices_inside.size(); i++)
	{
		float dist = glm::distance(avg_point, m->vertices[vertices_inside[i]]);
		if (nearest_diff > dist)
		{
			nearest_diff = dist; 
			nearest_vertex = i; 
		}
	}

	std::vector<float> point_distances(vertices_inside.size());
	for (size_t i = 0; i < vertices_inside.size(); i++)
	{
		point_distances[i] = glm::distance(m->vertices[vertices_inside[i]], avg_point );
	}
	//generate histogram

	// 1- get min and max
	float min = INFINITY;
	float max = -INFINITY;
	for (size_t i = 0; i < point_distances.size(); i++)
	{
		if (point_distances[i] > max)
		{
			max = point_distances[i];
		}
		if (point_distances[i] < min)
		{
			min = point_distances[i];
		}
	}
	float step = (max - min) / division_no;

	std::vector<float> histogram(division_no, 0);
	std::vector<int> histogram_sizes(division_no,0 );
	
	for (size_t i = 0; i < vertices_inside.size(); i++)
	{
		int index =  (point_distances[i] - min) / step;
		if (max == point_distances[i])
		{
			index--;
		}
		histogram[index] += shapeIndices[i] ;
		histogram_sizes[index]++;
	}

	//normalize histogram
	for (size_t i = 0; i < division_no; i++)
	{
		histogram[i] /= histogram_sizes[i];
	}
	return histogram; 
}

static Eigen::Matrix3d computeCurvatureTensor(const TrilateralMesh& mesh, int vertexIndex) 
{
	const glm::vec3& v = mesh.vertices[vertexIndex];
	Eigen::Matrix3d curvatureTensor = Eigen::Matrix3d::Zero();
	double area = mesh.areas[vertexIndex];

	for (int i = 0; i < mesh.adjacenies[vertexIndex].size(); i++) {
		glm::vec3 edgeVector = mesh.vertices[mesh.adjacenies[vertexIndex][i].first] - v;
		Eigen::Vector3d edgeVecEigen;
		edgeVecEigen(0) = edgeVector.x;
		edgeVecEigen(1) = edgeVector.y;
		edgeVecEigen(2) = edgeVector.z;
		edgeVecEigen.normalize();
		curvatureTensor += /*e.angle **/ edgeVecEigen * edgeVecEigen.transpose();
	}

	return curvatureTensor / area;
}

void trilateral_fuzzyGeodesic(MeshFactory& meshFac, int selectedIndex, int p1, int p2, int p3 , float fuzziness_sigma )
{
	FuzzyGeodesicList fuzzyLists[3];
	TrilateralMesh* m = &meshFac.mesh_vec[selectedIndex];
	fuzzyLists[0] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p2, fuzziness_sigma);
	fuzzyLists[1] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p2, p3, fuzziness_sigma);
	fuzzyLists[2] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p3, fuzziness_sigma);
   
	FuzzyGeodesic_FuzzyArea(m, fuzzyLists[0], true);
	FuzzyGeodesic_FuzzyArea(m, fuzzyLists[1], true);
	FuzzyGeodesic_FuzzyArea(m, fuzzyLists[2], true);
}

void trilateral_FPS_matching_w_fuzzy_geodesic(MeshFactory& mesh_fac, const int& selected_index, int sample_no , float fuzziness_sigma , bool recordTxt )
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];
	std::vector<unsigned int> sampled_points = furthest_point_sampling(m, sample_no, true);
	int N = m->vertices.size();
	int sample_size = sampled_points.size();

	// generate closest  trialterals
	std::vector<TrilateralDescriptor> trilateral_desc = get_trilateral_points_using_closest_pairs(mesh_fac, selected_index, sampled_points);

	std::vector<Eigen::Vector3f> area_vectors(sample_size);
	for (size_t i = 0; i < sample_size; i++)
	{
		int p1 = trilateral_desc[i].p1;
		int p2 = trilateral_desc[i].p2;
		int p3 = trilateral_desc[i].p3;
		FuzzyGeodesicList fuzzyLists[3];
		
		fuzzyLists[0] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p2, fuzziness_sigma);
		fuzzyLists[1] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p2, p3, fuzziness_sigma);
		fuzzyLists[2] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p3, fuzziness_sigma);

		area_vectors[i][0] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[0], true);
		area_vectors[i][1] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[1], true);
		area_vectors[i][2] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[2], true);

	}

	
	//compare each other 
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	for (size_t i = 0; i < sample_size; i++)
	{
		int min_index = -1;
		float min_float = INFINITY;
		for (size_t j = 0; j < sample_size; j++)
		{
			if (i == j)
			{
				continue;
			}
			Eigen::Vector3f dif_vec = area_vectors[i] - area_vectors[j];
			float dif = dif_vec.norm();
			if (min_float > dif)
			{
				min_float = dif; 
				min_index = j;
			}
		}

		resemblance_pairs.push_back(std::pair<int, int>(sampled_points[i], sampled_points[min_index]));

	}

	//buffer
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 0;
		m->colors[resemblance_pairs[i].second].b = 255;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;

	if (recordTxt)
	{
		Metric_write_to_file(m, "../../Results/Trilateral_w_FuzzyGeodesic.txt");
	}

}

void trilateral_w_skeleton_endpoints(MeshFactory& mesh_fac, const int& selected_index,
	float fuzziness_sigma, Skeleton& skeleton, bool recordTxt)
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point; 
	// 1 - get end points from skeleton
	std::vector<unsigned int> mesh_endpoints;
	skeleton_calculate_closest_mesh_points(skeleton, m, mesh_endpoints);
	int sample_no = mesh_endpoints.size();

	// generate closest  trialterals
	std::vector<TrilateralDescriptor> trilateral_desc = get_trilateral_points_using_closest_pairs(mesh_fac, selected_index, mesh_endpoints);
	std::vector<Eigen::Vector3f> area_vectors(sample_no);
	for (size_t i = 0; i < sample_no; i++)
	{
		int p1 = trilateral_desc[i].p1;
		int p2 = trilateral_desc[i].p2;
		int p3 = trilateral_desc[i].p3;
		FuzzyGeodesicList fuzzyLists[3];

		fuzzyLists[0] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p2, fuzziness_sigma);
		fuzzyLists[1] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p2, p3, fuzziness_sigma);
		fuzzyLists[2] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p3, fuzziness_sigma);

		area_vectors[i][0] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[0], true);
		area_vectors[i][1] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[1], true);
		area_vectors[i][2] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[2], true);

	}

	mesh_mid_point_index = skeleton_calculate_closest_mesh_point( skeleton,m , skeleton.mid_point_index);
	std::vector<float> distances_from_mesh_mid_point = Geodesic_dijkstra(*m, mesh_mid_point_index);
	mesh_mid_point = m->vertices[mesh_mid_point_index];
	//compare each other 
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	for (size_t i = 0; i < sample_no; i++)
	{
		int min_index = -1;
		float min_float = INFINITY;
		std::vector<float> dist_i = Geodesic_dijkstra(*m, mesh_endpoints[i]);
		for (size_t j = 0; j < sample_no; j++)
		{
			if (i == j)
			{
				continue;
			}
			// !! also check the distance between i and checkpoint and i and j
			/*float dist_i_mid = distances_from_mesh_mid_point[mesh_endpoints[i]];
			float dist_i_j = dist_i[mesh_endpoints[j]];
			if (dist_i_mid *15.0f/10.0f > dist_i_j)
			{
				continue; 
			}*/
			// also check if the distance
			//Eigen::Vector3f dif_vec = area_vectors[i] - area_vectors[j];
			float dif = permutation_return_smallest_dif( area_vectors[i], area_vectors[j], 3);
			if (min_float > dif )
			{
				min_float = dif;
				min_index = j;
			}
		}
		resemblance_pairs.push_back(std::pair<int, int>(mesh_endpoints[i], mesh_endpoints[min_index]));
	}
	//buffer
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 0;
		m->colors[resemblance_pairs[i].second].b = 255;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;
	if (recordTxt)
	{
		Metric_write_to_file(m, "../../Results/Trilateral_w_FuzzyGeodesic.txt");
	}

}
void trilateral_FPS_symmetry_using_skel_endpoints(MeshFactory& mesh_fac, const int& selected_index, int sample_no,
float fuzziness_sigma,Skeleton& skeleton,bool recordTxt)
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];
	int N = m->vertices.size();
	// 1 - get end points from skeleton
	std::vector<unsigned int> mesh_endpoints;
	skeleton_calculate_closest_mesh_points(skeleton, m, mesh_endpoints);
	// 2 - FPS 
	std::vector<unsigned int> sampled_points = furthest_point_sampling(m, sample_no, true);
	// 3 - generate closest  trialterals
	std::vector<TrilateralDescriptor> trilateral_desc = get_trilateral_points_using_closest_pairs(mesh_fac, selected_index, sampled_points);

	std::vector<Eigen::Vector3f> area_vectors(sample_no);
	for (size_t i = 0; i < sample_no; i++)
	{
		int p1 = trilateral_desc[i].p1;
		int p2 = trilateral_desc[i].p2;
		int p3 = trilateral_desc[i].p3;
		FuzzyGeodesicList fuzzyLists[3];

		fuzzyLists[0] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p2, fuzziness_sigma);
		fuzzyLists[1] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p2, p3, fuzziness_sigma);
		fuzzyLists[2] = FuzzyGeodesic_calculateFuzzyGedoesic(m, p1, p3, fuzziness_sigma);

		area_vectors[i][0] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[0], true);
		area_vectors[i][1] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[1], true);
		area_vectors[i][2] = FuzzyGeodesic_FuzzyArea(m, fuzzyLists[2], true);

	}
}

void trilateral_point_matching_with_skeleton_endpoints(TrilateralMesh* m , Skeleton& skeleton)
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get end points from skeleton
	std::vector<unsigned int> mesh_endpoints;
	skeleton_calculate_closest_mesh_points(skeleton, m, mesh_endpoints);
	int sample_no = mesh_endpoints.size();
	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();

	Plane plane; 
	/*if (!dom_sym_read_plane(m, plane))
	{
		plane = generate_dominant_symmetry_plane(m, 0.1);
	}*/

	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < sample_no; i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[mesh_endpoints[i]]) > 0  ) 
		{
			right_skeleton_indices.push_back(mesh_endpoints[i]);
		}
		else
		{
			left_skeleton_indices.push_back(mesh_endpoints[i]);
		}
	}

	//std::vector<TrilateralDescriptor> positive_mesh_trilateral_descriptor = get_trilateral_points_using_closest_pairs(&m, right_skeleton_indices);
	//std::vector<TrilateralDescriptor> negative_mesh_trilateral_descriptor = get_trilateral_points_using_closest_pairs(&m, left_skeleton_indices);

	//color
	/*for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		m->colors[right_skeleton_indices[i]] = glm::vec3(0, 0, 255);
	}
	for (size_t i = 0; i < left_skeleton_indices.size(); i++)
	{
		m->colors[left_skeleton_indices[i]] = glm::vec3(0, 255, 0);
	}

	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<int> global_is_visited(m->vertices.size(), OUTSIDE);
		std::vector<int> is_visited = trilateral_ROI(&m, positive_mesh_trilateral_descriptor[i].p1,
			positive_mesh_trilateral_descriptor[i].p2, positive_mesh_trilateral_descriptor[i].p3, 3, false);
		std::vector<float> histogram = histogramROi(mesh_fac, (int&)selected_index, positive_mesh_trilateral_descriptor[i].p1,
			positive_mesh_trilateral_descriptor[i].p2, positive_mesh_trilateral_descriptor[i].p3, 3, is_visited, global_is_visited);
		positive_mesh_trilateral_descriptor[i].histogram = histogram;
	}
	for (size_t i = 0; i < left_skeleton_indices.size(); i++)
	{
		std::vector<int> global_is_visited(m->vertices.size(), OUTSIDE);
		std::vector<int> is_visited = trilateral_ROI(&m, negative_mesh_trilateral_descriptor[i].p1,
			negative_mesh_trilateral_descriptor[i].p2, negative_mesh_trilateral_descriptor[i].p3, 3, false);
		std::vector<float> histogram = histogramROi(mesh_fac, (int&)selected_index, negative_mesh_trilateral_descriptor[i].p1,
			negative_mesh_trilateral_descriptor[i].p2, negative_mesh_trilateral_descriptor[i].p3, 3, is_visited, global_is_visited);
		negative_mesh_trilateral_descriptor[i].histogram = histogram;
	}
	// write a function for comparing two descriptor
	//irrelevant constants 
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = trilateral_unique_pairing_histogram(positive_mesh_trilateral_descriptor,
		negative_mesh_trilateral_descriptor);*/ 
	
		//point_match_trilateral_weights(&L_MDS_mesh, positive_mesh_trilateral_descriptor, negative_mesh_trilateral_descriptor, const1, const2, const3);

	//check their gaussian curvature
	/*std::vector<std::pair<unsigned int, unsigned int>> new_resemblance_pairs;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		float curv1 = gaussian_curvature(&m, resemblance_pairs[i].first);
		float curv2 = gaussian_curvature(&m, resemblance_pairs[i].second);
		if ( curv1/curv2 > 0.5 && curv2/curv1  < 2 )
		{
			new_resemblance_pairs.push_back(resemblance_pairs[i]);
		}

	}
	resemblance_pairs = new_resemblance_pairs;
	new_resemblance_pairs.clear();
	//normal angle criterion
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		std::vector<DvorakPairs> pair;
		DvorakPairs d1 ,d2;
		d1.p_index = resemblance_pairs[i].first;
		d2.p_index = resemblance_pairs[i].second;
		pair.push_back(d1);
		pair.push_back(d2);
		bool is_normal_angle = dvorak_normal_angle_criterion(&m, pair, 0, 1, 0.985);
		if (is_normal_angle)
		{
			new_resemblance_pairs.push_back(resemblance_pairs[i]);
		}
	}
	resemblance_pairs = new_resemblance_pairs;*/

	//forge it into two list
	/*std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 0;
		m->colors[resemblance_pairs[i].second].b = 255;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;

	Metric_write_to_file(m, "../../Results/Trilateral_W_LMDS_AND_DOMINANTPLANE.txt"); */

}
void trilateral_point_matching_with_skeleton_endpoints_w_HKS(TrilateralMesh* m, Skeleton& skeleton , std::vector<TrilateralDescriptor>&desc_left,
std::vector<TrilateralDescriptor>& desc_right , Plane& plane  )
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get end points from skeleton
	std::vector<unsigned int> mesh_endpoints;
	skeleton_calculate_closest_mesh_points(skeleton, m, mesh_endpoints);
	int sample_no = mesh_endpoints.size();
	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	//if (!dom_sym_read_plane(m, plane, plane_path))
	{
		plane = generate_dominant_symmetry_plane(m, 2);
	}

	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < sample_no; i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[mesh_endpoints[i]]) > 0)
		{
			right_skeleton_indices.push_back(mesh_endpoints[i]);
		}
		else
		{
			left_skeleton_indices.push_back(mesh_endpoints[i]);
		}
	}

	desc_left  = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	desc_right = get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);

	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = trilateral_hks_histogram(m, desc_left,
		desc_right , 5, 3 );//point_match_trilateral_weights(&L_MDS_mesh, positive_mesh_trilateral_descriptor, negative_mesh_trilateral_descriptor, const1, const2, const3);

	//check their gaussian curvature


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 0;
		m->colors[resemblance_pairs[i].second].b = 255;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;

	Metric_write_to_file(m, "../../Results/Trilateral_W_SKELETON_AND_HKS.txt");


}
void trilateral_point_matching_with_dvorak_endpoints(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_left,
	std::vector<TrilateralDescriptor>& desc_right, Plane& plane, int dvorak_enpoint_no , float convergence_ratio )
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<DvorakPairs> dvorak_pairs  = dvorak_extraction_of_significant_points(m, dvorak_enpoint_no);
	//sweep 
	dvorak_pairs = dvorak_distance_sweep(m, dvorak_pairs, 3.0f);
	Metric_set_gaussian(m, dvorak_enpoint_no, 3.0f);

	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	if (plane.isNull() )
	{
		plane = generate_dominant_symmetry_plane(m, convergence_ratio);
	}
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[dvorak_pairs[i].p_index]) > 0)
		{
			right_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
		else
		{
			left_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
	}

	desc_left = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	desc_right = get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);


	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = trilateral_hks_histogram(m, desc_left,
		desc_right, 5, 3);//point_match_trilateral_weights(&L_MDS_mesh, positive_mesh_trilateral_descriptor, negative_mesh_trilateral_descriptor, const1, const2, const3);

	//check their gaussian curvature


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;

	std::string path = "../../Results/";
	path = path + m->file_name + " Trilateral_W_Gaussian_curvature.txt ";
	Metric_write_to_file(m, path);

}

//OT being optimal transform
void trilateral_point_matching_with_gaussian_endpoints_and_OT(TrilateralMesh* m,Skeleton& skeleton, std::vector<TrilateralDescriptor>& desc_left,
	std::vector<TrilateralDescriptor>& desc_right, Plane& plane, int dvorak_enpoint_no, float convergence_ratio)
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<DvorakPairs> dvorak_pairs = dvorak_extraction_of_significant_points(m, dvorak_enpoint_no);
	/*std::vector<unsigned int> end_points;
	skeleton_get_end_points(skeleton, end_points);
	for (size_t i = 0; i < end_points.size(); i++)
	{
		DvorakPairs p;
		p.p_index = end_points[i];
		dvorak_pairs.push_back(p);
	}*/
	//sweep 
	dvorak_pairs = dvorak_distance_sweep(m, dvorak_pairs, 5.0f);
	Metric_set_gaussian(m, dvorak_enpoint_no, 5.0f);

	

	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	if (plane.isNull())
	{
		plane = generate_dominant_symmetry_plane(m, convergence_ratio);
	}
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[dvorak_pairs[i].p_index]) > 0)
		{
			right_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
		else
		{
			left_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
	}

	std::vector<unsigned int > mesh_vertices;
	SkeletonTree skelTree = skeleton_generate_skeleton_tree(m, skeleton);
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		mesh_vertices.push_back(dvorak_pairs[i].p_index);
	}
	std::vector<TrilateralDescriptor> descriptors;
    TrilateralDescriptor_get_trilateral_with_closest_points(m, skeleton, mesh_vertices, skelTree, descriptors);

	for (size_t i = 0; i < descriptors.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[descriptors[i].p1]) > 0)
		{
			desc_left.push_back(descriptors[i]);
		}
		else
		{
			desc_right.push_back(descriptors[i]);
		}
	}


	//desc_left = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	//desc_right = get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);


	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all(m, desc_left, desc_right ,true , 12 , 2);
	//std::vector<std::vector<float>> optimal_transforms_neg_pos = VarianceMin_compare_all(m, desc_right, desc_left ,true , 20 , 3) ;
	std::vector<std::pair<float,std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	

	for (size_t i = 0; i < desc_left.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_right.size(); j++)
		{
			int dvoak_index_i;
			int dvoak_index_j;
			for (size_t k = 0; k < dvorak_pairs.size(); k++)
			{
				if (dvorak_pairs[k].p_index == desc_left[i].p1)
				{
					dvoak_index_i = k;
				}
				if (dvorak_pairs[k].p_index == desc_right[j].p1)
				{
					dvoak_index_j = k;
				}
			}
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[desc_left[i].p1] - m->normalized_heat_kernel_signature[desc_right[j].p1]);
			bool is_hks = hks_dif < 0.3; 
			bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, 0.5, dvoak_index_i, dvoak_index_j);
			std::cout << " curv " << is_curv << std::endl;
			bool is_norm = dvorak_normal_angle_criterion(m,dvorak_pairs,  dvoak_index_i, dvoak_index_j,0.985);
			std::cout << " normal " << is_norm << std::endl;
			if (is_curv && is_norm )
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = optimal_transforms_pos_neg[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}
		}
	}
	std::vector<bool> used(N, false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ desc_left[i].p1 , desc_right[j].p1 });
			used[i] = used[j] = true;  // Mark these objects as used
		}
	}

	/*for (size_t i = 0; i < desc_right.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_left.size(); j++)
		{
			if (optimal_transforms_neg_pos[i][j] < smallest)
			{
				smallest = optimal_transforms_neg_pos[i][j];
				index = j;
			}
		}
		std::pair<unsigned int, unsigned int> pair;
		pair.first = desc_right[i].p1;
		pair.second = desc_left[index].p1;
		resemblance_pairs.push_back(pair);
		std::cout << " left == " << index << " right == " << desc_left.size() + i <<  " " << optimal_transforms_neg_pos[i][index] << std::endl;
	}*/
	
	//check their gaussian curvature


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;



	std::string path = "../../Results/";
	path = path + m->file_name + " Trilateral_W_Gaussian_curvature_and_Optimal_Transform.txt ";
	Metric_write_to_file(m, path);

}

void trilateral_point_matching_with_skeleton_endpoints_and_OT(TrilateralMesh* m, Skeleton& skeleton, std::vector<TrilateralDescriptor>& desc_left,
	std::vector<TrilateralDescriptor>& desc_right, Plane& plane, int dvorak_enpoint_no, float convergence_ratio)
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	//sweep 
	std::vector<unsigned int> end_points; 
	std::vector<DvorakPairs> dvorak_pairs; 
	skeleton_get_end_points(skeleton, end_points);
	for (size_t i = 0; i < end_points.size(); i++)
	{
		float curv = gaussian_curvature(m, end_points[i]);
		DvorakPairs p;
		p.p_index = end_points[i];
		p.gaussian_curv = curv;
		dvorak_pairs.push_back(p);
	}
	dvorak_pairs = dvorak_distance_sweep(m, dvorak_pairs, 3.0f);
	Metric_set_gaussian(m, dvorak_enpoint_no, 3.0f);

	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	if (plane.isNull())
	{
		plane = generate_dominant_symmetry_plane(m, convergence_ratio);
	}
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[dvorak_pairs[i].p_index]) > 0)
		{
			right_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
		else
		{
			left_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
	}

	desc_left = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	desc_right = get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);


	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all(m, desc_left, desc_right, true, 5, 1);
	//std::vector<std::vector<float>> optimal_transforms_neg_pos = VarianceMin_compare_all(m, desc_right, desc_left ,true , 20 , 3) ;
	std::vector<std::pair<unsigned int, unsigned int >> resemblance_pairs;
	for (size_t i = 0; i < desc_left.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_right.size(); j++)
		{
			int dvoak_index_i;
			int dvoak_index_j;
			for (size_t k = 0; k < dvorak_pairs.size(); k++)
			{
				if (dvorak_pairs[k].p_index == desc_left[i].p1)
				{
					dvoak_index_i = k;
				}
				if (dvorak_pairs[k].p_index == desc_right[j].p1)
				{
					dvoak_index_j = k;
				}
			}
			if (optimal_transforms_pos_neg[i][j] < smallest)
			{
				bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, 0.5, dvoak_index_i, dvoak_index_j);
				std::cout << " curv " << is_curv << std::endl;
				bool is_norm = dvorak_normal_angle_criterion(m, dvorak_pairs, dvoak_index_i, dvoak_index_j, 0.9);
				std::cout << " normal " << is_norm << std::endl;
				if (is_curv && is_norm)
				{
					smallest = optimal_transforms_pos_neg[i][j];
					index = j;
				}
			}
		}
		std::pair<unsigned int, unsigned int> pair;
		pair.first = desc_left[i].p1;
		pair.second = desc_right[index].p1;
		resemblance_pairs.push_back(pair);
		std::cout << " left == " << i << " right == " << desc_left.size() + index << " " << optimal_transforms_pos_neg[i][index] << std::endl;
	}
	/*for (size_t i = 0; i < desc_right.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_left.size(); j++)
		{
			if (optimal_transforms_neg_pos[i][j] < smallest)
			{
				smallest = optimal_transforms_neg_pos[i][j];
				index = j;
			}
		}
		std::pair<unsigned int, unsigned int> pair;
		pair.first = desc_right[i].p1;
		pair.second = desc_left[index].p1;
		resemblance_pairs.push_back(pair);
		std::cout << " left == " << index << " right == " << desc_left.size() + i <<  " " << optimal_transforms_neg_pos[i][index] << std::endl;
	}*/

	//check their gaussian curvature


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;



	std::string path = "../../Results/";
	path = path + m->file_name + " Trilateral_W_skeleton_endpoints_and_Optimal_Transform.txt ";
	Metric_write_to_file(m, path);
}
//OT being optimal transform
void trilateral_point_matching_with_gaussian_endpoints_and_OT_w_CDF(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_left,
	std::vector<TrilateralDescriptor>& desc_right, Plane& plane, int dvorak_enpoint_no, float convergence_ratio)
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<DvorakPairs> dvorak_pairs = dvorak_extraction_of_significant_points(m, dvorak_enpoint_no);
	//sweep 
	dvorak_pairs = dvorak_distance_sweep(m, dvorak_pairs, 3.0f);
	Metric_set_gaussian(m, dvorak_enpoint_no, 3.0f);
	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	if (plane.isNull())
	{
		plane = generate_dominant_symmetry_plane(m, convergence_ratio);
	}
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[dvorak_pairs[i].p_index]) > 0)
		{
			right_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
		else
		{
			left_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
	}

	desc_left = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	desc_right = get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);


	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all_w_CDF(m, desc_left, desc_right,50);
	std::vector<std::vector<float>> optimal_transforms_neg_pos = VarianceMin_compare_all_w_CDF(m, desc_right, desc_left,50);
	std::vector<std::pair<unsigned int, unsigned int >> resemblance_pairs;
	for (size_t i = 0; i < desc_left.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_right.size(); j++)
		{
			if (optimal_transforms_pos_neg[i][j] < smallest)
			{
				smallest = optimal_transforms_pos_neg[i][j];
				index = j;
			}
		}
		std::pair<unsigned int, unsigned int> pair;
		pair.first = desc_left[i].p1;
		pair.second = desc_right[index].p1;
		resemblance_pairs.push_back(pair);
	}
	for (size_t i = 0; i < desc_right.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_left.size(); j++)
		{
			if (optimal_transforms_neg_pos[i][j] < smallest)
			{
				smallest = optimal_transforms_neg_pos[i][j];
				index = j;
			}
		}
		std::pair<unsigned int, unsigned int> pair;
		pair.first = desc_right[i].p1;
		pair.second = desc_left[index].p1;
		resemblance_pairs.push_back(pair);
	}
	//check their gaussian curvature


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;

	std::string path = "../../Results/";
	path = path + m->file_name + " Trilateral_W_Gaussian_curvature_and_Optimal_Transform_w_CDF.txt ";
	Metric_write_to_file(m, path);

}
void trilateral_point_matching_with_skeleton_endpoints_anchors(MeshFactory& mesh_fac, const int& selected_index, Skeleton& skeleton,
	std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& desc_neg , Plane& plane )
{
	TrilateralMesh* m = &mesh_fac.mesh_vec[selected_index];
	TrilateralMesh m_pos; 
	TrilateralMesh m_neg; 
	std::vector<int> indices_m_pos;
	std::vector<int> indices_m_neg;
	std::vector<unsigned int> skeleton_indices_m_pos;
	std::vector<unsigned int> skeleton_indices_m_neg;
	const int NO_OF_SAMPLE_FPS = 50;
	// 1 - use the plane for divide the mesh into two
	dom_sym_generate_two_separate_mesh_using_dominant_symmetry_plane(plane, m, &m_pos, &m_neg, &indices_m_pos, &indices_m_neg);

	// 2 - use fps on them
	std::vector<unsigned int> fps_indices_pos = furthest_point_sampling(&m_pos, NO_OF_SAMPLE_FPS, true);
	std::vector<unsigned int> fps_indices_neg = furthest_point_sampling(&m_neg, NO_OF_SAMPLE_FPS, true);
	
	//convert to the actual mesh
	for (size_t i = 0; i < fps_indices_pos.size(); i++)
	{
		fps_indices_pos[i] = indices_m_pos[fps_indices_pos[i]];
	}
	for (size_t i = 0; i < fps_indices_neg.size(); i++)
	{
		fps_indices_neg[i] = indices_m_neg[fps_indices_neg[i]];
	}
	// now they should be indices of the original mesh
	//get the indices of skeleton 
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		skeleton_indices_m_pos.push_back(desc_pos[i].p1);
	}
	for (size_t i = 0; i < desc_neg.size(); i++)
	{
		skeleton_indices_m_neg.push_back(desc_neg[i].p1);
	}
	// Do trilateral with fps points and skeleton points
	std::vector<TrilateralDescriptor> trilateral_pos_vec;
	std::vector<TrilateralDescriptor> trilateral_neg_vec;
	for (size_t i = 0; i < NO_OF_SAMPLE_FPS; i++)
	{
		TrilateralDescriptor pos = trilateral_get_trilateral_using_closest_pairs_with_skeleton_indices(m, fps_indices_pos[i], skeleton_indices_m_pos);
		TrilateralDescriptor neg =  trilateral_get_trilateral_using_closest_pairs_with_skeleton_indices(m, fps_indices_neg[i], skeleton_indices_m_neg);
	}
}
TrilateralDescriptor trilateral_get_trilateral_using_closest_pairs_with_skeleton_indices(TrilateralMesh* m, unsigned int point_index,
	std::vector<unsigned int>& skeleton_indices)
{
	std::vector<float> distances = Geodesic_dijkstra(*m, point_index);
	int N = m->vertices.size();
	int smallest_index = CoreType_return_smallest_k(distances, 1);
	int second_smallest_index = CoreType_return_smallest_k(distances, 2);

	TrilateralDescriptor  desc = TrilateralDescriptor_create(m, (int)point_index , smallest_index , second_smallest_index, false );
	return desc;
}



void trilateral_point_matching_with_skeleton_endpoints_SpinImage(TrilateralMesh* m, Skeleton& skeleton, std::vector<TrilateralDescriptor>& desc_left,
	std::vector<TrilateralDescriptor>& desc_right, Plane& plane)
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get end points from skeleton
	std::vector<unsigned int> mesh_endpoints; 
	skeleton_calculate_closest_mesh_points(skeleton, m, mesh_endpoints);
	int sample_no = mesh_endpoints.size();
	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	/*if (!dom_sym_read_plane(m, plane))
	{
		plane = generate_dominant_symmetry_plane(m, 2);
	}*/
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < sample_no; i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[mesh_endpoints[i]]) > 0)
		{
			right_skeleton_indices.push_back(mesh_endpoints[i]);
		}
		else
		{
			left_skeleton_indices.push_back(mesh_endpoints[i]);
		}
	}

	desc_left = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	desc_right = get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);

	std::vector<std::vector<int> > desc_left_visited_vertices;
	std::vector<std::vector<int> > desc_right_visited_vertices;
	
	std::vector<std::pair<unsigned int, unsigned int>> resemblances;
	for (size_t i = 0; i < desc_left.size(); i++)
	{
		Histogram2D spin_image_i = SpinImage_generate_spin_image(m, desc_left[i].p1, desc_left[i].visited_indices
			, 4);
		float smallest = INFINITY;
		int smallest_index = -1;
		for (size_t j = 0; j < desc_right.size(); j++)
		{
			Histogram2D spin_image_j = SpinImage_generate_spin_image(m, desc_right[j].p1, desc_right[i].visited_indices
				, 4);
			float resemblance = Histogram2D_L2Norm(spin_image_i, spin_image_j);
			if (smallest > resemblance)
			{
				smallest = resemblance;
				smallest_index = j; 
			}
		}
		int index_i = desc_left[i].p1;
		int index_j = desc_right[smallest_index].p1;
		resemblances.push_back(std::pair<unsigned int, unsigned int>( index_i ,index_j));
	}
	for (size_t i = 0; i < desc_right.size(); i++)
	{
		Histogram2D spin_image_i = SpinImage_generate_spin_image(m, desc_right[i].p1, desc_right[i].visited_indices
			, 4);
		float smallest = INFINITY;
		int smallest_index = -1;
		for (size_t j = 0; j < desc_left.size(); j++)
		{
			Histogram2D spin_image_j = SpinImage_generate_spin_image(m, desc_left[j].p1, desc_left[i].visited_indices
				, 4);
			float resemblance = Histogram2D_L2Norm(spin_image_i, spin_image_j);
			if (smallest > resemblance)
			{
				smallest = resemblance;
				smallest_index = j;
			}
		}
		int index_i = desc_right[i].p1;
		int index_j = desc_left[smallest_index].p1;
		resemblances.push_back(std::pair<int, int>(index_j, index_i));
	}

	m->calculated_symmetry_pairs = resemblances; 

	Metric_write_to_file(m, "../../Results/Trilateral_W_SKELETON_AND_HKS.txt");


}

void trilateral_display_trilateral_from_skeleton_endpoints(TrilateralMesh* m, std::vector<TrilateralDescriptor>& positive_desc
, std::vector<TrilateralDescriptor>& negative_desc,Skeleton& skeleton, Plane& plane )
{
	int N = m->vertices.size();
	int mesh_mid_point_index = -1;
	glm::vec3 mesh_mid_point;
	// 1 - get end points from skeleton
	std::vector<unsigned int> mesh_endpoints; 
	skeleton_calculate_closest_mesh_points(skeleton, m, mesh_endpoints);
	int sample_no = mesh_endpoints.size();
	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	if( plane.isNull() )
	{
		plane = generate_dominant_symmetry_plane(m, 2);
	}
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < sample_no; i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[mesh_endpoints[i]]) > 0)
		{
			right_skeleton_indices.push_back(mesh_endpoints[i]);
		}
		else
		{
			left_skeleton_indices.push_back(mesh_endpoints[i]);
		}
	}

	positive_desc = get_trilateral_points_using_closest_pairs(m, right_skeleton_indices);
	negative_desc= get_trilateral_points_using_closest_pairs(m, left_skeleton_indices);

	m->update_raylib_mesh();
}

void trilateral_sampled_point_matching_with_descriptors(TrilateralMesh* m, std::vector<TrilateralDescriptor>& positive_desc
	, std::vector<TrilateralDescriptor>& negative_desc, Plane& plane, std::vector<unsigned int>& fps_points , float sweep_distance)
{
	// 1 - divide fps points into two
	std::vector<unsigned int> positive_fps_points;
	std::vector<unsigned int> negative_fps_points;
	for (size_t i = 0; i < fps_points.size(); i++)
	{
		int index = fps_points[i];
		if (get_point_status_from_plane(&plane, &m->vertices[index]) > 0)
		{
			positive_fps_points.push_back(index);
		}
		else
		{
			negative_fps_points.push_back(index);
		}

	}
	 // 1 - check if  fps points are closer than sweep_distances
	std::vector<unsigned int> vertices_deleted(positive_fps_points.size() , 0);
	for (size_t i = 0; i < positive_desc.size(); i++)
	{
		int desc_index = positive_desc[i].p1;
		glm::vec3 pi = m->vertices[desc_index];
		//do a geodesic 
		std::vector<float> distances = Geodesic_dijkstra(*m, desc_index);
		for (size_t j = 0; j < positive_fps_points.size(); j++)
		{
			int fps_index = positive_fps_points[j];
			float distance_i_j = distances[fps_index];
			if (distance_i_j < sweep_distance)
			{
				vertices_deleted[j] = 1;
			}
		}
	}
	//delete
	std::vector<unsigned int> temp_vec;
	for (size_t i = 0; i < positive_fps_points.size(); i++)
	{
		if (vertices_deleted[i] == 0)
		{
			temp_vec.push_back(positive_fps_points[i]);
		}
	}
	positive_fps_points = temp_vec; 

	//negative 
	vertices_deleted.resize(negative_fps_points.size(), 0);
	for (size_t i = 0; i < negative_desc.size(); i++)
	{
		int desc_index = negative_desc[i].p1;
		glm::vec3 pi = m->vertices[desc_index];
		//do a geodesic 
		std::vector<float> distances = Geodesic_dijkstra(*m, desc_index);
		for (size_t j = 0; j < negative_fps_points.size(); j++)
		{
			int fps_index = negative_fps_points[j];
			float distance_i_j = distances[fps_index];
			if (distance_i_j < sweep_distance)
			{
				vertices_deleted[j] = 1;
			}
		}
	}
	//delete
	temp_vec.clear();
	for (size_t i = 0; i < negative_fps_points.size(); i++)
	{
		if (vertices_deleted[i] == 0)
		{
			temp_vec.push_back(negative_fps_points[i]);
		}
	}
	negative_fps_points = temp_vec;


	// 2 - 
}


std::vector<NLateralDescriptor> NlateralMap_point_matching_with_skeleton_endpoints_and_OT(TrilateralMesh* m, Skeleton& skeleton,
	Plane& plane, int dvorak_enpoint_no, float convergence_ratio, int N )
{
	int size = m->vertices.size();
	int mesh_mid_point_index = -1;
	std::vector<NLateralDescriptor> desc_neg;
	std::vector<NLateralDescriptor> desc_pos;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<DvorakPairs> dvorak_pairs = dvorak_extraction_of_significant_points(m, dvorak_enpoint_no);
	/*std::vector<unsigned int> end_points;
	skeleton_get_end_points(skeleton, end_points);
	for (size_t i = 0; i < end_points.size(); i++)
	{
		DvorakPairs p;
		p.p_index = end_points[i];
		dvorak_pairs.push_back(p);
	}*/
	//sweep 
	dvorak_pairs = dvorak_distance_sweep(m, dvorak_pairs, 3.0f);
	Metric_set_gaussian(m, dvorak_enpoint_no, 3.0f);
	Metric_set_N(N);

	std::vector<unsigned int > skel_mid_point_indices;
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		skel_mid_point_indices.push_back(dvorak_pairs[i].p_index);
	}
	skeleton_distance_to_midpoint(m, skeleton, skel_mid_point_indices);

	// 2- generate plane 
	 //calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		plane_center += m->vertices[i];
	}
	plane_center /= m->vertices.size();
	if (plane.isNull())
	{
		plane = generate_dominant_symmetry_plane(m, convergence_ratio);
	}
	//divide the end points 
	std::vector<unsigned int> left_skeleton_indices;
	std::vector<unsigned int> right_skeleton_indices;
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[dvorak_pairs[i].p_index]) > 0)
		{
			right_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
		else
		{
			left_skeleton_indices.push_back(dvorak_pairs[i].p_index);
		}
	}

	std::vector<unsigned int > mesh_vertices;
	SkeletonTree skelTree = skeleton_generate_skeleton_tree(m, skeleton);
	for (size_t i = 0; i < dvorak_pairs.size(); i++)
	{
		mesh_vertices.push_back(dvorak_pairs[i].p_index);
	}
	std::vector<NLateralDescriptor> descriptors;
	descriptors = NLateral_generate_closest_points(m, mesh_vertices, N , 10 );
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		std::vector<unsigned int> skel_dist_vec = { descriptors[i].indices[0] };
		std::vector<unsigned int> skel_corresponding_point;
		std::vector<float> dist_mid = skeleton_distance_to_midpoint(m, skeleton, skel_dist_vec);
		//descriptors[i].skel_dist_mid = dist_mid[0];
		skeleton_get_closest_skeleton_endpoints(m, skeleton, skel_dist_vec, skel_corresponding_point);
		descriptors[i].skeleton_index = skel_corresponding_point[0];
	}
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &m->vertices[descriptors[i].indices[0]]) > 0)
		{
			desc_pos.push_back(descriptors[i]);
		}
		else
		{
			desc_neg.push_back(descriptors[i]);
		}
	}

	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all(m, desc_pos, desc_neg, true, 20, 2);
	//std::vector<std::vector<float>> optimal_transforms_neg_pos = VarianceMin_compare_all(m, desc_neg, desc_left ,true , 20 , 3) ;
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;



	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			int dvoak_index_i;
			int dvoak_index_j;
			for (size_t k = 0; k < dvorak_pairs.size(); k++)
			{
				if (dvorak_pairs[k].p_index == desc_pos[i].indices[0])
				{
					dvoak_index_i = k;
				}
				if (dvorak_pairs[k].p_index == desc_neg[j].indices[0])
				{
					dvoak_index_j = k;
				}
			}
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[desc_pos[i].indices[0]] - m->normalized_heat_kernel_signature[desc_neg[j].indices[0]]);
			bool is_hks = hks_dif < 0.3;
			bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, 0.5, dvoak_index_i, dvoak_index_j);
			std::cout << " curv " << is_curv << std::endl;
			bool is_norm = dvorak_normal_angle_criterion(m, dvorak_pairs, dvoak_index_i, dvoak_index_j, 0.985);
			std::cout << " normal " << is_norm << std::endl;
			//bool is_skel_dist_far = std::abs(desc_pos[i].skel_dist_mid - desc_neg[j].skel_dist_mid) < 10;
			//if (is_curv && is_norm && is_skel_dist_far)
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = optimal_transforms_pos_neg[i][j];
				res.second = std::make_pair(i, desc_pos.size() + j);
				compare_results.push_back(res);
			}
		}
	}
	for (size_t i = 0; i < desc_neg.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_pos.size(); j++)
		{
			int dvoak_index_i;
			int dvoak_index_j;
			for (size_t k = 0; k < dvorak_pairs.size(); k++)
			{
				if (dvorak_pairs[k].p_index == desc_neg[i].indices[0])
				{
					dvoak_index_i = k;
				}
				if (dvorak_pairs[k].p_index == desc_pos[j].indices[0])
				{
					dvoak_index_j = k;
				}
			}
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[desc_neg[i].indices[0]] - m->normalized_heat_kernel_signature[desc_pos[j].indices[0]]);
			bool is_hks = hks_dif < 0.3;
			bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, 0.5, dvoak_index_i, dvoak_index_j);
			std::cout << " curv " << is_curv << std::endl;
			bool is_norm = dvorak_normal_angle_criterion(m, dvorak_pairs, dvoak_index_i, dvoak_index_j, 0.985);
			std::cout << " normal " << is_norm << std::endl;
			//bool is_skel_dist_far = std::abs(desc_pos[j].skel_dist_mid - desc_neg[i].skel_dist_mid) < 10;
			//if (is_curv && is_norm && is_skel_dist_far)
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = optimal_transforms_pos_neg[j][i];
				res.second = std::make_pair(desc_pos.size() + i,  j);
				compare_results.push_back(res);
			}
		}
	}
	std::vector<NLateralDescriptor> total_desc;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		total_desc.push_back(desc_pos[i]);
	}
	for (size_t i = 0; i < desc_neg.size(); i++)
	{
		total_desc.push_back(desc_neg[i]);

	}
	std::vector<bool> used(desc_pos.size() + desc_neg.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i]  && !used[j]) {
			resemblance_pairs.push_back({ total_desc[i].indices[0], total_desc[j].indices[0]});
			used[i] = true;  // Mark these objects as used
			used[j] = true;  // Mark these objects as used
		}
	}


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, right_skeleton_indices[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}


	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors[i].vertices_inside.size(); j++)
		{
			int index = descriptors[i].vertices_inside[j];
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}

	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors[i].paths.size(); j++)
		{
			for (size_t k = 0; k < descriptors[i].paths[j].size(); k++)
			{
				for (size_t t = 0; t < descriptors[i].paths[j][k].size(); t++)
				{
					int index = descriptors[i].paths[j][k][t];
					m->raylib_mesh.colors[index*4 ]= 255;
					m->raylib_mesh.colors[index*4 + 1]= 0;
					m->raylib_mesh.colors[index*4 + 2]= 0;
					m->raylib_mesh.colors[index*4 + 3]= 255;
				}
			}
		}
		
	}
	m->calculated_symmetry_pairs = resemblance_pairs;
	m->update_raylib_mesh();


	std::string path = "../../Results/";
	path = path + m->file_name + " NLateral_W_Gaussian_curvature_and_Optimal_Transform.txt ";
	Metric_write_to_file(m, path);
	
	return descriptors;
}
//check if two descriptros include each other in indices, 
//if so return false 
static bool is_descriptors_comparable(NLateralDescriptor& desc1 , NLateralDescriptor& desc2)
{
	bool is_comparable = true; 
	for (size_t i = 0; i < desc1.indices.size(); i++)
	{
		int index_i = desc1.indices[i];
		for (size_t j = 0; j < desc2.indices.size(); j++)
		{
			int index_j = desc2.indices[j];
			if (index_i == index_j)
			{
				is_comparable = false;
				break;
			}
		}
	}
	return is_comparable;
}

std::vector<NLateralDescriptor> NlateralMap_point_matching_with_FPS_and_endpoints(TrilateralMesh* m, Skeleton& skeleton,
	int dvorak_enpoint_no, float sweep_distance, float hks_dif_param, float curv_param, float norm_angle_param, float skel_dist_param, int skel_depth_param
	, float proximity, int N)
{
	std::vector<NLateralDescriptor> nlateral_descs; 
	std::vector<std::pair<int, float>> hks_pairs;
	std::vector<DvorakPairs> dvorak_pairs;
	hks_pairs = HKS_extraction_significant_points(m, dvorak_enpoint_no);
	hks_pairs = HKS_sweep_distance(m, hks_pairs, sweep_distance);
	std::sort(hks_pairs.begin(), hks_pairs.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
		return a.second < b.second; // Compare by second element
		});
	dvorak_pairs = HKS_to_dvorak_pairs(m, hks_pairs);
	

	std::vector<unsigned int> fps_points = furthest_point_sampling(m, 50, true);
	for (size_t i = 0; i < fps_points.size(); i++)
	{
		int index = fps_points[i];
		//find closest two
		std::vector<float> distances = Geodesic_dijkstra(*m, index);
		
		//get smallest two
		float smallest = INFINITY;
		float second_smallest = INFINITY;
		int smallest_index = -1; 
		int second_smallest_index = -1;
		for (size_t j = 0; j < hks_pairs.size(); j++)
		{
			if (hks_pairs[j].second < smallest)
			{
				smallest = hks_pairs[j].second;
				smallest_index = hks_pairs[j].first;
			}
		}
		for (size_t j = 0; j < hks_pairs.size(); j++)
		{
			if (hks_pairs[j].second < second_smallest && smallest_index != hks_pairs[j].first)
			{
				second_smallest = hks_pairs[j].second;
				second_smallest_index = hks_pairs[j].first;
			}
		}
		std::vector<unsigned int> points = { fps_points[i] , (unsigned int)smallest_index ,(unsigned int)second_smallest_index };

		NLateralDescriptor desc = NLateral_generate_descriptor(m, points);
		std::vector<unsigned int> skel_dist_vec = { desc.indices[0] };
		std::vector<float> dist_mid = skeleton_distance_to_midpoint(m, skeleton, skel_dist_vec);
		//desc.skel_dist_mid = dist_mid[0];
		std::vector<unsigned int> skel_corresponding_point;
		skeleton_get_closest_skeleton_endpoints(m, skeleton, skel_dist_vec, skel_corresponding_point);
		desc.skeleton_index = skel_corresponding_point[0];
		nlateral_descs.push_back(desc);
	}

	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all(m, nlateral_descs, true, 5, 1);
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	for (size_t i = 0; i < nlateral_descs.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < nlateral_descs.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			int hks_index_i;
			int hks_index_j;
			for (size_t k = 0; k < hks_pairs.size(); k++)
			{
				if (hks_pairs[k].first == nlateral_descs[i].indices[0])
				{
					hks_index_i = k;
				}
				if (hks_pairs[k].first == nlateral_descs[j].indices[0])
				{
					hks_index_j = k;
				}
			}
			std::cout << i << "  " << j << std::endl;
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[nlateral_descs[i].indices[0]] - m->normalized_heat_kernel_signature[nlateral_descs[j].indices[0]]);
			std::cout << " hks diff " << hks_dif << std::endl;
			bool is_hks = hks_dif < hks_dif_param;
			std::cout << " is_hks " << is_hks << std::endl;
			//bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, curv_param, hks_index_i, hks_index_j);
			//std::cout << " curv " << is_curv << std::endl;
			//bool is_norm = dvorak_normal_angle_criterion(m, dvorak_pairs, hks_index_i, hks_index_j, norm_angle_param);
			//std::cout << " normal " << is_norm << std::endl;
			float skel_dist = 0;// std::abs(nlateral_descs[i].skel_dist_mid - nlateral_descs[j].skel_dist_mid);
			bool is_skel_dist_far = skel_dist < skel_dist_param;
			std::cout << " skeld dist " << skel_dist << std::endl;
			std::cout << " is skeld dist " << is_skel_dist_far << std::endl;
			int skel_depth;// = std::abs(nlateral_descs[i].depth - nlateral_descs[j].depth);
			std::cout << " skeld depth " << skel_depth << std::endl;
			bool is_skel_depth = skel_depth < skel_depth_param;
			std::cout << " is skel depth " << is_skel_depth << std::endl;
			bool is_descs_comp = is_descriptors_comparable(nlateral_descs[i], nlateral_descs[j]);
			bool is_proximity = Geodesic_proximity(*m, nlateral_descs[i], nlateral_descs[j], proximity);

			std::cout << " comparable  " << is_descs_comp << std::endl;
			if (/*is_curv && is_norm && */ is_skel_dist_far && is_descs_comp && is_hks && is_skel_depth && !is_proximity)
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = optimal_transforms_pos_neg[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}
		}
	}

	std::vector<bool> used(nlateral_descs.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ nlateral_descs[i].indices[0], nlateral_descs[j].indices[0] });
			used[i] = true;  // Mark these objects as used
			used[j] = true;  // Mark these objects as used
		}
	}


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	//float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}


	for (size_t i = 0; i < nlateral_descs.size(); i++)
	{
		for (size_t j = 0; j < nlateral_descs[i].vertices_inside.size(); j++)
		{
			int index = nlateral_descs[i].vertices_inside[j];
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}

	for (size_t i = 0; i < nlateral_descs.size(); i++)
	{
		for (size_t j = 0; j < nlateral_descs[i].paths.size(); j++)
		{
			for (size_t k = 0; k < nlateral_descs[i].paths[j].size(); k++)
			{
				for (size_t t = 0; t < nlateral_descs[i].paths[j][k].size(); t++)
				{
					int index = nlateral_descs[i].paths[j][k][t];
					m->raylib_mesh.colors[index * 4] = 255;
					m->raylib_mesh.colors[index * 4 + 1] = 0;
					m->raylib_mesh.colors[index * 4 + 2] = 0;
					m->raylib_mesh.colors[index * 4 + 3] = 255;
				}
			}
		}
	}
	m->calculated_symmetry_pairs = resemblance_pairs;
	m->update_raylib_mesh();

	return nlateral_descs;
}
std::vector<NLateralDescriptor> NlateralMap_point_matching_with_skeleton_endpoints_and_OT_without_sym_plane(TrilateralMesh* m, Skeleton& skeleton,
int dvorak_enpoint_no,float sweep_distance, float hks_dif_param , float curv_param , float norm_angle_param, float skel_dist_param , int skel_depth_param
,float proximity,  int N)
{
	int size = m->vertices.size();
	int mesh_mid_point_index = -1;
	std::vector<NLateralDescriptor> desc_neg;
	std::vector<NLateralDescriptor> desc_pos;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<std::pair<int, float>> hks_pairs;
	std::vector<DvorakPairs> dvorak_pairs;
	hks_pairs = HKS_extraction_significant_points(m, dvorak_enpoint_no);
	//sweep 
	hks_pairs = HKS_sweep_distance(m,hks_pairs, sweep_distance);
	dvorak_pairs = HKS_to_dvorak_pairs(m, hks_pairs);
	Metric_set_gaussian(m, dvorak_enpoint_no, sweep_distance);
	Metric_set_N(N);


	//check skel dist param
	float max_skel_dist = -INFINITY;
	for (size_t i = 0; i < skeleton.skeletonFormat.size(); i++)
	{
		std::vector<int> vertex_list;
		std::vector<float> distances;
		skeleton_calculate_dijkstra(skeleton, i, vertex_list, distances);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (distances[j] > max_skel_dist)
			{
				max_skel_dist = distances[j];
			}
		}
	}
	//skel_dist_param =  max_skel_dist / 32;

	std::vector<unsigned int > skel_mid_point_indices;
	for (size_t i = 0; i < hks_pairs.size(); i++)
	{
		skel_mid_point_indices.push_back(hks_pairs[i].first);
	}
	skeleton_distance_to_midpoint(m, skeleton, skel_mid_point_indices);

	std::vector<unsigned int > mesh_vertices;
	SkeletonTree skelTree = skeleton_generate_skeleton_tree(m, skeleton);
	for (size_t i = 0; i < hks_pairs.size(); i++)
	{
		mesh_vertices.push_back(hks_pairs[i].first);
	}
	std::vector<NLateralDescriptor> descriptors;
	descriptors = NLateral_generate_closest_points(m,mesh_vertices, N, 10);
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		std::vector<unsigned int> skel_dist_vec = { descriptors[i].indices[0] };
		std::vector<unsigned int> skel_corresponding_point;
		std::vector<float> dist_mid = skeleton_distance_to_midpoint(m, skeleton, skel_dist_vec);
		//descriptors[i].skel_dist_mid = dist_mid[0];
		skeleton_get_closest_skeleton_endpoints(m, skeleton, skel_dist_vec, skel_corresponding_point);
		descriptors[i].skeleton_index = skel_corresponding_point[0];
	}

	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all(m, descriptors, true, 5, 1);
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	for (size_t i = 0; i < descriptors.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < descriptors.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			int hks_index_i;
			int hks_index_j;
			for (size_t k = 0; k < hks_pairs.size(); k++)
			{
				if (hks_pairs[k].first == descriptors[i].indices[0])
				{
					hks_index_i = k;
				}
				if (hks_pairs[k].first == descriptors[j].indices[0])
				{
					hks_index_j = k;
				}
			}
			std::cout << i << "  " << j << std::endl;
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descriptors[i].indices[0]] - m->normalized_heat_kernel_signature[descriptors[j].indices[0]]);
			std::cout << " hks diff " << hks_dif << std::endl; 
			bool is_hks = hks_dif < hks_dif_param;
			std::cout << " is_hks " << is_hks << std::endl;
			bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, curv_param, hks_index_i, hks_index_j);
			std::cout << " curv " << is_curv << std::endl;
			bool is_norm = dvorak_normal_angle_criterion(m, dvorak_pairs, hks_index_i, hks_index_j, norm_angle_param);
			std::cout << " normal " << is_norm << std::endl;
			float skel_dist = 0;// std::abs(descriptors[i].skel_dist_mid - descriptors[j].skel_dist_mid);
			bool is_skel_dist_far = skel_dist < skel_dist_param;
			std::cout << " skeld dist " <<  skel_dist << std::endl;
			std::cout << " is skeld dist " << is_skel_dist_far  << std::endl;
			int skel_depth;// = std::abs(descriptors[i].depth - descriptors[j].depth);
			std::cout << " skeld depth " << skel_depth << std::endl;
			bool is_skel_depth = skel_depth < skel_depth_param;
			std::cout << " is skel depth " << is_skel_depth  << std::endl;
			bool is_descs_comp = is_descriptors_comparable( descriptors[i] , descriptors[j]);
			bool is_proximity = Geodesic_proximity(*m,descriptors[i], descriptors[j] , proximity);

			std::cout << " comparable  " << is_descs_comp << std::endl;
			if (is_curv && is_norm && is_skel_dist_far && is_descs_comp && is_hks && is_skel_depth && !is_proximity )
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = optimal_transforms_pos_neg[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}
		}
	}

	std::vector<bool> used(descriptors.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ descriptors[i].indices[0], descriptors[j].indices[0] });
			used[i] = true;  // Mark these objects as used
			used[j] = true;  // Mark these objects as used
		}
	}


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	//float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}


	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors[i].vertices_inside.size(); j++)
		{
			int index = descriptors[i].vertices_inside[j];
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}

	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors[i].paths.size(); j++)
		{
			for (size_t k = 0; k < descriptors[i].paths[j].size(); k++)
			{
				for (size_t t = 0; t < descriptors[i].paths[j][k].size(); t++)
				{
					int index = descriptors[i].paths[j][k][t];
					m->raylib_mesh.colors[index * 4] = 255;
					m->raylib_mesh.colors[index * 4 + 1] = 0;
					m->raylib_mesh.colors[index * 4 + 2] = 0;
					m->raylib_mesh.colors[index * 4 + 3] = 255;
				}
			}
		}

	}
	m->calculated_symmetry_pairs = resemblance_pairs;
	m->update_raylib_mesh();


	//std::string path = "../../Results/";
	//path = path + m->file_name + " NLateral_W_OT.txt ";
	//Metric_write_to_file(m, path);

	return descriptors;
}


std::vector<NLateralDescriptor> NlateralMap_point_matching_with_skeleton_endpoints_and_OT_without_sym_plane_FPS(TrilateralMesh* m, Skeleton& skeleton,
	int dvorak_enpoint_no, float sweep_distance, float hks_dif_param, float curv_param, float norm_angle_param, float skel_dist_param, float n_ring_param,
	float area_dif_param, int N)
{

	int size = m->vertices.size();
	int mesh_mid_point_index = -1;
	std::vector<NLateralDescriptor> desc_neg;
	std::vector<NLateralDescriptor> desc_pos;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<std::pair<int, float>> hks_pairs;
	std::vector<DvorakPairs> dvorak_pairs;

	std::vector<unsigned int> fps_indices = furthest_point_sampling(m, dvorak_enpoint_no, true);
	//hks_pairs = HKS_extraction_significant_points(m, dvorak_enpoint_no);
	//sweep 
	//hks_pairs = HKS_sweep_distance(m, hks_pairs, sweep_distance);
	for (size_t i = 0; i < fps_indices.size(); i++)
	{
		DvorakPairs p;
		p.p_index = fps_indices[i];
		p.gaussian_curv = gaussian_curvature(m ,p.p_index);
		dvorak_pairs.push_back(p);
	}
	//dvorak_pairs = HKS_to_dvorak_pairs(m, hks_pairs);
	Metric_set_gaussian(m, dvorak_enpoint_no, sweep_distance);
	Metric_set_N(N);


	//check skel dist param
	float max_skel_dist = -INFINITY;
	for (size_t i = 0; i < skeleton.skeletonFormat.size(); i++)
	{
		std::vector<int> vertex_list;
		std::vector<float> distances;
		skeleton_calculate_dijkstra(skeleton, i, vertex_list, distances);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (distances[j] > max_skel_dist)
			{
				max_skel_dist = distances[j];
			}
		}
	}
	//skel_dist_param =  max_skel_dist / 32;

	std::vector<unsigned int > skel_mid_point_indices;
	for (size_t i = 0; i < fps_indices.size(); i++)
	{
		skel_mid_point_indices.push_back(fps_indices[i]);
	}
	skeleton_distance_to_midpoint(m, skeleton, skel_mid_point_indices);

	SkeletonTree skelTree = skeleton_generate_skeleton_tree(m, skeleton);

	std::vector<NLateralDescriptor> descriptors;
	descriptors = NLateral_generate_closest_points(m, fps_indices, N,10); 
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		std::vector<unsigned int> skel_dist_vec = { descriptors[i].indices[0] };
		std::vector<unsigned int> skel_corresponding_point;
		std::vector<float> dist_mid = skeleton_distance_to_midpoint(m, skeleton, skel_dist_vec);
		//descriptors[i].skel_dist_mid = dist_mid[0];
		skeleton_get_closest_skeleton_endpoints(m, skeleton, skel_dist_vec, skel_corresponding_point);
		descriptors[i].skeleton_index = skel_corresponding_point[0];
	}
// maximum n ring 
	float maximum_n_ring = -INFINITY;
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors.size(); j++)
		{
			float dif = std::fabs(descriptors[i].n_ring_area - descriptors[j].n_ring_area);
			if (maximum_n_ring < dif)
			{
				maximum_n_ring = dif ;
			}
		}
	}
//maximum skel distance
	float maximum_skel_dist = -INFINITY;
	for (size_t i = 0; i < skeleton.skeletonFormat.size(); i++)
	{
		std::vector<int> vertex_list; 
		std::vector<float> vertex_dist; 
		skeleton_calculate_dijkstra(skeleton, i, vertex_list, vertex_dist);
		for (size_t j = 0; j < vertex_dist.size(); j++)
		{
			if (maximum_skel_dist < vertex_dist[j])
			{
				maximum_skel_dist = vertex_dist[j];
			}
		}
	}
//maximum area dif
	float maximum_area_dif = -INFINITY;
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		float area_i = descriptors[i].area;
		for (size_t j = i; j < descriptors.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			float area_j = descriptors[j].area;
			float dif = std::abs(area_i - area_j);
			if (dif > maximum_area_dif)
			{
				maximum_area_dif = dif; 
			}

		}
	}
	std::vector<std::vector<float>> optimal_transforms_pos_neg = VarianceMin_compare_all(m, descriptors, true, 20, 1);
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	std::ofstream file;
	file.open("../../Results/desc.txt");
	for (size_t i = 0; i < descriptors.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < descriptors.size(); j++)
		{
			if (i >= j)
			{
				continue;
			}
			int fps_index_i;
			int fps_index_j;
			for (size_t k = 0; k < fps_indices.size(); k++)
			{
				if (fps_indices[k] == descriptors[i].indices[0])
				{
					fps_index_i = k;
				}
				if (fps_indices[k] == descriptors[j].indices[0])
				{
					fps_index_j = k;
				}
			}
			file << i << "  " << j << std::endl;
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descriptors[i].indices[0]] - m->normalized_heat_kernel_signature[descriptors[j].indices[0]]);
			file << " hks diff " << hks_dif << std::endl;
			bool is_hks = hks_dif < hks_dif_param;
			file  << " is_hks " << is_hks << std::endl;
			bool is_curv = dvorak_curvature_similarity_criterion(dvorak_pairs, curv_param, fps_index_i, fps_index_j);
			file << " curv " << is_curv << std::endl;
			bool is_norm = dvorak_normal_angle_criterion(m, dvorak_pairs, fps_index_i, fps_index_j, norm_angle_param);
			file << " normal " << is_norm << std::endl;
			float skel_dist = 0;// std::abs(descriptors[i].skel_dist_mid - descriptors[j].skel_dist_mid);
			bool is_skel_dist_far = ( skel_dist / maximum_skel_dist ) < skel_dist_param;
			file << " skeld dist " << skel_dist << std::endl;
			file << " is skeld dist " << is_skel_dist_far << std::endl;
			bool is_descs_comp = is_descriptors_comparable(descriptors[i], descriptors[j]);
			file << " comparable  " << is_descs_comp << std::endl;
			float n_ring = std::abs(descriptors[i].n_ring_area - descriptors[j].n_ring_area);
			file << " n ring " << n_ring << std::endl;
			bool is_n_ring_close = (n_ring / maximum_n_ring ) < n_ring_param; 
			float area_dif = std::abs(descriptors[i].area- descriptors[j].area);
			bool is_area_dif = area_dif / maximum_area_dif < area_dif_param;
			
			if (/*is_curv && is_norm && */ is_skel_dist_far /* && is_descs_comp */ && is_hks && is_n_ring_close && is_area_dif)
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = optimal_transforms_pos_neg[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}
		}
	}

	std::vector<bool> used(descriptors.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] /* && !used[j]*/) {
			resemblance_pairs.push_back({ descriptors[i].indices[0], descriptors[j].indices[0] });
			used[i] = true;  // Mark these objects as used
			//used[j] = true;  // Mark these objects as used
		}
	}


	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	//float total_error = Metric_get_geodesic_cost_with_list(m, left_correspondences, right_correspondences);

	// color left red
	std::vector<unsigned int> is_selected(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 255;
		m->colors[resemblance_pairs[i].second].b = 0;
	}


	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors[i].vertices_inside.size(); j++)
		{
			int index = descriptors[i].vertices_inside[j];
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}

	for (size_t i = 0; i < descriptors.size(); i++)
	{
		for (size_t j = 0; j < descriptors[i].paths.size(); j++)
		{
			for (size_t k = 0; k < descriptors[i].paths[j].size(); k++)
			{
				for (size_t t = 0; t < descriptors[i].paths[j][k].size(); t++)
				{
					int index = descriptors[i].paths[j][k][t];
					m->raylib_mesh.colors[index * 4] = 255;
					m->raylib_mesh.colors[index * 4 + 1] = 0;
					m->raylib_mesh.colors[index * 4 + 2] = 0;
					m->raylib_mesh.colors[index * 4 + 3] = 255;
				}
			}
		}

	}
	m->calculated_symmetry_pairs = resemblance_pairs;
	m->update_raylib_mesh();


	return descriptors;
}