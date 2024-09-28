#pragma once

#include "../Include/MeshFactory.h"
#include <GL/glew.h>
#include <utility>
#include <algorithm>
#include <queue>
#include <stack>
#include <algorithm> 
#include <math.h>
#include <eigen/Eigen/Dense>
#include "Sampling.h"
#include "CoreTypeDefs.h"
#include "FuzzyGeodesic.h"
#include <src/Application/Include/CoreTypeDefs.h>
#include "../Include/Skeleton.h"

// 1 - use area error
// 2 - use geodesic error
// 3 - use euclidean error
// 4 - use curvature error
struct TrilateralError
{
	float areaError;
	float geodesicError;
	float euclideanError;
	float curvatureError;
};

//static int* trilateral_ROI(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no, bool& is_visited_interior);

std::vector<float> compute_geodesic_distances_min_heap_distances(Mesh& m, int point_index);

void trilateral_map_drawing_using_three_points(MeshFactory& mesh_fac, int& selected_index, int p1, int p2, int p3);


/*static float* trilateral_ROI(MeshFactory& mesh_fac  , int& selected_index, int point_index1, int point_index2, int point_index3, int division_no) // tau is the closeness division_no is the no of how much you want to separate
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
std::vector<TrilateralDescriptor> get_trilateral_points_using_closest_pairs(MeshFactory& mesh_fac, const int& selected_index, std::vector<unsigned int>& indices);
std::vector<TrilateralDescriptor> get_trilateral_points_using_closest_pairs(Mesh*m, std::vector<unsigned int>& indices);
std::vector<TrilateralDescriptor> get_trilateral_points_using_furthest_pairs(Mesh*m, std::vector<unsigned int>& indices);
TrilateralDescriptor trilateral_get_trilateral_using_closest_pairs_with_skeleton_indices(Mesh*m, unsigned int point_index, 
std::vector<unsigned int>& skeleton_indices);


std::vector<unsigned int> AverageGeodesicFunction(MeshFactory& mesh_fac, int& selected_index, int& number_of_points);
std::vector<unsigned int> minimumGeodesicFunction(MeshFactory& mesh_fac, int& selected_index, int& number_of_points, std::vector<unsigned int>& average_geodesic_function);


//static int* trilateral_ROI(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no, bool& is_visited_interior);
std::vector<float> histogramROi(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited);
std::vector<float> histogramROi_w_HKS(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited);
std::vector<float> histogram_roi_superior(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited);

// the ultimate trilateral descriptor generator 
TrilateralDescriptor  generate_trilateral_descriptor(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, bool is_simplified);
TrilateralDescriptor  generate_trilateral_descriptor(Mesh* m, int point_index1, int point_index2, int point_index3, bool is_simplified);


std::vector<std::pair<unsigned int, unsigned int>>  point_match_trilateral_weights(Mesh*m, std::vector<TrilateralDescriptor>& trilateralDescVec, const float& curvWeight,
	const float& geodesicWeight, const float& areaWeight);
std::vector<std::pair<unsigned int, unsigned int>>  point_match_trilateral_weights(Mesh* m, std::vector<TrilateralDescriptor>& trilateralDescVecLeft, std::vector<TrilateralDescriptor>& trilateralDescVecRight,const float& curvWeight,
	const float& geodesicWeight, const float& areaWeight);
void display_accuracy(Mesh* m, std::vector<std::pair<unsigned int, unsigned int>>& calculated_symmetry_pairs);

void point_matching_with_dominant_symmetry_plane(MeshFactory& mesh_fac, int& selected_index, Plane* plane, int sampling_no);

//region of interest
std::vector<int> trilateral_ROI(Mesh* m, int point_index1, int point_index2, int point_index3, int division_no, bool is_color);
void trilateral_ROI_area(Mesh* m, const std::vector<int>& trilateral_vertices,  float& total_area);


// sampling
void simple_sample(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int sample_size, int division_no);

void match_points_from2_mesh(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int division_no);

//from the paper Robust3DShapeCorrespondenceintheSpectralDomain 4.1 
std::vector<glm::vec3> generate_spectral_embedding(MeshFactory& meshFac, int mesh_index, std::vector<unsigned int> landmark_vertices);


float get_N_ring_area(Mesh* m, float point_index , int N );


 void reset_points(MeshFactory& mesh_fac, int meshIndex);

 // spin image creation 
 std::vector<float> trilateral_generate_spin_image(MeshFactory& mesh_fac, int selected_index, std::vector<int>& vertices_in_tri_area, int division_no);
 float chi_squre_distance(std::vector<float>& vec1, std::vector<float>& vec2);

 void trilateral_FPS_histogram_matching(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no, bool recordTxt);
 void trilateral_FPS_histogram_matching_w_spin_image(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no);
 void trilateral_FPS_histogram_matching_w_spin_image_MDS(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no);
 void trilateral_FPS_histogram_matching_w_principal_comp(MeshFactory& mesh_fac, const int& selected_index, int sample_no, int division_no);

 static std::vector<float> computePrincipalCurvatures(MeshFactory& meshFac, int selectedIndex, std::vector<int>& is_visited,
	 int division_no);

 void trilateral_self_matching_with_dominant_sym(MeshFactory& mesh_fac, int selected_index, const int number_of_n_lateral_points, const std::string& method_name
, int sym_iter_no,bool is_LMDS );

 //using fuzzy geodesics 

 void trilateral_fuzzyGeodesic(MeshFactory& meshFac, int selectedIndex, int p1, int p2, int p3, float fuzziness_sigma);
 void trilateral_FPS_matching_w_fuzzy_geodesic(MeshFactory& mesh_fac, const int& selected_index, int sample_no , float fuzziness_sigma, bool recordTxt );

 //skeleton and fuzzy
 void trilateral_w_skeleton_endpoints(MeshFactory& mesh_fac, const int& selected_index,
	 float fuzziness_sigma, Skeleton& skeleton, bool recordTxt);


 void trilateral_point_matching_with_skeleton_endpoints(MeshFactory& mesh_fac, const int& selected_index,Skeleton& skeleton );


 void trilateral_point_matching_with_skeleton_endpoints_w_HKS(MeshFactory& mesh_fac, const int& selected_index, Skeleton& skeleton,
std::vector<TrilateralDescriptor>& desc_left,std::vector<TrilateralDescriptor>& desc_right , Plane& plane);

 void trilateral_point_matching_with_skeleton_endpoints_anchors(MeshFactory& mesh_fac, const int& selected_index, Skeleton& skeleton,
	 std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& desc_neg, Plane& plane );