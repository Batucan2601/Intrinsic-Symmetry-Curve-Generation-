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
#include <src/Application/Include/CoreTypeDefs.h>

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

//static int* trialteral_ROI(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no, bool& is_visited_interior);

std::vector<float> compute_geodesic_distances_min_heap_distances(Mesh& m, int point_index);

void trilateral_map_drawing_using_three_points(MeshFactory& mesh_fac, int& selected_index, int p1, int p2, int p3);


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
std::vector<TrilateralDescriptor> get_trilateral_points_using_closest_pairs(MeshFactory& mesh_fac, const int& selected_index, std::vector<unsigned int>& indices);
std::vector<TrilateralDescriptor> get_trilateral_points_using_closest_pairs(Mesh*m, std::vector<unsigned int>& indices);
std::vector<unsigned int> AverageGeodesicFunction(MeshFactory& mesh_fac, int& selected_index, int& number_of_points);
std::vector<unsigned int> minimumGeodesicFunction(MeshFactory& mesh_fac, int& selected_index, int& number_of_points, std::vector<unsigned int>& average_geodesic_function);


//static int* trialteral_ROI(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no, bool& is_visited_interior);
std::vector<float>  histogramROi(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, int division_no, int* is_visited, bool is_visited_interior);

// the ultimate trilateral descriptor generator 
TrilateralDescriptor  generate_trilateral_descriptor(MeshFactory& mesh_fac, int& selected_index, int point_index1, int point_index2, int point_index3, bool is_simplified);
TrilateralDescriptor  generate_trilateral_descriptor(Mesh* m, int point_index1, int point_index2, int point_index3, bool is_simplified);


std::vector<std::pair<unsigned int, unsigned int>>  point_match_trilateral_weights(Mesh*m, std::vector<TrilateralDescriptor>& trilateralDescVec, const float& curvWeight,
	const float& geodesicWeight, const float& areaWeight);
std::vector<std::pair<unsigned int, unsigned int>>  point_match_trilateral_weights(Mesh* m, std::vector<TrilateralDescriptor>& trilateralDescVecLeft, std::vector<TrilateralDescriptor>& trilateralDescVecRight,const float& curvWeight,
	const float& geodesicWeight, const float& areaWeight);
void display_accuracy(Mesh* m, std::vector<std::pair<unsigned int, unsigned int>>& calculated_symmetry_pairs);

void point_matching_with_dominant_symmetry_plane(MeshFactory& mesh_fac, int& selected_index, Plane* plane, int sampling_no);

// sampling
void simple_sample(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int sample_size, int division_no);

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
void match_points_from2_mesh(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int division_no);

//from the paper Robust3DShapeCorrespondenceintheSpectralDomain 4.1 
std::vector<glm::vec3> generate_spectral_embedding(MeshFactory& meshFac, int mesh_index, std::vector<unsigned int> landmark_vertices);



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

 void reset_points(MeshFactory& mesh_fac, int meshIndex);