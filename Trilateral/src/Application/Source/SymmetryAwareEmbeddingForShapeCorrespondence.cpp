#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#pragma comment(linker, "/STACK:2000000")
#pragma comment(linker, "/HEAP:1000000000")

Plane generate_isomap_embedding(Mesh* mesh , bool simplify_mesh , float simplification_percentage )
{

	int M = 3 ; //dimensions 

	//simplification part
	Mesh simplified_mesh = *mesh; 
	if (simplify_mesh)
	{
		srand(time(0));

		std::vector<std::pair<unsigned int , unsigned int>> random_indices_list;
		int size = simplified_mesh.vertices.size() * simplification_percentage / 100.0f ;
		//new symmetry pairs
		while( random_indices_list.size() < size / 2 )
		{
			unsigned int rand_no = rand() % mesh->vertices.size();
			bool is_exists = false;
			for (size_t i = 0; i < random_indices_list.size() ; i++)
			{
				if (rand_no == random_indices_list[i].first || rand_no == random_indices_list[i].second)
				{
					is_exists = true; 
					break; 
				}
			}
			if (!is_exists)
			{
				std::pair<unsigned int, unsigned int> pair;
				pair.first = rand_no; 
				for (size_t i = 0; i < mesh->symmetry_pairs.size(); i++)
				{
					if (pair.first == mesh->symmetry_pairs[i].first)
					{
						pair.second = mesh->symmetry_pairs[i].second;
						break;
					}
					else if (pair.first == mesh->symmetry_pairs[i].second)
					{
						pair.second = mesh->symmetry_pairs[i].first;
						break;
					}
				}
				random_indices_list.push_back(pair);
			}

		}
		simplified_mesh.symmetry_pairs = random_indices_list;
	}

	int N = simplified_mesh.symmetry_pairs.size() * 2;

	Eigen::MatrixXf delta(N, N);
	// generate geodesic distances matrix ( delta)
	if (!simplify_mesh)
	{
		for (size_t i = 0; i < N; i++)
		{
			std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, i);
			for (size_t j = 0; j < N; j++)
			{
				delta(i, j) = distances[j];
			}
		}
	}
	else
	{
		for (size_t i = 0; i < N/2; i++)
		{
			std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, simplified_mesh.symmetry_pairs[i].first);
			int delta_index = 0; 
			for (size_t j = 0; j < mesh->vertices.size(); j++)
			{
				bool is_j_index_allowed = false;
				for (size_t k = 0; k < N/2; k++)
				{
					if (simplified_mesh.symmetry_pairs[k].first == j || simplified_mesh.symmetry_pairs[k].second == j)
					{
						is_j_index_allowed = true; 
						break;
					}
				}
				if (is_j_index_allowed)
				{
					delta(i, delta_index++) = distances[j];
				}
			}
			//repeat for secone 
			distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, simplified_mesh.symmetry_pairs[i].second);
			delta_index = 0;
			for (size_t j = 0; j < mesh->vertices.size(); j++)
			{
				bool is_j_index_allowed = false;
				for (size_t k = 0; k < N / 2; k++)
				{
					if (simplified_mesh.symmetry_pairs[k].first == j || simplified_mesh.symmetry_pairs[k].second == j)
					{
						is_j_index_allowed = true;
						break;
					}
				}
				if (is_j_index_allowed)
				{
					delta(i, delta_index++) = distances[j];
				}
			}
		}
	}

	//generate J
	Eigen::MatrixXf J(N, N);
	J.fill(-1.0f / (float)N);
	J.diagonal(1 - (1.0f / (float)N) );

	//generate K
	Eigen::MatrixXf K(N, N);
	// K = -1/2 J delta^2 J
	K = -1 / 2 * J * (delta * delta) * J;

	//eigendecompose isomap

	//solves the eignevector
	Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
	eigensolver.compute(K);
	Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
	Eigen::MatrixXf eigen_vectors = eigensolver.eigenvectors().real();

	//get best 3 eigenvalues
	std::vector<unsigned int> best_top_eigen_value_indices; 
	int best_top_eigen_values = 4; 
	for (size_t i = 0; i < best_top_eigen_values; i++)
	{
		unsigned int best_index = -1;
		float best_value = -INFINITE;
		for (size_t j = 0; j < eigen_values.rows(); j++)
		{
			bool is_best_selected_already = false; 
			for (size_t k = 0; k < best_top_eigen_value_indices.size(); k++)
			{
				if (best_top_eigen_value_indices[k] == j)
				{
					is_best_selected_already = true; 
					break;
				}
			}
			if (is_best_selected_already)
			{
				continue;
			}
			if (best_value < eigen_values(j) )
			{
				best_value = eigen_values(j);
				best_index = j;
			}
		}
		best_top_eigen_value_indices.push_back(best_index);
	}
	// diagonalize and take sqrt
	Eigen::MatrixXf eigen_values_diagonal(M , M );
	eigen_values_diagonal.fill(0);
	for (size_t i = 0; i < M; i++)
	{
		eigen_values_diagonal(i, i) = sqrt(eigen_values(i));
	}

	glm::vec3 vectors[4];
	vectors[0] = glm::vec3(eigen_vectors.col(best_top_eigen_value_indices[0])(0), eigen_vectors.col(best_top_eigen_value_indices[0])(1), eigen_vectors.col(best_top_eigen_value_indices[0])(2));
	vectors[1] = glm::vec3(eigen_vectors.col(best_top_eigen_value_indices[1])(0) , eigen_vectors.col(best_top_eigen_value_indices[1])(1) , eigen_vectors.col(best_top_eigen_value_indices[1])(2));
	vectors[2] = glm::vec3(eigen_vectors.col(best_top_eigen_value_indices[2])(0) , eigen_vectors.col(best_top_eigen_value_indices[2])(1) , eigen_vectors.col(best_top_eigen_value_indices[2])(2));
	vectors[3] = glm::vec3(eigen_vectors.col(best_top_eigen_value_indices[3])(0) , eigen_vectors.col(best_top_eigen_value_indices[3])(1) , eigen_vectors.col(best_top_eigen_value_indices[3])(2));
	// take eigen value pairs
	// stated 1-2 1-3 1-4  2-1 3-1 4-1
	//generate 6 plane 

	Plane planes[6];
	planes[0] = generate_plane_from_two_vectors(vectors[0]  , vectors[1]);
	planes[1] = generate_plane_from_two_vectors(vectors[0], vectors[2]);
	planes[2] = generate_plane_from_two_vectors(vectors[0], vectors[3]);

	planes[3] = generate_plane_from_two_vectors(vectors[1], vectors[0]);
	planes[4] = generate_plane_from_two_vectors(vectors[2], vectors[0]);
	planes[5] = generate_plane_from_two_vectors(vectors[3], vectors[0]);

	float scores[6];
	
	scores[0] = generate_symmetry_score(*mesh, &planes[0]);
	scores[1]= generate_symmetry_score(*mesh,  &planes[1]);
	scores[2] = generate_symmetry_score(*mesh, &planes[2]);

	scores[3]= generate_symmetry_score(*mesh, &planes[3]);
	scores[4]= generate_symmetry_score(*mesh, &planes[4]);
	scores[5]= generate_symmetry_score(*mesh, &planes[5]);


	float best_score = -INFINITY; 
	int best_index = -1; 
	for (size_t i = 0; i < 6 ; i++)
	{
		if (scores[i] > best_score)
		{
			best_score = scores[i];
			best_index = i;
		}
	}

	return planes[best_index];
}

float generate_symmetry_score(Mesh mesh, Plane* p1 )
{
	// 1-  project the mesh triangles into the plane p1 
	std::vector<unsigned int> triangles_positive;
	std::vector<unsigned int> triangles_negative;
	// project the triangles only that are on the same side
	for (size_t i = 0; i < mesh.triangles.size(); i+=3 )
	{
		float status_p1 = get_point_status_from_plane(p1 , &mesh.vertices[mesh.triangles[i]]);
		float status_p2 = get_point_status_from_plane(p1 , &mesh.vertices[mesh.triangles[i+1]]);
		float status_p3 = get_point_status_from_plane(p1 , &mesh.vertices[mesh.triangles[i+2]]);
		if ( (status_p1 >= 0 && status_p2 >= 0 && status_p3 >= 0 )
		|| (status_p1 < 0 && status_p2 < 0 && status_p3 < 0) ) 
		{
			if (status_p1 >= 0)
			{
				triangles_positive.push_back(mesh.triangles[i]);
				triangles_positive.push_back(mesh.triangles[i+1]);
				triangles_positive.push_back(mesh.triangles[i+2]);

			}
			else
			{
				triangles_negative.push_back(mesh.triangles[i]);
				triangles_negative.push_back(mesh.triangles[i+1]);
				triangles_negative.push_back(mesh.triangles[i+2]);
			}
		}
	}
	
	//project
	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		mesh.vertices[i] = project_point_to_plane(p1, &mesh.vertices[i]);
	}

	float total_area_positive = 0;
	float total_area_negative = 0;
	float biggest_total_area = 0;
	//calculate total area for each of the sides
	for (size_t i = 0; i < triangles_positive.size(); i+=3)
	{
		//????????????????????
		total_area_positive += compute_triangle_area(mesh.vertices[triangles_positive[i]], mesh.vertices[triangles_positive[i+1]], mesh.vertices[triangles_positive[i+2]]);
	}
	for (size_t i = 0; i < triangles_negative.size(); i+=3)
	{
		total_area_negative += compute_triangle_area(mesh.vertices[triangles_negative[i]], mesh.vertices[triangles_negative[i + 1]], mesh.vertices[triangles_negative[i + 2]]);
	}
	
	if (total_area_positive > total_area_negative)
	{
		biggest_total_area = total_area_positive;
		return total_area_negative / biggest_total_area;
	}
	else
	{
		biggest_total_area = total_area_negative;
		return total_area_positive / biggest_total_area;

	}


	//for now much more simplified
	//for now just get  the number of intersections  rather than area 
	//TODO: will be improved later with the areas
	

	// find the total area intersection
	/*float intersection_area = 0;
	for (size_t i = 0; i < triangles_positive.size(); i+=3)
	{
		for (size_t j = 0; j < triangles_negative.size(); j+3)
		{
			bool is_intersect = is_triangles_intersect(mesh.vertices[triangles_positive[i]], mesh.vertices[triangles_positive[i + 1]], mesh.vertices[triangles_positive[i + 2]],
				mesh.vertices[triangles_negative[j]], mesh.vertices[triangles_negative[j + 1]], mesh.vertices[triangles_negative[j + 2]]);
			if (is_intersect)
			{
				intersection_area += 1; 
			}
		}
	}*/

	
	//return intersection_area / biggest_total_area;
}