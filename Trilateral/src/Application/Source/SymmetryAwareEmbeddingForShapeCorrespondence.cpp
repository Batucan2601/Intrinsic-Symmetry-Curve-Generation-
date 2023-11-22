#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/TrilateralMap.h"
#include "../Include/MetricCalculations.h"
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
		//int size = simplified_mesh.vertices.size() * simplification_percentage / 100.0f ;
		int size = 50;
		//new symmetry pairs
		/*while( random_indices_list.size() < size / 2 )
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

		}*/
		
		std::vector<unsigned int> fps_indices = furthest_point_sampling(&simplified_mesh, size);
		for (size_t i = 0; i < size; i++)
		{
			std::pair<unsigned int, unsigned int> pair;
			pair.first = fps_indices[i];
			for (size_t j = 0; j < size; j++)
			{
				if (pair.first == mesh->symmetry_pairs[j].first)
				{
					pair.second = mesh->symmetry_pairs[j].second;
					break;
				}
				else if (pair.first == mesh->symmetry_pairs[j].second)
				{
					pair.second = mesh->symmetry_pairs[j].first;
					break;
				}
			}
			random_indices_list.push_back(pair);

		}
		simplified_mesh.symmetry_pairs = random_indices_list;
	}

	int N = simplified_mesh.symmetry_pairs.size() ;

	Eigen::MatrixXf delta(N, N);
	Eigen::MatrixXf J(N, N);
	Eigen::MatrixXf K(N, N);

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
	J.fill(-1.0f / (float)N);
	for (size_t i = 0; i < N; i++)
	{
		J(i, i) = 1 - (1.0f / (float)N);
	}

	//generate K
	// K = -1/2 J delta^2 J
	Eigen::MatrixXf delta_sqr = delta.array().square();
	K = -0.5 * J * delta_sqr * J;

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
	planes[0] = generate_plane_from_two_vectors(vectors[0] , vectors[1]);
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


// Extract the N-largest eigen values and eigen vectors
inline void ExtractNLargestEigens(unsigned n, Eigen::VectorXd& S, Eigen::MatrixXd& V)
{
	// Note: m is the original dimension
	const unsigned m = S.rows();

	// Copy the original matrix
	const Eigen::MatrixXd original_V = V;

	// Sort by eigenvalue
	constexpr double                         epsilon = 1e-16;
	std::vector<std::pair<double, unsigned>> index_value_pairs(m);
	for (unsigned i = 0; i < m; ++i)
	{
		index_value_pairs[i] = std::make_pair(std::max(S(i), epsilon), i);
	}
	std::partial_sort(index_value_pairs.begin(),
		index_value_pairs.begin() + n,
		index_value_pairs.end(),
		std::greater<std::pair<double, unsigned>>());

	// Resize matrices
	S.resize(n);
	V.resize(m, n);

	// Set values
	for (unsigned i = 0; i < n; ++i)
	{
		S(i) = index_value_pairs[i].first;
		V.col(i) = original_V.col(index_value_pairs[i].second);
	}
	/*std::vector<unsigned int> indices;
	for (size_t i = 0; i < n; i++)
	{
		float best_value = -INFINITY;
		int best_index = -1;
		for (size_t j = 0; j < S.rows(); j++)
		{
			if (best_value < S(j))
			{
				bool is_already_selected = false; 
				for (size_t k = 0; k < indices.size(); k++)
				{
					if (indices[k] == j)
					{
						is_already_selected = true; 
						break; 
					}
				}
				if (!is_already_selected)
				{
					best_index = j;
					best_value = S(j);
				}
			}
			
		}
		indices.push_back(best_index);
	}
	
	for (size_t i = 0; i < n; i++)
	{
		float temp = S(indices[i]);
		S(indices[i]) = S(i);
		S(i) = temp;

		Eigen::VectorXd temp_mat = V.col(indices[i]);
		V.col(indices[i]) = V.col(i);
		V.col(i) = temp_mat;
	}*/

}

Eigen::MatrixXd ComputeClassicalMds(const Eigen::MatrixXd& D, const unsigned target_dim)
{
	assert(D.rows() == D.cols());
	assert(D.rows() >= target_dim);

	const auto n = D.rows();
	const auto ones = Eigen::VectorXd::Ones(n);
	const auto I = Eigen::MatrixXd::Identity(n, n);

	const auto H = I - (1.0 / static_cast<double>(n)) * ones * ones.transpose();
	const auto K = -0.5 * H * D.cwiseAbs2() * H;

	const Eigen::EigenSolver<Eigen::MatrixXd> solver(K);

	Eigen::VectorXd S = solver.eigenvalues().real();
	Eigen::MatrixXd V = solver.eigenvectors().real();

	ExtractNLargestEigens(target_dim, S, V);

	const Eigen::MatrixXd X = Eigen::DiagonalMatrix<double, Eigen::Dynamic>(S.cwiseSqrt()) * V.transpose();

	return X;
}

Plane generate_symmetry_plane_dividing_classical_MDS(Mesh* mesh)
{
	// 1 - create classical MDS
	// 1.1 - generate NXN distance matrix
	int N = mesh->vertices.size();
	Eigen::MatrixXd delta(N, N);

	// 1.2generate geodesic distances matrix ( delta)
	for (size_t i = 0; i < N; i++)
	{
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, i);
		for (size_t j = 0; j < N; j++)
		{
			delta(i, j) = distances[j];
		}
	}
	//2 MDS 
	Eigen::MatrixXd mds_matrix = ComputeClassicalMds(delta, 3);

	// 3 give plane to dominant symmetry plane
	Plane plane = generate_dominant_symmetry_plane(*mesh);
	return plane; 
}


//computes landmark MDS and returns the embedded mesh
Mesh compute_landmark_MDS(Mesh* mesh, const unsigned target_dim, const int no_of_landmarks )
{
	// 1 - get fps wit point determination
	std::vector<unsigned int> landmark_vertex_indices = furthest_point_sampling(mesh, no_of_landmarks);

	// 2 - create landmark matrix 
	Eigen::MatrixXd landmark_delta(no_of_landmarks, no_of_landmarks);
	for (size_t i = 0; i < no_of_landmarks; i++)
	{
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, landmark_vertex_indices[i]);
		for (size_t j = 0; j < no_of_landmarks; j++)
		{
			landmark_delta(i, j) = distances[landmark_vertex_indices[j]] * distances[landmark_vertex_indices[j]];
		}
	}

	Eigen::MatrixXd H = Eigen::MatrixXd::Identity(no_of_landmarks, no_of_landmarks) - (1.0 / no_of_landmarks) * Eigen::VectorXd::Ones(no_of_landmarks)
		* Eigen::VectorXd::Ones(no_of_landmarks
		).transpose();
	Eigen::MatrixXd B = -0.5 * H * landmark_delta * H;


	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(B);
	Eigen::VectorXd S = eigensolver.eigenvalues();
	Eigen::MatrixXd V = eigensolver.eigenvectors();

	ExtractNLargestEigens(target_dim, S, V);

	Eigen::MatrixXd L(target_dim, no_of_landmarks);
	for (size_t i = 0; i < target_dim; i++)
	{
		L.row(i) = sqrt(S(i)) * V.col(i).transpose();

	}
	Mesh landmark_mesh = *mesh;
	landmark_mesh.vertices.clear();
	for (size_t i = 0; i < no_of_landmarks; i++)
	{
		//landmark_mesh.vertices.push_back(glm::vec3(L(0,i) , L(1, i) , L(2, i)) );
		landmark_mesh.vertices.push_back(mesh->vertices[landmark_vertex_indices[i]]);
	}

	// distance based triangulation

	//find n x n matrix
	Eigen::MatrixXd delta_n(no_of_landmarks, no_of_landmarks);
	for (size_t i = 0; i < no_of_landmarks; i++)
	{
		for (size_t j = 0; j < no_of_landmarks; j++)
		{
			double distance = glm::distance(landmark_mesh.vertices[i], landmark_mesh.vertices[j]);
			delta_n(i, j) = distance * distance;
		}
	}
	//compute mean 
	Eigen::VectorXd delta_n_mean(no_of_landmarks);
	for (size_t i = 0; i < no_of_landmarks; i++)
	{
		delta_n_mean(i) = 0;
	}
	for (size_t i = 0; i < no_of_landmarks; i++)
	{
		delta_n_mean = delta_n_mean + delta_n.col(i);
	}
	delta_n_mean /= no_of_landmarks;

	//pseudo inverse transpose
	Eigen::MatrixXd L_k(target_dim, no_of_landmarks);
	for (size_t i = 0; i < target_dim; i++)
	{
		L_k.row(i) = V.col(i).transpose() / sqrt(S(i));
	}
	std::cout << " V dim " << V.rows() << " " << V.cols() << std::endl;

	std::vector<glm::vec3> temp_vertices;
	for (size_t i = 0; i < mesh->vertices.size(); i++)
	{
		Eigen::VectorXd delta_a(no_of_landmarks);
		bool is_point_landmark = false;
		size_t j = 0;
		for (j = 0; j < no_of_landmarks; j++)
		{
			float dist = glm::distance(mesh->vertices[i], landmark_mesh.vertices[j]);
			delta_a(j) = dist * dist;
		}
		Eigen::VectorXd x_a_eigen = -1.0 / 2.0 * L_k * (delta_a - delta_n_mean);
		glm::vec3 x_a(x_a_eigen(0), x_a_eigen(1), x_a_eigen(2));
		temp_vertices.push_back(x_a);
	}
	landmark_mesh.vertices = temp_vertices;
	//get the center point from mesh

	return landmark_mesh; 
}
Plane trilateral_symmetry_with_landmark_MDS_with_plane(Mesh* mesh, const unsigned target_dim, const int no_of_landmarks , const int no_of_trilateral_points)
{
	Mesh L_MDS_mesh = compute_landmark_MDS(mesh, target_dim);
	//calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < L_MDS_mesh.vertices.size(); i++)
	{
		plane_center += L_MDS_mesh.vertices[i];
	}
	plane_center /= mesh->vertices.size();
	Plane plane = generate_dominant_symmetry_plane(plane_center , L_MDS_mesh);

	std::vector<unsigned int> points_plane_positive;
	std::vector<unsigned int> points_plane_negative;
	//now separate the points into two sides of the plane 
	for (size_t i = 0; i < L_MDS_mesh.vertices.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &L_MDS_mesh.vertices[i]) >= 0 )
		{
			points_plane_positive.push_back(i);
		}
		else
		{
			points_plane_negative.push_back(i);
		}
	}
	// now do two distinct fps
	std::vector<unsigned int > fps_positive =  furthest_point_sampling_on_partial_points(&L_MDS_mesh, no_of_trilateral_points, points_plane_positive);
	std::vector<unsigned int > fps_negative =  furthest_point_sampling_on_partial_points(&L_MDS_mesh, no_of_trilateral_points, points_plane_negative);

	/*std::vector<unsigned int > fps_points =  furthest_point_sampling(&L_MDS_mesh, no_of_trilateral_points * 2 , true  );
	fps_positive.clear();
	fps_negative.clear();
	for (size_t i = 0; i < fps_points.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &L_MDS_mesh.vertices[i]) >= 0)
		{
			fps_positive.push_back(i);
		}
		else
		{
			fps_negative.push_back(i);
		}
	}*/

	// trilateral computation
	std::vector<TrilateralDescriptor> positive_mesh_trilateral_descriptor = get_trilateral_points_using_closest_pairs(&L_MDS_mesh ,fps_positive);
	std::vector<TrilateralDescriptor> negative_mesh_trilateral_descriptor = get_trilateral_points_using_closest_pairs(&L_MDS_mesh, fps_negative);

	// write a function for comparing two descriptor
	//irrelevant constants 
	float const1 = 0;
	float const2 = 0;
	float const3 = 0;
	std::vector<std::pair<unsigned int , unsigned int>> resemblance_pairs =  point_match_trilateral_weights(&L_MDS_mesh, positive_mesh_trilateral_descriptor, negative_mesh_trilateral_descriptor, const1, const2, const3);
	
	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = get_geodesic_cost_with_list(&L_MDS_mesh, left_correspondences, right_correspondences);
	
	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < fps_positive.size(); i++)
	{
		std::vector<float> distances  = compute_geodesic_distances_fibonacci_heap_distances(*mesh, fps_positive[i]);
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (maximum_geodesic_distance < distances[j])
			{
				maximum_geodesic_distance = distances[j];
			}
		}
	}

	// color left red
	std::vector<unsigned int> is_selected(mesh->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		mesh->colors[resemblance_pairs[i].first].r = 255;
		mesh->colors[resemblance_pairs[i].first].g = 0;
		mesh->colors[resemblance_pairs[i].first].b = 0;

		mesh->colors[resemblance_pairs[i].second].r = 0;
		mesh->colors[resemblance_pairs[i].second].g = 0;
		mesh->colors[resemblance_pairs[i].second].b = 255;
	}
	//L_MDS_mesh.colors = mesh->colors;
	//*mesh = L_MDS_mesh;
	//color right  blue 

	/*L_MDS_mesh.colors.clear();
	for (size_t i = 0; i < L_MDS_mesh.vertices.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &L_MDS_mesh.vertices[i]) >= 0)
		{
			L_MDS_mesh.colors.push_back(glm::vec3(255.0 , 0.0 , 0.0 ) );
		}
		else
		{
			L_MDS_mesh.colors.push_back(glm::vec3(0, 255, 0.0));
		}
	}
	*mesh = L_MDS_mesh;*/

	std::cout << " total average error is " << total_error << " maximum geodesic distance is " << maximum_geodesic_distance << std::endl; 
	return plane;
}


