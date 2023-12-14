#include "../Include/NLateralDescriptor.h"
#include "../Include/TrilateralMap.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#include "../Include/MetricCalculations.h"
#pragma region Nlateral struct

NLateralDescriptor::NLateralDescriptor(Mesh& mesh, const std::vector<unsigned int>& point_indices, int N)
{
	this->point_indices = point_indices;
	this->N = N;
	int point_size = point_indices.size();
	this->mesh = &mesh;
}

void NLateralDescriptor::get_euclidian_distances()
{
	//euclidian
	int point_size = this->point_indices.size();
	for (size_t i = 0; i < point_size; i++)
	{
		this->euclidian_distances.push_back(std::vector<double>());
		for (size_t j = 0; j < point_size; j++)
		{
			if (i == j)
			{
				continue;
			}
			float euclidian_dist = 0;
			euclidian_dist = glm::distance(mesh->vertices[point_indices[i]], mesh->vertices[point_indices[j]]);
			this->euclidian_distances[i].push_back(euclidian_dist);
		}
	}
}
void NLateralDescriptor::get_geodesic_distances()
{
	//geodesic
	int point_size = this->point_indices.size();
	for (size_t i = 0; i < point_size; i++)
	{
		this->geodesic_distances.push_back(std::vector<double>());
		std::vector<float> geodesic_dist = compute_geodesic_distances_min_heap_distances((Mesh&)mesh, this->point_indices[i]);

		for (size_t j = 0; j < point_size; j++)
		{
			if (i == j)
			{
				continue;
			}
			this->geodesic_distances[i].push_back(geodesic_dist[this->point_indices[j]]);
		}
	}
}
void NLateralDescriptor::get_curvatures()
{
	//curvature
	int point_size = this->point_indices.size();
	for (size_t i = 0; i < point_size; i++)
	{
		curvatures.push_back(std::vector<double>());

		for (size_t j = 0; j < point_size; j++)
		{
			if (i == j)
			{
				continue;
			}
			else
			{
				curvatures[i].push_back(this->euclidian_distances[i][j] / this->geodesic_distances[i][j]);
			}
		}
	}
}
void NLateralDescriptor::get_k_ring_areas(int K)
{
	//k-ring-area
	int point_size = this->point_indices.size();
	for (size_t i = 0; i < point_size; i++)
	{
		float k_ring_area = get_N_ring_area(mesh, this->point_indices[i], K);
		this->k_ring_areas.push_back(k_ring_area);
	}
}

#pragma endregion



NLateralDescriptor generate_NLateralDescriptor(Mesh* m, const std::vector<unsigned int>& mesh_indices  , const std::vector<bool>& parameter_checkbox
	, const std::vector<float>& parameter_weights, const std::vector<std::string>& parameter_names)
{
	NLateralDescriptor nlateralDescriptor( *m ,  mesh_indices , mesh_indices.size());

	for (size_t i = 0; i < parameter_checkbox.size(); i++)
	{
		if (parameter_checkbox[i]) //if true
		{
			if (parameter_names[i].find("area") != std::string::npos)
			{
				//nlateralDescriptor.get
			}
			else if (parameter_names[i].find("euclidian") != std::string::npos)
			{
				nlateralDescriptor.get_euclidian_distances();
			}
			else if (parameter_names[i].find("geodesic") != std::string::npos)
			{
				nlateralDescriptor.get_geodesic_distances();
			}
			else if (parameter_names[i].find("ring") != std::string::npos)
			{
				// get n ring number 
				std::vector<int> num_vec = getNumberFromString(parameter_names[i]);
				nlateralDescriptor.get_k_ring_areas(num_vec[0]);
			}
			else if (parameter_names[i].find("curvature") != std::string::npos)
			{
				// get n ring number 
				nlateralDescriptor.get_curvatures();
			}
		}
	}
	return nlateralDescriptor;
}

std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_furthest_pairs(Mesh* m, std::vector<unsigned int>& indices, int N  )
{
	std::vector<NLateralDescriptor> nLateralDescVec;
	for (size_t i = 0; i < indices.size(); i++)
	{
		//get n-1 of the furthest indexed points
		std::vector<float> geodesic_distances = compute_geodesic_distances_fibonacci_heap_distances(*m, indices[i]);
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
		
		//get n-1 furthest
		std::vector<bool> is_already_selected(  m->vertices.size() , false);
		std::vector<unsigned int> selected_indices;
		is_already_selected[indices[i]] = true; // itself is selected 
		selected_indices.push_back(indices[i]);
		for (size_t j = 0; j < N-1; j++)
		{
			float maxVal = -INFINITY;
			float maxIndexFirst = -1;
			for (size_t k = 0; k < distances.size(); k++)
			{
				if (maxVal < distances[k].first && !is_already_selected[k])
				{
					maxIndexFirst = distances[k].second;
					maxVal = distances[k].first;
				}
			}
			selected_indices.push_back(maxIndexFirst);
		}
		
		NLateralDescriptor desc = generate_NLateralDescriptor(m, selected_indices, N_LATERAL_PARAMETERS.parameter_checkbox,
		N_LATERAL_PARAMETERS.parameter_weights, N_LATERAL_PARAMETERS.parameter_names);
		nLateralDescVec.push_back(desc);
	}
	return nLateralDescVec;
}
std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_closest_pairs(Mesh* m, std::vector<unsigned int>& indices, int N)
{
	std::vector<NLateralDescriptor> nLateralDescVec;
	for (size_t i = 0; i < indices.size(); i++)
	{
		//get n-1 of the closest indexed points
		std::vector<float> geodesic_distances = compute_geodesic_distances_fibonacci_heap_distances(*m, indices[i]);
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

		//get n-1 closest
		std::vector<bool> is_already_selected(m->vertices.size(), false);
		std::vector<unsigned int> selected_indices;
		is_already_selected[indices[i]] = true; // itself is selected 
		selected_indices.push_back(indices[i]);
		for (size_t j = 0; j < N - 1; j++)
		{
			float maxVal = INFINITY;
			float maxIndexFirst = -1;
			for (size_t k = 0; k < distances.size(); k++)
			{
				if (maxVal > distances[k].first && !is_already_selected[k])
				{
					maxIndexFirst = distances[k].second;
					maxVal = distances[k].first;
				}
			}
			selected_indices.push_back(maxIndexFirst);
		}

		NLateralDescriptor desc = generate_NLateralDescriptor(m, selected_indices, N_LATERAL_PARAMETERS.parameter_checkbox,
			N_LATERAL_PARAMETERS.parameter_weights, N_LATERAL_PARAMETERS.parameter_names);
		nLateralDescVec.push_back(desc);
	}
	return nLateralDescVec;
}
NLateralParameters::NLateralParameters()
{
	this->parameter_weights = std::vector<float>(this->NO_OF_PARAMETERS);
	this->parameter_names = std::vector<std::string>(this->NO_OF_PARAMETERS);
	this->parameter_checkbox = std::vector<bool>(this->NO_OF_PARAMETERS);
	this->n_lateral_construction_methods = std::vector<std::string>(this->N_LATERAL_CONSTRUCTION_METHOD_NO);

	parameter_names[0] = "area";
	parameter_names[1] = "euclidian distance";
	parameter_names[2] = "geodesic distance";
	parameter_names[3] = "curvature";
	parameter_names[4] = "Heat Kernel Signature";
	parameter_names[this->K_RING_POS] = "k ring area = ";
	parameter_names[6] = "X";
	parameter_names[7] = "X";
	parameter_names[8] = "X";


	n_lateral_construction_methods[0] = "closest points";
	n_lateral_construction_methods[1] = "furthest points";

	//init arrays
	for (size_t i = 0; i < this->NO_OF_PARAMETERS; i++)
	{
		parameter_checkbox[i] = false;
	}
	for (size_t i = 0; i < this->NO_OF_PARAMETERS; i++)
	{
		parameter_weights[i] = 1.0;
	}
}


void start_n_lateral_algorithm(Mesh* mesh)
{

	Mesh L_MDS_mesh = compute_landmark_MDS(mesh, 3); // 3 is as always 
	//calculate center of the plane 
	glm::vec3 plane_center(0, 0, 0);
	for (size_t i = 0; i < L_MDS_mesh.vertices.size(); i++)
	{
		plane_center += L_MDS_mesh.vertices[i];
	}
	plane_center /= mesh->vertices.size();
	Plane plane = generate_dominant_symmetry_plane(plane_center, L_MDS_mesh);

	std::vector<unsigned int> points_plane_positive;
	std::vector<unsigned int> points_plane_negative;
	//now separate the points into two sides of the plane 
	for (size_t i = 0; i < L_MDS_mesh.vertices.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &L_MDS_mesh.vertices[i]) >= 0)
		{
			points_plane_positive.push_back(i);
		}
		else
		{
			points_plane_negative.push_back(i);
		}
	}
	// now do two distinct fps
	std::vector<unsigned int > fps_positive = furthest_point_sampling_on_partial_points(&L_MDS_mesh, N_LATERAL_PARAMETERS.no_of_N_lateral_pairs, points_plane_positive);
	std::vector<unsigned int > fps_negative = furthest_point_sampling_on_partial_points(&L_MDS_mesh, N_LATERAL_PARAMETERS.no_of_N_lateral_pairs, points_plane_negative);

	// trilateral computation
	std::vector<NLateralDescriptor> positive_mesh_N_lateral_descriptor;
	std::vector<NLateralDescriptor> negative_mesh_N_lateral_descriptor;
	if (N_LATERAL_PARAMETERS.current_n_lateral_construction_method.find("closest") != std::string::npos)
	{
		positive_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_closest_pairs(&L_MDS_mesh, fps_positive, N_LATERAL_PARAMETERS.N);
		negative_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_closest_pairs(&L_MDS_mesh, fps_positive, N_LATERAL_PARAMETERS.N);
	}
	else if (N_LATERAL_PARAMETERS.current_n_lateral_construction_method.find("furthest") != std::string::npos)
	{
		positive_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_furthest_pairs(&L_MDS_mesh, fps_positive, N_LATERAL_PARAMETERS.N);
		negative_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_furthest_pairs(&L_MDS_mesh, fps_positive, N_LATERAL_PARAMETERS.N);
	}

	// write a function for comparing two descriptor
	//irrelevant constants 
	float const1 = 0;
	float const2 = 0;
	float const3 = 0;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = point_match_n_lateral_descriptors(&L_MDS_mesh, positive_mesh_N_lateral_descriptor, negative_mesh_N_lateral_descriptor);

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
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, fps_positive[i]);
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

	mesh->calculated_symmetry_pairs = resemblance_pairs;

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

	std::cout << " total average error is " << total_error << " maximum geodesic distance is " << maximum_geodesic_distance <<  " " << total_error / maximum_geodesic_distance <<
	std::endl;
}

std::vector <std::pair<unsigned int, unsigned int>> point_match_n_lateral_descriptors(Mesh* m, const std::vector<NLateralDescriptor>& nlateral_vec_left, const std::vector<NLateralDescriptor>& n_lateral_vec_right)
{
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;


	//first get the vector size with checkboxed parameters 
	int size_of_vector = 0;
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; i++)
	{
		if (N_LATERAL_PARAMETERS.parameter_checkbox[i])
		{
			if (N_LATERAL_PARAMETERS.parameter_names[i].find("euclidian"))
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("geodesic"))
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("curvature"))
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("ring"))
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("area"))
			{
				size_of_vector += 1;
			}
		}
	}
	// need to get all of the permuations as vectors
	std::vector<unsigned int> permutation_vector(N_LATERAL_PARAMETERS.N - 1);
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.N - 1; i++)
	{
		permutation_vector.push_back(i);
	}

	//all permutations for indices
	std::vector<std::vector<unsigned int>> all_permutations;
	do {
		all_permutations.push_back(permutation_vector);
	} while (std::next_permutation(permutation_vector.begin(), permutation_vector.end()));


	
	for (size_t i = 0; i < nlateral_vec_left.size(); i++)
	{
		//for comparisons in the end
		float smallest_dif = INFINITY;
		std::pair<unsigned int, unsigned int> smallest_pair; // in terms of indices


		NLateralDescriptor desc_i = nlateral_vec_left[i];
		std::vector<Eigen::VectorXd>desc_i_vectors(all_permutations.size()); //size is number of permutations 

		// generate vectors for all of the permutations
		for (size_t j = 0; j < all_permutations.size(); j++)
		{
			desc_i_vectors[j] = Eigen::VectorXd(size_of_vector);
		}
		// fill all of the vectors
		for (size_t j = 0; j < all_permutations.size(); j++)
		{
			int current_size = 0;
			for (size_t k = 0; k < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; k++)
			{
				if (N_LATERAL_PARAMETERS.parameter_checkbox[k])
				{
					if (N_LATERAL_PARAMETERS.parameter_names[k].find("euclidian"))
					{
						for (int t = 0; t < desc_i.euclidian_distances.size(); t++)
						{
							desc_i_vectors[j](current_size++) = desc_i.euclidian_distances[0][all_permutations[j][t]];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[k].find("geodesic"))
					{
						for (int t = 0; t < desc_i.euclidian_distances.size(); t++)
						{
							desc_i_vectors[j](current_size++) = desc_i.geodesic_distances[0][all_permutations[j][t]];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[k].find("curvature"))
					{
						for (int t = 0; t < desc_i.euclidian_distances.size(); t++)
						{
							desc_i_vectors[j](current_size++) = desc_i.curvatures[0][all_permutations[j][t]];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[k].find("ring"))
					{
						for (int t = 0; t < desc_i.k_ring_areas.size(); t++)
						{
							desc_i_vectors[j](current_size++) = desc_i.k_ring_areas[all_permutations[j][t]];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[k].find("area"))
					{
						desc_i_vectors[j](current_size++) = desc_i.area;
					}
				}
			}
		
		}

		for (size_t j = 0; j < n_lateral_vec_right.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			NLateralDescriptor desc_j = n_lateral_vec_right[j];

			std::vector<Eigen::VectorXd>desc_j_vectors(all_permutations.size()); //size is number of permutations 

			// generate vectors for all of the permutations
			for (size_t k = 0; k < all_permutations.size(); k++)
			{
				desc_j_vectors[k] = Eigen::VectorXd(size_of_vector);
			}


			for (size_t k = 0; k < all_permutations.size(); k++)
			{
				int current_size = 0;
				for (size_t t = 0; t < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; t++)
				{
					if (N_LATERAL_PARAMETERS.parameter_checkbox[t])
					{
						if (N_LATERAL_PARAMETERS.parameter_names[t].find("euclidian"))
						{
							for (int t = 0; t < desc_i.euclidian_distances.size(); t++)
							{
								desc_j_vectors[k](current_size++) = desc_i.euclidian_distances[0][all_permutations[k][t]];
							}
						}
						else if (N_LATERAL_PARAMETERS.parameter_names[t].find("geodesic"))
						{
							for (int t = 0; t < desc_i.euclidian_distances.size(); t++)
							{
								desc_j_vectors[k](current_size++) = desc_i.geodesic_distances[0][all_permutations[k][t]];
							}
						}
						else if (N_LATERAL_PARAMETERS.parameter_names[t].find("curvature"))
						{
							for (int t = 0; t < desc_i.euclidian_distances.size(); t++)
							{
								desc_j_vectors[k](current_size++) = desc_i.curvatures[0][all_permutations[k][t]];
							}
						}
						else if (N_LATERAL_PARAMETERS.parameter_names[t].find("ring"))
						{
							for (int t = 0; t < desc_i.k_ring_areas.size(); t++)
							{
								desc_j_vectors[k](current_size++) = desc_i.k_ring_areas[all_permutations[k][t]];
							}
						}
						else if (N_LATERAL_PARAMETERS.parameter_names[t].find("area"))
						{
							desc_j_vectors[k](current_size++) = desc_i.area;
						}
					}
				}

			}


			
			//compare all
			for (size_t k = 0; k < all_permutations.size(); k++)
			{
				for (size_t t = 0; t < all_permutations.size(); t++)
				{
					float dist = (desc_i_vectors[k] - desc_j_vectors[t]).norm();
					if (dist < smallest_dif)
					{
						dist = smallest_dif;
						smallest_pair.first = i;
						smallest_pair.second = j;
					}
				}
			}

		}

		resemblance_pairs.push_back(std::pair<unsigned int, unsigned int >(nlateral_vec_left[smallest_pair.first].point_indices[0], nlateral_vec_left[smallest_pair.second].point_indices[0]));

	}
	return resemblance_pairs;
}