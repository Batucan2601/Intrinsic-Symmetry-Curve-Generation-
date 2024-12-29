#include "../Include/NLateralDescriptor.h"
#include "../Include/TrilateralMap.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#include "../Include/MetricCalculations.h"
#include "../Include/Geodesic.h"
#include "../Include/ROI.h""
#include "../Include/ShapeDiameter.h""
#pragma region Nlateral struct

NLateralDescriptor::NLateralDescriptor(TrilateralMesh& mesh, const std::vector<unsigned int>& point_indices, int N)
{
	this->point_indices = point_indices;
	this->N = N;
	int point_size = point_indices.size();
	this->mesh = mesh;
}

NLateralDescriptor::NLateralDescriptor()
{
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
			
			float euclidian_dist = 0;
			euclidian_dist = glm::distance(mesh.vertices[point_indices[i]], mesh.vertices[point_indices[j]]);
			if (i == j)
			{
				euclidian_dist = 0;
			}
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
		std::vector<float> geodesic_dist = compute_geodesic_distances_min_heap_distances(mesh, this->point_indices[i]);

		for (size_t j = 0; j < point_size; j++)
		{
			if (i == j)
			{
				this->geodesic_distances[i].push_back(0);

			}
			else
			{
				this->geodesic_distances[i].push_back(geodesic_dist[this->point_indices[j]]);
			}
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
				curvatures[i].push_back(0);
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
		float k_ring_area = get_N_ring_area(&mesh, this->point_indices[i], K);
		this->k_ring_areas.push_back(k_ring_area);
	}
}

#pragma endregion



NLateralDescriptor generate_NLateralDescriptor(TrilateralMesh* m, const std::vector<unsigned int>& mesh_indices  , const std::vector<bool>& parameter_checkbox
	, const std::vector<float>& parameter_weights, const std::vector<std::string>& parameter_names)
{
	NLateralDescriptor nlateralDescriptor( *m ,  mesh_indices , mesh_indices.size());

	for (size_t i = 0; i < parameter_checkbox.size(); i++)
	{
		if (parameter_checkbox[i]) //if true
		{
			if (parameter_names[i].find("area") != std::string::npos)
			{
				nlateralDescriptor.get_ROI();
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

std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_furthest_pairs(TrilateralMesh* m, std::vector<unsigned int>& indices
, NLateralParameters N_LATERAL_PARAMETERS)
{
	int N = N_LATERAL_PARAMETERS.N;
	std::vector<NLateralDescriptor> nLateralDescVec;
	for (size_t i = 0; i < indices.size(); i++)
	{
		//get n-1 of the furthest indexed points
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
				if (maxVal < distances[k].first && !is_already_selected[distances[k].second])
				{
					maxIndexFirst = distances[k].second;
					maxVal = distances[k].first;
				}
			}
			selected_indices.push_back(maxIndexFirst);
			is_already_selected[(int)maxIndexFirst] = true;
		}
		
		NLateralDescriptor desc = generate_NLateralDescriptor(m, selected_indices, N_LATERAL_PARAMETERS.parameter_checkbox,
		N_LATERAL_PARAMETERS.parameter_weights, N_LATERAL_PARAMETERS.parameter_names);
		nLateralDescVec.push_back(desc);
	}
	return nLateralDescVec;
}
std::vector<NLateralDescriptor> get_N_lateral_descriptor_using_closest_pairs(TrilateralMesh* m, std::vector<unsigned int>& indices,  NLateralParameters N_LATERAL_PARAMETERS )
{
	int N = N_LATERAL_PARAMETERS.N;
	std::vector<NLateralDescriptor> nLateralDescVec;
	for (size_t i = 0; i < indices.size(); i++)
	{
		//get n-1 of the closest indexed points
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
				if (maxVal > distances[k].first && !is_already_selected[ distances[k].second])
				{
					maxIndexFirst = distances[k].second;
					maxVal = distances[k].first;
				}
			}
			selected_indices.push_back(maxIndexFirst);
			is_already_selected[(int)maxIndexFirst] = true;
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


void start_n_lateral_algorithm(MeshFactory& mesh_fac, int selected_mesh, NLateralParameters N_LATERAL_PARAMETERS)
{
	TrilateralMesh* mesh = &mesh_fac.mesh_vec[selected_mesh];

	TrilateralMesh L_MDS_mesh = compute_landmark_MDS(mesh, 3); // 3 is as always 
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
		positive_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_closest_pairs(&L_MDS_mesh, fps_positive, N_LATERAL_PARAMETERS);
		negative_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_closest_pairs(&L_MDS_mesh, fps_negative, N_LATERAL_PARAMETERS);
	}
	else if (N_LATERAL_PARAMETERS.current_n_lateral_construction_method.find("furthest") != std::string::npos)
	{
		positive_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_furthest_pairs(&L_MDS_mesh, fps_positive, N_LATERAL_PARAMETERS);
		negative_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_furthest_pairs(&L_MDS_mesh, fps_negative, N_LATERAL_PARAMETERS);
	}

	// calculate the maximums
	NLateral_parameters_calculate_maximums(&L_MDS_mesh, N_LATERAL_PARAMETERS, fps_positive, fps_negative);


	// write a function for comparing two descriptor
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = point_match_n_lateral_descriptors(&L_MDS_mesh, positive_mesh_N_lateral_descriptor, negative_mesh_N_lateral_descriptor
	, N_LATERAL_PARAMETERS);

	//forge it into two list
	std::vector<unsigned int> left_correspondences;
	std::vector<unsigned int> right_correspondences;
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		left_correspondences.push_back(resemblance_pairs[i].first);
		right_correspondences.push_back(resemblance_pairs[i].second);
	}
	float total_error = Metric_get_geodesic_cost_with_list(&L_MDS_mesh, left_correspondences, right_correspondences);

	// now use fps points to get maximum distance in order to compare to 
	float maximum_geodesic_distance = 0;
	for (size_t i = 0; i < fps_positive.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*mesh, fps_positive[i]);
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

	// file creating
	ofstream txtFile( "C:/Users/BATU/source/repos/Trilateral/Results/" + mesh->file_name + ".txt");
	// 
	txtFile << "===============================================================================" << std::endl;
	txtFile << " N lateral " << std::endl;
	//first write the N_lateral parameters
	txtFile << " N ===== " << N_LATERAL_PARAMETERS.N << std::endl;
	txtFile << "number of sampled pairs  " << N_LATERAL_PARAMETERS.no_of_N_lateral_pairs  <<   "  " << N_LATERAL_PARAMETERS.no_of_N_lateral_pairs * 2.0  / 
	mesh->vertices.size() << std::endl;
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; i++)
	{
		if (N_LATERAL_PARAMETERS.parameter_checkbox[i])
		{
			txtFile << N_LATERAL_PARAMETERS.parameter_names[i] << "  " << N_LATERAL_PARAMETERS.parameter_weights[i] << std::endl; 
		}
	}
	txtFile << "===============================================================================" << std::endl;
	float error_percentage = total_error / maximum_geodesic_distance;
	txtFile <<  " geodesic error " + std::to_string(error_percentage) + "\n";
	txtFile.close();

}

std::vector <std::pair<unsigned int, unsigned int>> point_match_n_lateral_descriptors(TrilateralMesh* m, const std::vector<NLateralDescriptor>& nlateral_vec_left, const std::vector<NLateralDescriptor>& n_lateral_vec_right,
	NLateralParameters N_LATERAL_PARAMETERS)
{
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	std::vector<std::vector<Eigen::VectorXd>> n_lateral_right_vectors;
	//first get the vector size with checkboxed parameters 
	int size_of_vector = 0;
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; i++)
	{
		if (N_LATERAL_PARAMETERS.parameter_checkbox[i])
		{
			if (N_LATERAL_PARAMETERS.parameter_names[i].find("euclidian") != std::string::npos)
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("geodesic") != std::string::npos)
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("curvature") != std::string::npos)
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("ring") != std::string::npos)
			{
				size_of_vector += N_LATERAL_PARAMETERS.N;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("area") != std::string::npos )
			{
				size_of_vector += 1;
			}
		}
	}
	
	for (size_t i = 0; i < nlateral_vec_left.size(); i++)
	{
		NLateralDescriptor desc_i = nlateral_vec_left[i];
		float minimum_dist = INFINITY;

		Eigen::VectorXf vector_i(size_of_vector);
		int dynamic_vec_size = 0;
		int closest_index_j = -1;

		for (size_t j = 0; j < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; j++)
		{
			if (N_LATERAL_PARAMETERS.parameter_checkbox[j])
			{
				if (N_LATERAL_PARAMETERS.parameter_names[j].find("euclidian") != std::string::npos)
				{
					for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
					{
						vector_i(dynamic_vec_size++) = desc_i.euclidian_distances[0][k];
					}
				}
				else if (N_LATERAL_PARAMETERS.parameter_names[j].find("geodesic") != std::string::npos)
				{
					for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
					{
						vector_i(dynamic_vec_size++) = desc_i.geodesic_distances[0][k];
					}
				}
				else if (N_LATERAL_PARAMETERS.parameter_names[j].find("curvature") != std::string::npos)
				{
					for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
					{
						vector_i(dynamic_vec_size++) = desc_i.euclidian_distances[0][k];
					}
				}
				else if (N_LATERAL_PARAMETERS.parameter_names[j].find("ring") != std::string::npos)
				{
					for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
					{
						vector_i(dynamic_vec_size++) = desc_i.euclidian_distances[0][k];
					}
				}
				else if (N_LATERAL_PARAMETERS.parameter_names[j].find("area") != std::string::npos)
				{
						vector_i(dynamic_vec_size++) = desc_i.area;
				}
			}
		}


		for (size_t j = 0; j < n_lateral_vec_right.size(); j++)
		{
			NLateralDescriptor desc_j = n_lateral_vec_right[j];
			Eigen::VectorXf vector_j(size_of_vector);
			int dynamic_vec_size = 0;
			for (size_t j = 0; j < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; j++)
			{
				if (N_LATERAL_PARAMETERS.parameter_checkbox[j])
				{
					if (N_LATERAL_PARAMETERS.parameter_names[j].find("euclidian") != std::string::npos)
					{
						for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
						{
							vector_j(dynamic_vec_size++) = desc_j.euclidian_distances[0][k];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[j].find("geodesic") != std::string::npos)
					{
						for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
						{
							vector_j(dynamic_vec_size++) = desc_j.geodesic_distances[0][k];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[j].find("curvature") != std::string::npos)
					{
						for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
						{
							vector_j(dynamic_vec_size++) = desc_j.euclidian_distances[0][k];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[j].find("ring") != std::string::npos)
					{
						for (size_t k = 0; k < N_LATERAL_PARAMETERS.N; k++)
						{
							vector_j(dynamic_vec_size++) = desc_j.euclidian_distances[0][k];
						}
					}
					else if (N_LATERAL_PARAMETERS.parameter_names[j].find("area") != std::string::npos)
					{
						vector_j(dynamic_vec_size++) = desc_j.area;
					}
				}
			}


			Eigen::VectorXf dif = vector_i - vector_j;
			float norm = dif.norm();
			if (minimum_dist > norm)
			{
				minimum_dist = norm;
				closest_index_j = j;
			}

		}

		resemblance_pairs.push_back(std::pair<unsigned int, unsigned int >(nlateral_vec_left[i].point_indices[0], n_lateral_vec_right[closest_index_j].point_indices[0]));
	}

	return resemblance_pairs;
}

void NLateral_parameters_calculate_maximums(TrilateralMesh* m, NLateralParameters& N_LATERAL_PARAMETERS , std::vector<unsigned int>& left , std::vector<unsigned int>& right)
{
	// for each parameter calculate the maximum in order to normalize the parameters to give them meaningfull weights
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; i++)
	{
		if (N_LATERAL_PARAMETERS.parameter_checkbox[i])
		{
			if (N_LATERAL_PARAMETERS.parameter_names[i].find("euclidian") != std::string::npos )
			{
				float maximum_dist = -INFINITY; 
				//calculate the maximum euclidian distances
				for (size_t j = 0; j < left.size(); j++)
				{
					for (size_t k = 0; k < right.size(); k++)
					{
						if (maximum_dist < glm::distance(m->vertices[left[i]], m->vertices[right[j]]))
						{
							maximum_dist = glm::distance(m->vertices[left[i]], m->vertices[right[j]]);
						}
					}
				}
				N_LATERAL_PARAMETERS.parameter_maximums["euclidian"] = maximum_dist;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("geodesic") != std::string::npos )
			{
				float maximum_dist = -INFINITY;
				//calculate the maximum euclidian distances
				std::vector<float> distances = Geodesic_dijkstra(*m, left[i]);
				for (size_t j = 0; j < distances.size(); j++)
				{
					if (maximum_dist < distances[j])
					{
						maximum_dist = distances[j];
					}
				}
				N_LATERAL_PARAMETERS.parameter_maximums["geodesic"] = maximum_dist;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("curvature") != std::string::npos )
			{
				N_LATERAL_PARAMETERS.parameter_maximums["curvature"] = 1.0;
			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("ring") != std::string::npos)
			{
				float maximum_area = -INFINITY;
				for (size_t j = 0; j < left.size(); j++)
				{
					float area = get_N_ring_area(m, left[j], N_LATERAL_PARAMETERS.N_RING_NO);
					if (maximum_area < area)
					{
						maximum_area = area; 
					}
				}
				for (size_t j = 0; j < right.size(); j++)
				{
					float area = get_N_ring_area(m, right[j], N_LATERAL_PARAMETERS.N_RING_NO);
					if (maximum_area < area)
					{
						maximum_area = area;
					}
				}
				N_LATERAL_PARAMETERS.parameter_maximums["ring"] = 1.0;

			}
			else if (N_LATERAL_PARAMETERS.parameter_names[i].find("area") != std::string::npos)
			{
				float maximum_area = -INFINITY;
				
				maximum_area = 1.0f;
			}
		}
	}
}
// lets try connecting to n-closest endpoint for now?
void start_n_lateral_algorithm_for_mesh(Skeleton& skeleton, NLateralParameters N_LATERAL_PARAMETERS)
{
	int N = skeleton.skeletonFormat.size();
	std::vector<unsigned int> end_point_indices;

	for (size_t i = 0; i < N; i++)
	{
		if (skeleton.skeletonFormat[i].label == END)
		{
			end_point_indices.push_back(i);
		}
	}
	int end_point_size = end_point_indices.size();
	//calcualte distances between 
	std::vector<std::vector<float>>distances_between_end_points(end_point_size);
	for ( int i = 0; i < end_point_size; i++)
	{
		std::vector<float> distances_i(N);
		
		//irrelevant for this operation
		float dist = 0;
		std::vector<unsigned int> vertex_list; 
		//irrelevant for this operation

		//skeleton_calculate_distances_and_vertex_list(skeleton, end_point_indices[i], 0, dist, vertex_list, distances_i);
		
		std::vector<float> distances_between_end_points_i;
		for (size_t j = 0; j < end_point_size; j++)
		{
			distances_between_end_points_i.push_back(end_point_indices[j]);
		}

	}


	//step 2: generate trilateral descriptors
	// need to get all of the permuations as vectors
	std::vector<unsigned int> permutation_vector;
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.N; i++)
	{
		permutation_vector.push_back(i);
	}
	//all permutations for indices
	std::vector<std::vector<unsigned int>> all_permutations;
	do {
		all_permutations.push_back(permutation_vector);
	} while (std::next_permutation(permutation_vector.begin(), permutation_vector.end()));


	/*for (size_t i = 0; i < length; i++)
	{

	}*/

}


void start_n_lateral_algorithm_with_skeleton_end_points(TrilateralMesh* m, NLateralParameters& N_LATERAL_PARAMETERS,
std::vector<unsigned int>&mesh_left_endpoints, std::vector<unsigned int>&mesh_right_endpoints)
{
	// trilateral computation
	std::vector<NLateralDescriptor> positive_mesh_N_lateral_descriptor;
	std::vector<NLateralDescriptor> negative_mesh_N_lateral_descriptor;
	if (N_LATERAL_PARAMETERS.current_n_lateral_construction_method.find("closest") != std::string::npos)
	{
		positive_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_closest_pairs(m, mesh_left_endpoints, N_LATERAL_PARAMETERS);
		negative_mesh_N_lateral_descriptor = get_N_lateral_descriptor_using_closest_pairs(m, mesh_right_endpoints, N_LATERAL_PARAMETERS);
	}
	NLateral_parameters_calculate_maximums(m, N_LATERAL_PARAMETERS, mesh_left_endpoints, mesh_right_endpoints);

	// write a function for comparing two descriptor
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs = point_match_n_lateral_descriptors(m, positive_mesh_N_lateral_descriptor, negative_mesh_N_lateral_descriptor
		, N_LATERAL_PARAMETERS);

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
	for (size_t i = 0; i < mesh_right_endpoints.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, mesh_right_endpoints[i]);
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
	/*for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->colors[resemblance_pairs[i].first].r = 255;
		m->colors[resemblance_pairs[i].first].g = 0;
		m->colors[resemblance_pairs[i].first].b = 0;

		m->colors[resemblance_pairs[i].second].r = 0;
		m->colors[resemblance_pairs[i].second].g = 0;
		m->colors[resemblance_pairs[i].second].b = 255;
	}*/
	/*for (size_t i = 0; i < positive_mesh_N_lateral_descriptor.size(); i++)
	{
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[0]].r = 255;
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[0]].g = 0;
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[0]].b = 0;

		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[1]].r = 255;
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[1]].g = 0;
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[1]].b = 0;

		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[2]].r = 255;
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[2]].g = 0;
		m->colors[positive_mesh_N_lateral_descriptor[i].point_indices[2]].b = 0;
	}
	for (size_t i = 0; i < negative_mesh_N_lateral_descriptor.size(); i++)
	{
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[0]].r = 0;
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[0]].g = 0;
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[0]].b = 255;

		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[1]].r = 0;
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[1]].g = 0;
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[1]].b = 255;

		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[2]].r = 0;
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[2]].g = 0;
		m->colors[negative_mesh_N_lateral_descriptor[i].point_indices[2]].b = 255;
	}*/

	for (size_t i = 0; i < mesh_left_endpoints.size(); i++)
	{
		m->colors[mesh_left_endpoints[i]].r = 0;
		m->colors[mesh_left_endpoints[i]].g = 0;
		m->colors[mesh_left_endpoints[i]].b = 255;
	}

	for (size_t i = 0; i < mesh_right_endpoints.size(); i++)
	{
		m->colors[mesh_right_endpoints[i]].r = 0;
		m->colors[mesh_right_endpoints[i]].g = 255;
		m->colors[mesh_right_endpoints[i]].b = 0;
	}

	m->calculated_symmetry_pairs = resemblance_pairs;






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

	// file creating
	ofstream txtFile("C:/Users/BATU/source/repos/Trilateral/Results/" + m->file_name + ".txt");
	// 
	txtFile << "===============================================================================" << std::endl;
	txtFile << " N lateral " << std::endl;
	//first write the N_lateral parameters
	txtFile << " N ===== " << N_LATERAL_PARAMETERS.N << std::endl;
	txtFile << "number of sampled pairs  " << N_LATERAL_PARAMETERS.no_of_N_lateral_pairs << "  " << N_LATERAL_PARAMETERS.no_of_N_lateral_pairs * 2.0 /
		m->vertices.size() << std::endl;
	for (size_t i = 0; i < N_LATERAL_PARAMETERS.NO_OF_PARAMETERS; i++)
	{
		if (N_LATERAL_PARAMETERS.parameter_checkbox[i])
		{
			txtFile << N_LATERAL_PARAMETERS.parameter_names[i] << "  " << N_LATERAL_PARAMETERS.parameter_weights[i] << std::endl;
		}
	}
	txtFile << "===============================================================================" << std::endl;
	float error_percentage = total_error / maximum_geodesic_distance;
	txtFile << " geodesic error " + std::to_string(error_percentage) + "\n";
	txtFile.close();

}

void  NLateralDescriptor::get_ROI()
{
	float total_area = 0;
	bool is_visited_interior = false;
	//std::vector<int> roi_indices =  trilateral_ROI(&this->mesh, this->point_indices[0], this->point_indices[1], this->point_indices[2], 1, is_visited_interior);
	//trilateral_ROI_area(&this->mesh, roi_indices, total_area);

	this->area = total_area;
}


//this is trilteral
//actually select most far away  from midpoitns
std::vector<NLateralDescriptor> NLateral_select_farthest_to_midpoint(TrilateralMesh* m, std::vector<unsigned int>& indices,
int no_of_points, unsigned int smallest_agd_index , int histogram_size )
{
	std::vector<NLateralDescriptor> descs; 
	std::vector<unsigned int> nlateral_indices;
	std::vector<float> midpoint_distances = Geodesic_dijkstra(*m , smallest_agd_index);
	for (size_t i = 0; i < indices.size(); i++)
	{
		std::vector<float> distances_from_i = Geodesic_dijkstra(*m, indices[i]);
		std::vector<std::pair<float, unsigned int>> selected_indices_dist; 
		for (size_t j = 0; j < indices.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			int index = indices[j];
			selected_indices_dist.push_back(std::make_pair(distances_from_i[index], index));
		}
		//sort and pick no_of_points
		std::sort(selected_indices_dist.begin(), selected_indices_dist.end() );

		//now pick the two which are closest to midpoint
		std::vector<std::pair<float, unsigned int>>  closest_points_midpoint;
		for (size_t j = 0; j < no_of_points; j++)
		{
			closest_points_midpoint.push_back(std::make_pair(midpoint_distances[selected_indices_dist[j].second], selected_indices_dist[j].second));
		}
		std::sort(closest_points_midpoint.begin(), closest_points_midpoint.end());
		
		nlateral_indices.push_back(indices[i]);
		nlateral_indices.push_back(closest_points_midpoint[no_of_points - 1 ].second);
		nlateral_indices.push_back(closest_points_midpoint[no_of_points - 2].second);

		NLateralDescriptor desc;
		desc = NLateral_generate_descriptor(m, nlateral_indices);

		desc.n_ring_area = get_N_ring_area(m, desc.indices[0], 1);
		for (size_t i = 0; i < 3; i++)
		{
			desc.create_histogram_area(m, histogram_size, i);
		}
		for (size_t i = 0; i < 3; i++)
		{
			desc.create_histogram_HKS(m, histogram_size,i);
		}
		desc.paths_ratio = NLateral_get_paths_ratio(m, desc);
		descs.push_back(desc);


		nlateral_indices.clear();
	}
	return descs;
}




std::vector<NLateralDescriptor> NLateral_generate_closest_points(TrilateralMesh* m, std::vector<unsigned int>& indices,int N
,int histogram_size)

{
	// 1 - go gaussian to skeleton
	std::vector<unsigned int> skeleton_end_points;
	//skeleton_get_closest_skeleton_endpoints(m, skel, indices,skeleton_end_points);
	// get every skeltreeNode
	std::vector<NLateralDescriptor> nLateralDescVec;
	for (size_t i = 0; i < indices.size(); i++)
	{
		//get n-1 of the closest indexed points
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



		//get n-1 closest
		std::vector<bool> is_already_selected(m->vertices.size(), false);
		std::vector<unsigned int> selected_indices;
		std::vector<unsigned int> indices_for_nlateral_construction;
		is_already_selected[indices[i]] = true; // itself is selected 
		selected_indices.push_back(i);
		indices_for_nlateral_construction.push_back(indices[i]);
		for (size_t j = 0; j < N - 1; j++)
		{
			float maxVal = INFINITY;
			float maxIndexFirst = -1;
			for (size_t k = 0; k < distances.size(); k++)
			{
				if (maxVal > distances[k].first && !is_already_selected[distances[k].second])
				{
					maxIndexFirst = distances[k].second;
					maxVal = distances[k].first;
				}
			}
			is_already_selected[(int)maxIndexFirst] = true;
			indices_for_nlateral_construction.push_back(maxIndexFirst);
			//check which index corresponds to maxIndexFirst
			for (size_t k = 0; k <indices.size() ; k++)
			{
				if (indices[k] == maxIndexFirst)
				{
					selected_indices.push_back(k);
				}
			}
		}
		
		NLateralDescriptor desc;
		desc = NLateral_generate_descriptor(m, indices_for_nlateral_construction);
		//desc.depth = skeleton_end_point_node_list[selected_indices[0]].depth;
		
		desc.n_ring_area = get_N_ring_area(m, desc.indices[0], 1);
		desc.create_histogram_area(m, histogram_size,0);
		desc.create_histogram_HKS(m, histogram_size,0);
		//desc.create_histogram_SDF(m, histogram_size);
		desc.paths_ratio = NLateral_get_paths_ratio(m, desc);
		nLateralDescVec.push_back(desc);
	}
	return nLateralDescVec;
}

std::vector<NLateralDescriptor> NLateral_generate_with_midpoint(TrilateralMesh* m, std::vector<unsigned int>& indices,
	int N, int histogram_size)
{
	//get mid point
	glm::vec3 mid_point(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		mid_point = mid_point + m->vertices[i];
	}
	mid_point /= m->vertices.size();
	float closest = INFINITY;
	int index = -1; 
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		float dist = glm::distance(mid_point, m->vertices[i]);
		if (dist < closest)
		{
			closest = dist; 
			index = i; 
		}
	}
	std::vector<NLateralDescriptor> nLateralDescVec;
	for (size_t i = 0; i < indices.size(); i++)
	{
		for (size_t j = 0; j < indices.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			std::vector<unsigned int > nlateral_points;
			nlateral_points.push_back(indices[i]);
			nlateral_points.push_back(indices[j]);
			nlateral_points.push_back(index);
			NLateralDescriptor desc = NLateral_generate_descriptor(m, nlateral_points);
			desc.n_ring_area = get_N_ring_area(m, desc.indices[0], 1);
			desc.create_histogram_area(m, histogram_size,0);
			desc.create_histogram_HKS(m, histogram_size ,0);
			desc.paths_ratio = NLateral_get_paths_ratio(m, desc);
			nLateralDescVec.push_back(desc);
		}
	}
	return nLateralDescVec;
}
NLateralDescriptor NLateral_generate_descriptor(TrilateralMesh* m, const std::vector<unsigned int>& mesh_indices )
{
	NLateralDescriptor desc;
	int N = mesh_indices.size();

	std::vector<std::vector<std::vector<int>>> paths(N,std::vector<std::vector<int>>(N));
	desc.paths = paths;
	//generate paths
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				continue;
			}
			std::vector<int> path_i_j = Geodesic_between_two_points(*m, mesh_indices[i], mesh_indices[j]);
			desc.paths[i][j] = path_i_j;
		}
	}

	//check visited vertices
	desc.vertices_inside =  Nlateral_check_vertices_visited(m, desc);
	desc.triangles_inside =  Nlateral_check_triangles_visited(m, desc);
	desc.indices = mesh_indices; 
	desc.N = N;
	desc.area = 0;
	for (size_t i = 0; i < desc.vertices_inside.size(); i++)
	{
		desc.area += m->areas[desc.vertices_inside[i]];
	}//get max distance
	desc.max_distance = Nlateral_get_maximum_dist(m, desc);
	std::vector<float> distances = Geodesic_dijkstra(*m, desc.indices[0]);
	return desc;
}
float Nlateral_get_maximum_dist(TrilateralMesh* m, NLateralDescriptor& desc)
{
	std::vector<float> distances = Geodesic_dijkstra(*m, desc.indices[0]);
	float max = -INFINITY; 
	for (size_t i = 0; i < desc.paths.size(); i++)
	{
		for (size_t j = 0; j < desc.paths[i].size(); j++)
		{
			for (size_t k = 0; k < desc.paths[i][j].size(); k++)
			{
				int index = desc.paths[i][j][k];
				if (distances[index] > max)
				{
					max = distances[index];
				}
			}
		}
	}
	for (size_t i = 0; i < desc.vertices_inside.size(); i++)
	{
		int index = desc.vertices_inside[i];
		if (distances[index] > max)
		{
			max = distances[index];
		}
	}
	return max; 
}

std::vector<unsigned int> Nlateral_check_vertices_visited(TrilateralMesh* m, NLateralDescriptor &desc )
{

	std::vector<int> is_visited(m->vertices.size(), OUTSIDE);


	std::vector<int> edge_vertices;
	for (size_t i = 0; i < desc.paths.size(); i++)
	{
		for (size_t j = 0; j < desc.paths[i].size(); j++)
		{
			for (size_t k = 0; k < desc.paths[i][j].size(); k++)
			{
				edge_vertices.push_back(desc.paths[i][j][k]);
				is_visited[desc.paths[i][j][k]] = EDGE; 
			}
		}
	}


	std::vector<int> edge_vertices_neighbours;

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

	std::vector<unsigned int> vertices_inside_n_lateral; 
	//the minimum sized batch is our inside 
	for (size_t i = 0; i < visited_vertices_list.size(); i++)
	{
		if (visited_vertices_list[i].size() < 50.0 / 100.0 * m->vertices.size())
		{
			for (size_t j = 0; j < visited_vertices_list[i].size(); j++)
			{
				vertices_inside_n_lateral.push_back(visited_vertices_list[i][j]);
			}
		}
	}

	// Sort the vector to bring duplicates together
	std::sort(vertices_inside_n_lateral.begin(), vertices_inside_n_lateral.end());

	// Remove duplicates using std::unique, which moves unique elements to the front
	auto last = std::unique(vertices_inside_n_lateral.begin(), vertices_inside_n_lateral.end());

	// Erase the "leftover" elements after the unique elements
	vertices_inside_n_lateral.erase(last, vertices_inside_n_lateral.end());

	// if every index has same length they are colinear
	return vertices_inside_n_lateral;
}
std::vector<unsigned int> Nlateral_check_triangles_visited(TrilateralMesh* m, NLateralDescriptor& desc)
{
	std::vector<unsigned int> triangles;
	for (size_t i = 0; i < m->triangles.size(); i+=3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i+1];
		int index3 = m->triangles[i+2];
		for (size_t j = 0; j < desc.vertices_inside.size(); j++)
		{
			int index_j = desc.vertices_inside[j];
			if (index_j == index1 || index_j == index3 || index_j == index2)
			{
				triangles.push_back(index1);
				triangles.push_back(index2);
				triangles.push_back(index3);
			}
		}
	}
	return triangles; 
}
void NLateralDescriptor_write(std::string filename,TrilateralMesh* m , std::vector<NLateralDescriptor>& desc)
{
	std::ofstream file;                // Create an ofstream object for file output
	// Open the file in write mode
	file.open(filename);
	// Check if the file was opened successfully
	if (!file) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}
	
	// Write some data to the file
	//desc format 1 - 
	for (size_t i = 0; i < desc.size(); i++)
	{
		file << " NEW DESC " << std::endl;

		file << " N == " << desc[0].N << std::endl;

		file << " indices ";
		for (size_t j = 0; j < desc[i].indices.size(); j++)
		{
			file << desc[i].indices[j] << " ";
		}
		file << std::endl;
		file << " vertices inside ";
		for (size_t j = 0; j < desc[i].vertices_inside.size(); j++)
		{
			file << desc[i].vertices_inside[j] << " ";
		}
		file << std::endl;
		for (size_t j = 0; j < desc[i].paths.size(); j++)
		{
			file << " new path " << desc[i].paths.size();
			file << std::endl; 
			for (size_t k = 0; k < desc[i].paths[j].size(); k++)
			{
				file << " path of path ";
				file << std::endl;
				file << "path of points ";
				for (size_t t = 0; t < desc[i].paths[j][k].size(); t++)
				{
					file << desc[i].paths[j][k][t] << " ";
				}
				file << std::endl;
			}
		}
		file << std::endl;
		file << " descriptor end " << std::endl;
	}
	file << " resemblance pairs ";
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		std::pair<unsigned int, unsigned int> pairs;
		pairs.first = m->calculated_symmetry_pairs[i].first;
		pairs.second = m->calculated_symmetry_pairs[i].second;
		file << pairs.first << " " << pairs.second <<  " ";
	}

	// Close the file
	file.close();

}

void NLateralDescriptor_read(std::string filename, TrilateralMesh* m, std::vector<NLateralDescriptor>& desc)
{
	std::ifstream file(filename);                // Create an ofstream object for file output

	// Check if the file was opened successfully
	if (!file) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}
	// Write some data to the file
	//desc format 1 - 
	int negative_start_index = -1;
	bool is_new_desc = false;
	bool is_new_path = false;
	std::string line;
	NLateralDescriptor new_desc;
	while (std::getline(file, line))
	{
		if (line.find(" NEW DESC ") != std::string::npos)
		{
			is_new_desc = true;
			new_desc = NLateralDescriptor();
		}
		if (line.find(" N == ") != std::string::npos)
		{
			line = line.substr(5);
			std::stringstream ss(line);
			std::vector<int> nums;
			int num;
			while (ss >> num)
			{
				nums.push_back(num);
			}
			new_desc.N = nums[0];
		}
		if (line.find("indices") != std::string::npos)
		{
			line = line.substr(8);
			std::stringstream ss(line);
			std::vector<unsigned int> nums;
			int num;
			while (ss >> num)
			{
				nums.push_back(num);
			}
			new_desc.indices = nums;
		}
		if (line.find("vertices inside") != std::string::npos)
		{
			line = line.substr(16);
			std::stringstream ss(line);
			std::vector<unsigned int> visited;
			unsigned int num;
			while (ss >> num)
			{
				visited.push_back(num);
			}
			new_desc.vertices_inside = visited;
		}
		if (line.find("new path") != std::string::npos)
		{
			line = line.substr(8);
			std::stringstream ss(line);
			is_new_path = true; 
			std::vector<std::vector<int>> path;
			new_desc.paths.push_back(path);
		}
		if (line.find("path of path") != std::string::npos )
		{
			line = line.substr(12);
			std::vector<int> path;
			new_desc.paths[new_desc.paths.size()-1].push_back(path);
		}
		if (line.find("path of points") != std::string::npos)
		{
			line = line.substr(14);
			std::stringstream ss(line);
			std::vector<int> vectors;
			int num;
			while (ss >> num)
			{
				vectors.push_back(num);
			}
			std::vector<int> path;
			int index_i = new_desc.paths.size() - 1;
			int index_j = new_desc.paths[index_i].size() - 1;
			for (size_t i = 0; i < vectors.size(); i++)
			{
				new_desc.paths[index_i][index_j].push_back(vectors[i]);
			}
		}
		if (line.find(" descriptor end ") != std::string::npos)
		{
			desc.push_back(new_desc);
		}
		if (line.find("resemblance pairs") != std::string::npos)
		{
			line = line.substr(18);
			std::stringstream ss(line);
			std::vector<int> vectors;
			int num;
			while (ss >> num)
			{
				vectors.push_back(num);
			}
			for (size_t i = 0; i < vectors.size(); i+=2)
			{
				std::pair<unsigned int, unsigned int> pairs;
				pairs.first = vectors[i];
				pairs.second = vectors[i + 1];
				m->calculated_symmetry_pairs.push_back(pairs);
			}
		}
	}
	// Close the file
	file.close();
	for (size_t i = 0; i < desc.size(); i++)
	{
		for (size_t j = 0; j < desc[i].paths.size(); j++)
		{
			for (size_t k = 0; k < desc[i].paths[j].size(); k++)
			{
				for (size_t t = 0; t < desc[i].paths[j][k].size(); t++)
				{
					int index = desc[i].paths[j][k][t];
					m->raylib_mesh.colors[index * 4] = 255;
					m->raylib_mesh.colors[index * 4 + 1] = 0;
					m->raylib_mesh.colors[index * 4 + 2] = 0;
					m->raylib_mesh.colors[index * 4 + 3] = 255;
				}
			}
		}
		for (size_t j = 0; j < desc[i].vertices_inside.size(); j++)
		{
			int index = desc[i].vertices_inside[j];
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}

	m->update_raylib_mesh();
	return;
}
void Nlateral_display_desc(TrilateralMesh* m , NLateralDescriptor& desc)
{
	for (size_t i = 0; i < desc.paths.size(); i++)
	{
		for (size_t j = 0; j < desc.paths[i].size(); j++)
		{
			for (size_t k = 0; k < desc.paths[i][j].size(); k++)
			{
				int index = desc.paths[i][j][k];
				m->raylib_mesh.colors[index * 4] = 255;
				m->raylib_mesh.colors[index * 4 + 1] = 0;
				m->raylib_mesh.colors[index * 4 + 2] = 0;
				m->raylib_mesh.colors[index * 4 + 3] = 255;
			}
		}
	}

	std::vector<float> distances_from_desc = Geodesic_dijkstra(*m, desc.indices[0]);
	float dist_step = desc.max_distance / desc.area_histogram[0].size();
	for (size_t i = 0; i < desc.vertices_inside.size(); i++)
	{
		int index = desc.vertices_inside[i];
		float dist = distances_from_desc[index];
		int step = dist / dist_step;
		if (step % 2 == 0) //green; 
		{
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
		else // (step % 2 == 1) //blue
		{
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 0;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}
	for (size_t i = 0; i < desc.indices.size(); i++)
	{
		m->raylib_mesh.colors[desc.indices[i] * 4] = 0;
		m->raylib_mesh.colors[desc.indices[i] * 4 + 1] = 255;
		m->raylib_mesh.colors[desc.indices[i] * 4 + 2] = 255;
		m->raylib_mesh.colors[desc.indices[i] * 4 + 3] = 255;
	}

	m->raylib_mesh.colors[desc.indices[0] * 4] = 255;
	m->raylib_mesh.colors[desc.indices[0] * 4 + 1] = 255;
	m->raylib_mesh.colors[desc.indices[0] * 4 + 2] = 255;
	m->raylib_mesh.colors[desc.indices[0] * 4 + 3] = 255;
	m->update_raylib_mesh();
	UploadMesh(&m->raylib_mesh, false);
}
void Nlateral_display_desc(TrilateralMesh* m, std::pair<std::vector<NLateralDescriptor>, std::vector<NLateralDescriptor>>& descs, int index)
{
	NLateralDescriptor* desc = &descs.first[index];

	//make all black
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		m->raylib_mesh.colors[i * 4] = 0;
		m->raylib_mesh.colors[i * 4 + 1] = 0;
		m->raylib_mesh.colors[i * 4 + 2] = 0;
		m->raylib_mesh.colors[i * 4 + 3] = 255;
	}

	for (size_t i = 0; i < desc->paths.size(); i++)
	{
		for (size_t j = 0; j < desc->paths[i].size(); j++)
		{
			for (size_t k = 0; k < desc->paths[i][j].size(); k++)
			{
				int index = desc->paths[i][j][k];
				m->raylib_mesh.colors[index * 4] = 255;
				m->raylib_mesh.colors[index * 4 + 1] = 0;
				m->raylib_mesh.colors[index * 4 + 2] = 0;
				m->raylib_mesh.colors[index * 4 + 3] = 255;
			}
		}
	}
	/*for (size_t i = 0; i < desc->vertices_inside.size(); i++)
	{
		int index = desc->vertices_inside[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 255;
		m->raylib_mesh.colors[index * 4 + 2] = 0;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}*/
	std::vector<float> distances_from_desc  = Geodesic_dijkstra(*m, desc->indices[0]);
	float dist_step = desc->max_distance / desc->area_histogram[0].size();
	for (size_t i = 0; i < desc->vertices_inside.size(); i++)
	{
		int index = desc->vertices_inside[i];
		float dist = distances_from_desc[index];
		int step = dist / dist_step;
		if (step % 2 == 0) //green; 
		{
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
		else // (step % 2 == 1) //blue
		{
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 0;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}
	for (size_t i = 0; i < desc->indices.size(); i++)
	{
		m->raylib_mesh.colors[desc->indices[i] * 4] = 0;
		m->raylib_mesh.colors[desc->indices[i] * 4 + 1] = 255;
		m->raylib_mesh.colors[desc->indices[i] * 4 + 2] = 255;
		m->raylib_mesh.colors[desc->indices[i] * 4 + 3] = 255;
	}

	m->raylib_mesh.colors[desc->indices[0] * 4] = 255;
	m->raylib_mesh.colors[desc->indices[0] * 4 + 1] = 255;
	m->raylib_mesh.colors[desc->indices[0] * 4 + 2] = 255;
	m->raylib_mesh.colors[desc->indices[0] * 4 + 3] = 255;




	NLateralDescriptor* desc_neg = &descs.second[index];

	for (size_t i = 0; i < desc_neg->paths.size(); i++)
	{
		for (size_t j = 0; j < desc_neg->paths[i].size(); j++)
		{
			for (size_t k = 0; k < desc_neg->paths[i][j].size(); k++)
			{
				int index = desc_neg->paths[i][j][k];
				m->raylib_mesh.colors[index * 4] = 255;
				m->raylib_mesh.colors[index * 4 + 1] = 0;
				m->raylib_mesh.colors[index * 4 + 2] = 0;
				m->raylib_mesh.colors[index * 4 + 3] = 255;
			}
		}
	}
	for (size_t i = 0; i < desc_neg->vertices_inside.size(); i++)
	{
		int index = desc_neg->vertices_inside[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 0;
		m->raylib_mesh.colors[index * 4 + 2] = 255;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}

	for (size_t i = 0; i < desc_neg->indices.size(); i++)
	{
		m->raylib_mesh.colors[desc_neg->indices[i] * 4] = 0;
		m->raylib_mesh.colors[desc_neg->indices[i] * 4 + 1] = 255;
		m->raylib_mesh.colors[desc_neg->indices[i] * 4 + 2] = 255;
		m->raylib_mesh.colors[desc_neg->indices[i] * 4 + 3] = 255;
	}

	m->raylib_mesh.colors[desc_neg->indices[0] * 4] = 255;
	m->raylib_mesh.colors[desc_neg->indices[0] * 4 + 1] = 255;
	m->raylib_mesh.colors[desc_neg->indices[0] * 4 + 2] = 255;
	m->raylib_mesh.colors[desc_neg->indices[0] * 4 + 3] = 255;

	m->update_raylib_mesh();
}
void Nlateral_display_desc(TrilateralMesh* m, std::vector<NLateralDescriptor>& descs, int index)
{
	NLateralDescriptor* desc = &descs[index];

	m->color_all(LIGHTGRAY);
	for (size_t i = 0; i < desc->paths.size(); i++)
	{
		for (size_t j = 0; j < desc->paths[i].size(); j++)
		{
			for (size_t k = 0; k < desc->paths[i][j].size(); k++)
			{
				int index = desc->paths[i][j][k];
				m->raylib_mesh.colors[index * 4] = 255;
				m->raylib_mesh.colors[index * 4 + 1] = 0;
				m->raylib_mesh.colors[index * 4 + 2] = 0;
				m->raylib_mesh.colors[index * 4 + 3]= 255;
			}
		}
	}
	for (size_t i = 0; i < desc->vertices_inside.size(); i++)
	{
		int index = desc->vertices_inside[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 255;
		m->raylib_mesh.colors[index * 4 + 2] = 0;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}

	for (size_t i = 0; i < desc->indices.size(); i++)
	{
		m->raylib_mesh.colors[desc->indices[i] * 4] = 0;
		m->raylib_mesh.colors[desc->indices[i] * 4 + 1] = 255;
		m->raylib_mesh.colors[desc->indices[i] * 4 + 2] = 255;
		m->raylib_mesh.colors[desc->indices[i] * 4 + 3] = 255;
	}

	m->raylib_mesh.colors[desc->indices[0] * 4] = 255;
	m->raylib_mesh.colors[desc->indices[0] * 4 + 1] = 255;
	m->raylib_mesh.colors[desc->indices[0] * 4 + 2] = 255;
	m->raylib_mesh.colors[desc->indices[0] * 4 + 3] = 255;

	//check if matched
	if (m->calculated_symmetry_pairs.size() > 0)
	{
		int new_index = -1; 
		//get match
		for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
		{
			int index_i1 = m->calculated_symmetry_pairs[i].first;
			int index_i2 = m->calculated_symmetry_pairs[i].second;
			if (index_i1 == desc->indices[0])
			{
				new_index = index_i2;
				break; 
			}
			else if (index_i2 == desc->indices[0])
			{
				new_index = index_i1;
				break;
			}
		}
		if (new_index == -1)
		{
			m->update_raylib_mesh();
		}

		if (new_index != -1)
		{
			//get the corresponding index
			for (size_t d = 0; d < descs.size(); d++)
			{
				if (descs[d].indices[0] == new_index)
				{
					for (size_t i = 0; i < descs[d].paths.size(); i++)
					{
						for (size_t j = 0; j < descs[d].paths[i].size(); j++)
						{
							for (size_t k = 0; k < descs[d].paths[i][j].size(); k++)
							{
								int index = descs[d].paths[i][j][k];
								m->raylib_mesh.colors[index * 4] = 255;
								m->raylib_mesh.colors[index * 4 + 1] = 0;
								m->raylib_mesh.colors[index * 4 + 2] = 0;
								m->raylib_mesh.colors[index * 4 + 3] = 255;
							}
						}
					}
					for (size_t i = 0; i < descs[d].vertices_inside.size(); i++)
					{
						int index = descs[d].vertices_inside[i];
						m->raylib_mesh.colors[index * 4] = 0;
						m->raylib_mesh.colors[index * 4 + 1] = 0;
						m->raylib_mesh.colors[index * 4 + 2] = 255;
						m->raylib_mesh.colors[index * 4 + 3] = 255;
					}

					for (size_t i = 0; i < descs[d].indices.size(); i++)
					{
						m->raylib_mesh.colors[descs[d].indices[i] * 4] = 0;
						m->raylib_mesh.colors[descs[d].indices[i] * 4 + 1] = 255;
						m->raylib_mesh.colors[descs[d].indices[i] * 4 + 2] = 255;
						m->raylib_mesh.colors[descs[d].indices[i] * 4 + 3] = 255;
					}

					m->raylib_mesh.colors[descs[d].indices[0] * 4] = 255;
					m->raylib_mesh.colors[descs[d].indices[0] * 4 + 1] = 255;
					m->raylib_mesh.colors[descs[d].indices[0] * 4 + 2] = 255;
					m->raylib_mesh.colors[descs[d].indices[0] * 4 + 3] = 255;

					std::cout << "match ==  i " << index << " j " << d << std::endl;
				}
			}
		}
	}
	std::vector<float> distances_from_desc = Geodesic_dijkstra(*m, desc->indices[0]);
	float dist_step = desc->max_distance / desc->area_histogram[0].size();
	for (size_t i = 0; i < desc->vertices_inside.size(); i++)
	{
		int index = desc->vertices_inside[i];
		float dist = distances_from_desc[index];
		int step = dist / dist_step;
		if (step % 2 == 0) //green; 
		{
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 255;
			m->raylib_mesh.colors[index * 4 + 2] = 0;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
		else // (step % 2 == 1) //blue
		{
			m->raylib_mesh.colors[index * 4] = 0;
			m->raylib_mesh.colors[index * 4 + 1] = 0;
			m->raylib_mesh.colors[index * 4 + 2] = 255;
			m->raylib_mesh.colors[index * 4 + 3] = 255;
		}
	}
	m->update_raylib_mesh();
	UploadMesh(&m->raylib_mesh, false);

}

void NLateralDescriptor::create_histogram_HKS(TrilateralMesh* m, int hist_no, int descriptor_point_index)
{
	Histogram hist;
	hist.init(hist_no);
	std::vector<float> distances = Geodesic_dijkstra(*m, this->indices[descriptor_point_index]);
	float longest_path = -INFINITY;
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		if (longest_path < distances[this->vertices_inside[i]])
		{
			longest_path = distances[this->vertices_inside[i]];
		}
	}
	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				if (longest_path < distances[this->paths[i][j][k]])
				{
					longest_path = distances[this->paths[i][j][k]];
				}
			}
		}
	}
	float step_size = longest_path / hist_no;
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		int index = this->vertices_inside[i];
		int step_1 = distances[index] / step_size; 
		if (step_1 == hist_no)
		{
			step_1--; 
		}
		hist.histogram[step_1] += m->normalized_heat_kernel_signature[index];
	}
	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				int step_i = distances[this->paths[i][j][k]] / step_size;

				if (step_i == hist_no)
				{
					step_i--;
				}

				hist.histogram[step_i] += m->normalized_heat_kernel_signature[this->paths[i][j][k]];
			}
		}
	}
	this->hks_histogram[descriptor_point_index] = hist;
}
void NLateralDescriptor::create_histogram_area_HKS_combined(TrilateralMesh* m, int hist_no)
{
	Histogram hist;
	hist.init(hist_no);
	//get tlongest path
	std::vector<float> distances = Geodesic_dijkstra(*m, this->indices[0]);
	float longest_path = -INFINITY;
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		if (longest_path < distances[this->vertices_inside[i]])
		{
			longest_path = distances[this->vertices_inside[i]];
		}
	}
	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				if (longest_path < distances[this->paths[i][j][k]])
				{
					longest_path = distances[this->paths[i][j][k]];
				}
			}
		}
	}
	float step_size = longest_path / hist_no;
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		int index = this->vertices_inside[i];
		int step_1 = distances[index] / step_size;
		if (step_1 == hist_no)
		{
			step_1--;
		}
		hist.histogram[step_1] += m->normalized_heat_kernel_signature[index] * m->areas[index];
	}

	//lets ignore the paths only
	this->area_hks_histogram = hist; 
}
void NLateralDescriptor::create_histogram_SDF(TrilateralMesh* m, int hist_no)
{
	Histogram hist;
	hist.init(hist_no);
	std::vector<float> distances = Geodesic_dijkstra(*m, this->indices[0]);
	float longest_path = -INFINITY;
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		if (longest_path < distances[this->vertices_inside[i]])
		{
			longest_path = distances[this->vertices_inside[i]];
		}
	}

	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				if (longest_path < distances[this->paths[i][j][k]])
				{
					longest_path = distances[this->paths[i][j][k]];
				}
			}
		}
	}
	float step_size = longest_path / hist_no;
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		int index = this->vertices_inside[i];
		int step_1 = distances[index] / step_size;
		if (step_1 == hist_no)
		{
			step_1--;
		}
		hist.histogram[step_size] += m->sdf[index];
	}
	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				int step_i = distances[this->paths[i][j][k]] / step_size;

				if (step_i == hist_no)
				{
					step_i--;
				}

				hist.histogram[step_i] += m->sdf[this->paths[i][j][k]];
			}
		}
	}
	this->sdf_histogram = hist;
}
void NLateralDescriptor::create_histogram_area(TrilateralMesh* m, int hist_no ,int descriptor_point_index)
{
	Histogram hist;
	hist.init(hist_no);
	std::vector<float> distances = Geodesic_dijkstra(*m, this->indices[descriptor_point_index]);
	float longest_path = -INFINITY; 
	for (size_t i = 0; i < this->vertices_inside.size(); i++)
	{
		if (longest_path < distances[this->vertices_inside[i]])
		{
			longest_path = distances[this->vertices_inside[i]];
		}
	}
	
	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				if (longest_path < distances[this->paths[i][j][k]])
				{
					longest_path = distances[this->paths[i][j][k]];
				}
			}
		}
	}

	float step_size = longest_path / hist_no;
	for (size_t i = 0; i < this->triangles_inside.size(); i+=3)
	{
		int index_1 = this->triangles_inside[i];
		int index_2 = this->triangles_inside[i + 1];
		int index_3 = this->triangles_inside[i + 2];
		int step_1 = distances[index_1] / step_size;
		int step_2 = distances[index_2] / step_size;
		int step_3 = distances[index_3] / step_size;
		float triangle_area = compute_triangle_area(m->vertices[index_1], m->vertices[index_2], m->vertices[index_3]);
		if (step_1 == hist_no)
		{
			step_1--;
		}
		if (step_2 == hist_no)
		{
			step_2--;
		}
		if (step_3 == hist_no)
		{
			step_3--;
		}
		// compute the areas per triangle
		int min = std::min( step_1, std::min(step_2, step_3));
		int max = std::max( step_1, std::max(step_2, step_3));
		if (min == max)
		{
			hist.histogram[min] += triangle_area;
		}
		for (size_t j = min; j < max; j++)
		{
			hist.histogram[j] += ( triangle_area / (max-min) );
		}
	}

	for (size_t i = 0; i < this->paths.size(); i++)
	{
		for (size_t j = 0; j < this->paths[i].size(); j++)
		{
			for (size_t k = 0; k < this->paths[i][j].size(); k++)
			{
				int step_i = distances[this->paths[i][j][k]] / step_size;

				if (step_i == hist_no)
				{
					step_i--;
				}

				hist.histogram[step_i] += m->areas[this->paths[i][j][k]];
			}
		}
	}
	this->area_histogram[descriptor_point_index] = hist;
}


void NLateral_compute_skel_point_dist(TrilateralMesh* m,Skeleton& skel, NLateralDescriptor& desc)
{
	int index = desc.indices[0];
	int skel_index = desc.skeleton_index;
	desc.skel_point_dist = glm::distance(skel.skeletonFormat[skel_index].point , m->vertices[index]);
}


std::vector<unsigned int> NLateral_sweepDistance(TrilateralMesh* m, std::vector<unsigned int> indices, float sweep_distance)
{
	std::vector<unsigned int> extracted_indices;
	while (indices.size() > 0)
	{
		unsigned int p = indices[0];
		extracted_indices.push_back(p);
		indices.erase(indices.begin(), indices.begin() + 1);
		int index = p;
		std::vector<float> distances = Geodesic_dijkstra(*m, index);
		std::vector<unsigned int> close_vertices;
		for (size_t i = 0; i < indices.size(); i++)
		{
			if (distances[indices[i]] < sweep_distance)
			{
				close_vertices.push_back(i);
			}
		}

		std::sort(close_vertices.begin(), close_vertices.end(), std::greater<unsigned int>());

		// Remove elements from the vector starting from the highest index
		for (int index : close_vertices) {
			if (index >= 0 && index < indices.size()) {
				indices.erase(indices.begin() + index);
			}
			else {
				std::cerr << "Invalid index: " << index << std::endl;
			}
		}

	}
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		m->raylib_mesh.colors[i * 4] = 0;
		m->raylib_mesh.colors[i * 4 + 1] = 0;
		m->raylib_mesh.colors[i * 4 + 2] = 0;
		m->raylib_mesh.colors[i * 4 + 3] = 255;
	}
	for (size_t i = 0; i < extracted_indices.size(); i++)
	{
		int index = extracted_indices[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 0;
		m->raylib_mesh.colors[index * 4 + 2] = 255;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}
	m->update_raylib_mesh();

	return extracted_indices;
}

//for now only for N = 3 
bool NLateral_check_path_lengths(TrilateralMesh* m , NLateralDescriptor& desc1, NLateralDescriptor& desc2, float similarity)
{
	float desc1_1_2_length = 0;
	for (size_t i = 0; i < desc1.paths[0][1].size()-1; i++)
	{
		int index = desc1.paths[0][1][i];
		int index_1 = desc1.paths[0][1][i+1];
		float length = glm::distance(m->vertices[index], m->vertices[index_1]);
		desc1_1_2_length += length; 
	}
	float desc1_1_3_length = 0;
	for (size_t i = 0; i < desc1.paths[0][2].size() - 1; i++)
	{
		int index = desc1.paths[0][2][i];
		int index_1 = desc1.paths[0][2][i + 1];
		float length = glm::distance(m->vertices[index], m->vertices[index_1]);
		desc1_1_3_length += length;
	}

	float desc1_ratio = desc1_1_2_length / desc1_1_3_length;
	if (desc1_ratio > 1)
	{
		desc1_ratio = 1.0f / desc1_ratio;
	}

	float desc2_1_2_length = 0;
	for (size_t i = 0; i < desc2.paths[0][1].size() - 1; i++)
	{
		int index = desc2.paths[0][1][i];
		int index_1 = desc2.paths[0][1][i + 1];
		float length = glm::distance(m->vertices[index], m->vertices[index_1]);
		desc2_1_2_length += length;
	}
	float desc2_1_3_length = 0;
	for (size_t i = 0; i < desc2.paths[0][2].size() - 1; i++)
	{
		int index = desc2.paths[0][2][i];
		int index_1 = desc2.paths[0][2][i + 1];
		float length = glm::distance(m->vertices[index], m->vertices[index_1]);
		desc2_1_3_length += length;
	}

	float desc2_ratio = desc2_1_2_length / desc2_1_3_length;
	if (desc2_ratio > 1)
	{
		desc2_ratio = 1.0f / desc2_ratio;
	}

	if (std::fabsf(desc1_ratio - desc2_ratio) > similarity)
	{
		return false;
	}
	return true; 


}
bool Nlateral_compare_angles(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float angle_param)
{
	//first 
	glm::vec3 vec1 = m->vertices[desc1.paths[0][1][desc1.paths[0][1].size() - 2]] - m->vertices[desc1.indices[0]];
	glm::vec3 vec2 = m->vertices[desc1.paths[0][2][desc1.paths[0][2].size() - 2]] - m->vertices[desc1.indices[0]];
	float angle = glm::dot(vec1, vec2) / (glm::length(vec1) * glm::length(vec2));
	float cosTheta = glm::clamp(angle, -1.0f, 1.0f);
	float angleRadians = std::acos(cosTheta);
	float angleDegrees_1 = glm::degrees(angleRadians);

	vec1 = m->vertices[desc2.paths[0][1][desc2.paths[0][1].size() - 2]] - m->vertices[desc2.indices[0]];
	vec2 = m->vertices[desc2.paths[0][2][desc2.paths[0][2].size() - 2]] - m->vertices[desc2.indices[0]];
	angle = glm::dot(vec1, vec2) / (glm::length(vec1) * glm::length(vec2));
	cosTheta = glm::clamp(angle, -1.0f, 1.0f);
	angleRadians = std::acos(cosTheta);
	float angleDegrees_2 = glm::degrees(angleRadians);
	if (std::abs(angleDegrees_1 - angleDegrees_2) > angle_param)
	{
		return false; 
	}
	return true; 
}
bool NLateral_compare_HKS(TrilateralMesh* m,NLateralDescriptor& desc1, NLateralDescriptor& desc2 , float hks_perc, std::ofstream& file)
{

	file << " HKS " << std::endl;
	file << " hks comp indices 0 0 " << std::abs(m->normalized_heat_kernel_signature[desc1.indices[0]] - m->normalized_heat_kernel_signature[desc2.indices[0]])<< std::endl; 
	// check indices[0]
	if (std::abs(m->normalized_heat_kernel_signature[desc1.indices[0]] - m->normalized_heat_kernel_signature[desc2.indices[0]]) > hks_perc)
	{
		return false; 
	}
	//check others
	int index1_2 = desc1.indices[1];
	int index1_3 = desc1.indices[2];
	int index2_2 = desc2.indices[1];
	int index2_3 = desc2.indices[2];

	file << " hks comp indices 1 1 " << std::abs(m->normalized_heat_kernel_signature[index1_2] - m->normalized_heat_kernel_signature[index2_2]) << std::endl;
	file << " hks comp indices 1 2 " << std::abs(m->normalized_heat_kernel_signature[index1_2] - m->normalized_heat_kernel_signature[index2_3]) << std::endl;

	if (std::abs(m->normalized_heat_kernel_signature[index1_2] - m->normalized_heat_kernel_signature[index2_2]) > hks_perc &&
		std::abs(m->normalized_heat_kernel_signature[index1_2] - m->normalized_heat_kernel_signature[index2_3]) > hks_perc)
	{
		return false; 
	}
	file << " hks comp indices 2 1 " << std::abs(m->normalized_heat_kernel_signature[index1_3] - m->normalized_heat_kernel_signature[index2_2]) << std::endl;
	file << " hks comp indices 2 2 " << std::abs(m->normalized_heat_kernel_signature[index1_3] - m->normalized_heat_kernel_signature[index2_3]) << std::endl;

	if (std::abs(m->normalized_heat_kernel_signature[index1_3] - m->normalized_heat_kernel_signature[index2_2]) > hks_perc &&
		std::abs(m->normalized_heat_kernel_signature[index1_3] - m->normalized_heat_kernel_signature[index2_3]) > hks_perc)
	{
		return false;
	}
	return true;
}

bool NLateral_compare_skeldist_mid(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float skel_percentage, float maximum_skel,
	std::ofstream& file )
{
	file << " skel dist comp " << std::endl;
	float skel_dist = std::abs(desc1.skel_dist_mid[0] - desc2.skel_dist_mid[0]);
	file << " skel dist 0 0 " << skel_dist / maximum_skel;
	if (skel_dist / maximum_skel >  skel_percentage)
	{
		return false;
	}

	float skel_dist_2_2 = std::abs(desc1.skel_dist_mid[1] - desc2.skel_dist_mid[1]);
	float skel_dist_2_3 = std::abs(desc1.skel_dist_mid[1] - desc2.skel_dist_mid[2]);
	float skel_dist_3_3 = std::abs(desc1.skel_dist_mid[2] - desc2.skel_dist_mid[2]);

	file << " skel dist 1  1 " << skel_dist_2_2 / maximum_skel  <<  std::endl;
	file << " skel dist 1  2 " << skel_dist_2_2 / maximum_skel  <<  std::endl;
	file << " skel dist 2  2 " << skel_dist_2_2 / maximum_skel  <<  std::endl;

	if (skel_dist_2_2 / maximum_skel > skel_percentage && skel_dist_2_3 / maximum_skel > skel_percentage)
	{
		return false; 
	}
	
	if (skel_dist_3_3 / maximum_skel > skel_percentage && skel_dist_2_3 / maximum_skel > skel_percentage)
	{
		return false; 
	}

	return true; 
}

bool NLateral_compare_Nring(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float maximum_n_ring,std::ofstream& file)
{

	return true; 
}

bool NLateral_compare_area(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float maximum_area , float area_percentage
,std::ofstream& file )
{

	float area_1_2 = m->areas[desc1.indices[1]];
	float area_1_3 = m->areas[desc1.indices[2]];
	float area_2_2 = m->areas[desc2.indices[1]];
	float area_2_3 = m->areas[desc2.indices[2]];

	file << " 1 2 " << std::abs(area_1_2 - area_2_2) / maximum_area << std::endl;
	file << " 1 3 " << std::abs(area_1_2 - area_2_3) / maximum_area << std::endl;
	file << " 2 3 " << std::abs(area_1_3 - area_2_3) / maximum_area << std::endl;

	if (std::abs(area_1_2 - area_2_2) / maximum_area > area_percentage &&
		std::abs(area_1_2 - area_2_3) / maximum_area > area_percentage)
	{
		return false; 
	}
	if (std::abs(area_1_3 - area_2_2) / maximum_area > area_percentage &&
		std::abs(area_1_3 - area_2_3) / maximum_area > area_percentage)
	{
		return false;
	}

	return true; 

}

double NLateral_get_paths_ratio(TrilateralMesh* m, NLateralDescriptor& desc)
{
	float dist1 = 0;
	for (size_t i = 0; i < desc.paths[0][1].size()-1; i++)
	{
		int index1 = desc.paths[0][1][i];
		int index2 = desc.paths[0][1][i+1];
		dist1 += glm::distance(m->vertices[index1], m->vertices[index2]);
	}
	float dist2 = 0;
	for (size_t i = 0; i < desc.paths[0][2].size() - 1; i++)
	{
		int index1 = desc.paths[0][2][i];
		int index2 = desc.paths[0][2][i + 1];
		dist2 += glm::distance(m->vertices[index1], m->vertices[index2]);
	}

	if (dist1 / dist2 <= 1)
	{
		return dist1 / dist2;
	}
	return dist2 / dist1; 
}
bool NLatera_compare_depth(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, int depth_dif_param, std::ofstream& file)
{

	// 1 - display depths
	for (size_t i = 0; i < desc1.depth.size(); i++)
	{
		file << " i " << i << " " << desc1.depth[i] << std::endl;
	}
	for (size_t i = 0; i < desc2.depth.size(); i++)
	{
		file << " j " << i << " " << desc2.depth[i] << std::endl;
	}
	std::vector<int> is_dif_ok(desc1.indices.size(), 0);
	for (size_t i = 0; i < desc1.depth.size(); i++)
	{
		int depth_i = desc1.depth[i];
		for (size_t j = 0; j < desc2.depth.size(); j++)
		{
			int depth_j = desc2.depth[j];
			if (std::abs(depth_i - depth_j) < depth_dif_param)
			{
				is_dif_ok[i] = 1; 
			}
		}
	}

	for (size_t i = 0; i < is_dif_ok.size(); i++)
	{
		if (is_dif_ok[i] != 1)
		{
			return false;
		}
	}
	return true; 
}

unsigned int NLateral_get_closest_index_to_midpoint(TrilateralMesh* m, std::vector<unsigned int>& points )
{
	//midpoint 
	glm::vec3 mid_point(0, 0, 0);
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		mid_point = mid_point + m->vertices[i];
	}
	mid_point /=  m->vertices.size();

	float closest_dist = INFINITY;
	int closest_index = -1;
	for (size_t i = 0; i < points.size(); i++)
	{
		float dist = glm::distance(m->vertices[points[i]], mid_point);
		if (dist < closest_dist)
		{
			dist = closest_dist; 
			closest_index = points[i];
		}
	}
	return closest_index;
}
bool NLateral_compare_trilateral_with_midpoint(TrilateralMesh* m, unsigned int p1, unsigned int p2, unsigned int p_middle, float dissimilarity
,std::ofstream& file )
{
	NLateralDescriptor desc1;
	NLateralDescriptor desc2;
	std::vector<unsigned int> indices_p1_desc = { p1,p2,p_middle };
	std::vector<unsigned int> indices_p2_desc = { p2,p1,p_middle };
	desc1 = NLateral_generate_descriptor(m, indices_p1_desc);
	desc2 = NLateral_generate_descriptor(m, indices_p2_desc);
	desc1.create_histogram_area(m, 10,0);
	desc2.create_histogram_area(m, 10,0);


	desc1.area_histogram[0].normalize(1);
	desc2.area_histogram[0].normalize(1);
	float dif = Histogram_ChiSquareDistance(desc1.area_histogram[0], desc2.area_histogram[0]);

	file << " trilateral dissimilarity " << dif <<  " " << (dif < dissimilarity) << std::endl;

	if (dif > dissimilarity)
	{
		return false;
	}
	return true;
}

bool NLateral_compare_distance_to_midpoint(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2 , unsigned int midpoint_index
,float distance_to_mid_param,std::ofstream& file )
{
	std::vector<float> distances = Geodesic_dijkstra(*m,midpoint_index);
	float dist1 = distances[desc1.indices[0]];
	float dist2 = distances[desc2.indices[0]];
	
	if (dist1 == 0 || dist2 == 0)
	{
		return true; 
	}
	float ratio = dist1 / dist2;
	if (ratio > 1)
	{
		ratio = 1 / ratio; 
	}
	file << " distance to modpoint " << ratio << std::endl;
	
	if (ratio < distance_to_mid_param)
	{
		return false;
	}
	return true; 
}

bool Nlateral_check_endpoint(TrilateralMesh* m, Skeleton& skel, NLateralDescriptor& desc1, NLateralDescriptor& desc2)
{
	std::vector<unsigned int> mesh_endpoints;
	skeleton_calculate_closest_mesh_points(skel, m, mesh_endpoints);
	bool is_endpoints[] = { false , false };
	int index1 = desc1.indices[0];
	int index2 = desc2.indices[1];
	for (size_t i = 0; i < mesh_endpoints.size(); i++)
	{
		if (index1 == mesh_endpoints[i])
		{
			is_endpoints[0] = true;
		}
		if (index1 == mesh_endpoints[i])
		{
			is_endpoints[1] = true; 
		}
	}
	if (is_endpoints[0] == is_endpoints[1])
	{
		return true; 
	}
	return false; 
}

bool NLateral_compare_SDF(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float maximum_sdf, 
float param , std::ofstream& file)
{
	float sdf_1 = ShapeDiameter_calculate_simple(m, desc1.indices[0]);
	float sdf_2 = ShapeDiameter_calculate_simple(m, desc2.indices[0]);
	
	float dif = std::abs(sdf_1 - sdf_2);
	file << " Sdf " << dif / maximum_sdf << std::endl;
	if ( (dif /maximum_sdf)  < param )
	{
		return true; 
	}
	return false; 
}

bool NLateral_compare_position_to_midpoint(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, unsigned int midpoint_index,
	float distances_from_mid, float distances_between_desc,  std::ofstream& file)
{
	int index1 = desc1.indices[0];
	int index2 = desc2.indices[0];

	std::vector<float> distances = Geodesic_dijkstra(*m, midpoint_index);
	std::vector<float> distances_form_desc1 = Geodesic_dijkstra(*m, index1);

	// get maximum distance
	float max_dist_from_mid = -INFINITY; 
	for (size_t i = 0; i < distances.size(); i++)
	{
		if (max_dist_from_mid < distances[i])
		{
			max_dist_from_mid = distances[i];
		}
	}
	float max_dist_from_index1 = -INFINITY;
	for (size_t i = 0; i < distances.size(); i++)
	{
		if (max_dist_from_index1 < distances[i])
		{
			max_dist_from_index1 = distances[i];
		}
	}
	file << " position to mid point " << std::endl;
	file << "dist from desc1 " << distances_form_desc1[index2] << std::endl;
	file << "max dist from desc1 " << max_dist_from_index1 << std::endl;
	file << "distances from  index1 " << distances[index1] << std::endl;
	file << "maximum distance from mid  " << max_dist_from_mid << std::endl;
	//check if the descriptors are close 
	if (distances_form_desc1[index2] / max_dist_from_index1 <  distances_between_desc &&
		distances[index1] / max_dist_from_mid  > distances_from_mid)
	{
		return false; 
	}
	return true;

}

void Nlateral_write_matching_points(TrilateralMesh* m)
{
	std::ofstream file("matching_points.txt");
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int first = m->calculated_symmetry_pairs[i].first;
		int second = m->calculated_symmetry_pairs[i].second;
		file << " " << first << " " << second << "\n";
	}
	file.close();

}
void Nlateral_read_matching_points(TrilateralMesh* m)
{
	std::ifstream file("matching_points.txt");
	// Read the file line by line
	std::vector<std::pair<unsigned int, unsigned int>> sym_pair; 
	std::string line; 
	while (std::getline(file, line)) {
		std::stringstream ss(line); // Use stringstream to parse the line
		int number;
		std::vector<int> nums;
		// Extract numbers from the line
		while (ss >> number) {
			nums.push_back(number);
		}
		sym_pair.push_back(std::make_pair(nums[0] ,nums[1]));
	}
	m->calculated_symmetry_pairs = sym_pair; 
}

bool NLateral_compare_HKS(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float hks_param)
{
	float dif_0_0 = std::abs(m->normalized_heat_kernel_signature[desc1.indices[0]] - m->normalized_heat_kernel_signature[desc2.indices[0]]);
	if (dif_0_0 > hks_param)
	{
		return false; 
	}
	float dif_1_1 = std::abs(m->normalized_heat_kernel_signature[desc1.indices[1]] - m->normalized_heat_kernel_signature[desc2.indices[1]]);
	float dif_1_2 = std::abs(m->normalized_heat_kernel_signature[desc1.indices[1]] - m->normalized_heat_kernel_signature[desc2.indices[2]]);
	if (dif_1_1 > hks_param )
	{
		return false;
	}
	float dif_2_1 = std::abs(m->normalized_heat_kernel_signature[desc1.indices[2]] - m->normalized_heat_kernel_signature[desc2.indices[1]]);
	float dif_2_2 = std::abs(m->normalized_heat_kernel_signature[desc1.indices[2]] - m->normalized_heat_kernel_signature[desc2.indices[2]]);
	if (dif_2_2 > hks_param )
	{
		return false;
	}
	return true; 
}

bool NLateral_compare_FuzzyGeodesics(TrilateralMesh* m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float fuzzy_param)
{

	if (std::abs(desc1.fuzzy - desc2.fuzzy) < fuzzy_param )
	{
		return true; 
	}
	return false; 
}


void NLateral_calculate_fuzzyGeodesics(TrilateralMesh* m, std::vector<NLateralDescriptor>& descs)
{
	//compute every fuzzygeodesic
	for (size_t i = 0; i < descs.size(); i++)
	{
		float total_fuzziness_d1 = 0;
		for (size_t j = 0; j < descs[i].vertices_inside.size(); j++)
		{
			int index = descs[i].vertices_inside[j];
			float fuzziness = FuzzyGeodesic_fuzzyGeo(m, index, descs[i].indices[0], descs[i].indices[1]);
			fuzziness += FuzzyGeodesic_fuzzyGeo(m, index, descs[i].indices[0], descs[i].indices[2]);
			fuzziness += FuzzyGeodesic_fuzzyGeo(m, index, descs[i].indices[1], descs[i].indices[2]);

			total_fuzziness_d1 = fuzziness / 3;
		}
		total_fuzziness_d1 = total_fuzziness_d1 / descs[i].vertices_inside.size();
		descs[i].fuzzy = total_fuzziness_d1;
	}
	
}