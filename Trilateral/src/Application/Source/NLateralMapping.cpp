#include "../Include/NLateralMapping.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/Sampling.h"
#include "../Include/MetricCalculations.h"
#include "../Include/VarianceMinimizingTransportPlan.h"
#include "../Include/Geodesic.h"
#include "../Include/ShapeDiameter.h"
#include "../Include/CurvatureGeneration.h"

std::vector<NLateralDescriptor> NlateralMap_point_matching_w_average_geodesic(TrilateralMesh* m, 
	int dvorak_enpoint_no, float sweep_distance, float hks_dif_param,  float closeness_param, 
	float area_dif_param,float min_geo_tau,int avg_n_ring,
	float distance_to_mid_param , float sdf_param , int N , std::vector<unsigned int>& agd_point_indices)
{
	int size = m->vertices.size();
	int mesh_mid_point_index = -1;
	std::vector<NLateralDescriptor> descs;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<std::pair<int, float>> hks_pairs;
	std::vector<DvorakPairs> dvorak_pairs;
	std::vector<unsigned int> points_pos;
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	Metric_set_gaussian(m, dvorak_enpoint_no, sweep_distance);
	Metric_set_N(N);

	//for sdf calculation
	//m->calculate_sdf();
	if (agd_point_indices.size() == 0)
	{
		float biggest_num = 0;
		agd_point_indices = Geodesic_avg_dijkstra_modified(m, sweep_distance, avg_n_ring, false, biggest_num);
		for (size_t i = 0; i < 3; i++)
		{
			agd_point_indices = Geodesic_min_dijkstra(m, agd_point_indices, sweep_distance, min_geo_tau, false);
		}
	}
	
	
	//point_indices = NLateral_sweepDistance(m, point_indices, sweep_distance);
	for (size_t i = 0; i < agd_point_indices.size(); i++)
	{
		int index = agd_point_indices[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 255;
		m->raylib_mesh.colors[index * 4 + 2] = 0;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}
	int hist_size = 5;
	//descs = NLateral_generate_closest_points(m, agd_point_indices,  N, hist_size);
	unsigned int mid_point_index = -1;
	unsigned int mid_point_index_2 = -1;
	Geodesic_mid_point_w_AGD(m, mid_point_index, mid_point_index_2);
	descs = NLateral_select_farthest_to_midpoint(m, agd_point_indices,5, mid_point_index, 10);

	for (size_t i = 0; i < agd_point_indices.size(); i++)
	{
		int index = agd_point_indices[i];
		DvorakPairs p;
		p.p_index = index;
		p.gaussian_curv = gaussian_curvature(m, index);
		dvorak_pairs.push_back(p);
	}


	std::vector<std::vector<float>> hist_diffs;
	for (size_t i = 0; i < descs.size(); i++)
	{
		std::vector<float> hist_i;
		for (size_t j = 0; j < 3; j++)
		{
			descs[i].area_histogram[j].normalize(1);
			descs[i].hks_histogram[j].normalize(1);
		}
		for (size_t j = 0; j < descs.size(); j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				descs[j].area_histogram[k].normalize(1);
				descs[j].hks_histogram[k].normalize(1);
			}
			float dif_area_arr[3] = {0,0,0};
			float dif_hks_arr[3] = { 0,0,0 };
			for (size_t k = 0; k < 3; k++)
			{
				dif_area_arr[k] = Histogram_ChiSquareDistance(descs[i].area_histogram[k], descs[j].area_histogram[k]);
				dif_hks_arr[k] = Histogram_ChiSquareDistance(descs[i].hks_histogram[k], descs[j].hks_histogram[k]);
				dif_hks_arr[k] = 0; 
				dif_area_arr[k] = VarianceMin_compare(m, descs[i], descs[j], true, 10, 1);
			}
			glm::vec3 dif_area_arr_vec = glm::vec3(dif_area_arr[0], dif_area_arr[1], dif_area_arr[2]);
			glm::vec3 dif_hks_arr_vec = glm::vec3(dif_hks_arr[0], dif_hks_arr[1], dif_hks_arr[2]);
			float dif_area = glm::length(dif_area_arr_vec);
			float dif_hks = glm::length(dif_hks_arr_vec);
			hist_i.push_back(  std::sqrtf( ( dif_area*dif_area)  + (dif_hks * dif_hks)) );
		}
		hist_diffs.push_back(hist_i);
	}

	// maximum n ring 
	float maximum_n_ring = -INFINITY;
	for (size_t i = 0; i < descs.size(); i++)
	{
		for (size_t j = 0; j < descs.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			float dif = std::fabs(descs[i].n_ring_area - descs[j].n_ring_area);
			if (maximum_n_ring < dif)
			{
				maximum_n_ring = dif;
			}
		}
	}
	//maximum area dif
	float maximum_area_dif = -INFINITY;
	for (size_t i = 0; i < descs.size(); i++)
	{
		float area_i = descs[i].area;
		for (size_t j = i; j < descs.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			float area_j = descs[j].area;
			float dif = std::abs(area_i - area_j);
			if (dif > maximum_area_dif)
			{
				maximum_area_dif = dif;
			}

		}
	}

	//maximum skel point dif
	float maximum_skel_point_dist = -INFINITY;
	for (size_t i = 0; i < descs.size(); i++)
	{
		float skel_i = descs[i].skel_point_dist;
		for (size_t j = i; j < descs.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			float skel_j = descs[j].skel_point_dist;
			float dif = std::abs(skel_i - skel_j);
			if (dif > maximum_area_dif)
			{
				maximum_skel_point_dist = dif;
			}

		}
	}

	//maximum gaussian curvature
	/*float maximum_gaussian_curve = -INFINITY;
	for (size_t i = 0; i < descs.size(); i++)
	{
		float gaussian_i = dvorak_pairs[i].gaussian_curv;
		for (size_t j = i; j < descs.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			float gaussian_j = dvorak_pairs[j].gaussian_curv;
			float div = std::abs(gaussian_i / gaussian_j);
			if (div > maximum_gaussian_curve)
			{
				maximum_gaussian_curve = div;
			}

		}
	}*/
	//maximum sdf
	std::vector<float> sdf = computeSDF(m, 20,20);
	auto best_sdf_auto = std::max_element(sdf.begin(), sdf.end());
	int best_sdf_index = (std::distance(sdf.begin(), best_sdf_auto));
	float maximum_sdf = sdf[best_sdf_index];
	
	std::ofstream file("../../Trilateral/Mesh/descriptor.txt");
	for (size_t i = 0; i < descs.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < descs.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}

			int fps_index_i;
			int fps_index_j;
			file << " i " << i << " j " << j << std::endl;
			bool is_hks;
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descs[i].indices[0]] - m->normalized_heat_kernel_signature[descs[j].indices[0]]);
			//is_hks = NLateral_compare_HKS(m , descs[i] , descs[j] , hks_dif_param);
			is_hks = hks_dif < hks_dif_param;
			file << " is hks " << is_hks <<  std::endl; 
			
			
			bool is_area_dif = NLateral_compare_area(m,descs[i],descs[j],maximum_area_dif , area_dif_param,file);
			file << " is area dif " << is_area_dif <<  std::endl;

			//bool is_fuzzy = NLateral_compare_FuzzyGeodesics(m, descs[i], descs[j], fuzzy_param);
 			//float gaussian_curve = std::abs(dvorak_pairs[i].gaussian_curv / dvorak_pairs[j].gaussian_curv);
			//bool is_gaussian = gaussian_curve / maximum_gaussian_curve < curv_param;
			//file << " area dif " << area_dif << " " << is_area_dif << std::endl;
			//bool is_endpoint = Nlateral_check_endpoint(m, skeleton, descs[i], descs[j]);
			bool is_sdf = NLateral_compare_SDF(m, descs[i], descs[j], sdf, sdf_param, file);
			bool is_points_close_to_midpoint = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], mid_point_index, distance_to_mid_param, file);
			bool is_points_far_from_each_other = Nlateral_compare_closeness(m, descs[i], descs[j], mid_point_index, closeness_param, file);
			//bool is_ratio = NLateral_compare_path_ratio(m, descs[i], descs[j], ratio_dif_param, file);
			//file << " is depth " << is_depth << std::endl;
			file << " histogram diff " << hist_diffs[i][j] << std::endl;
			if ( is_hks && is_area_dif && is_points_close_to_midpoint  && is_sdf && is_points_far_from_each_other) /* && is_hks  && is_gaussian && && is_points_close && is_area_dif*/
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = hist_diffs[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}

		}
	}
	std::vector<bool> used(descs.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i]   && !used[j]) {
			resemblance_pairs.push_back({ descs[i].indices[0], descs[j].indices[0] });
			used[i] = true;  // Mark these objects as used
			//used[j] = true;  // Mark these objects as used
		}
	}

	//NLateralMapping_get_best_pairs(m , descs,resemblance_pairs,mid_point_index,mid_point_index_2, 10);
	

	/*for (size_t i = 0; i < descs.size(); i++)
	{
		float smallest_dif = INFINITY;
		int smallest_index = -1; 
		for (size_t j = 0; j < compare_results.size(); j++)
		{
			std::pair<unsigned int, unsigned int> p = compare_results[j].second;
			int i_ = p.first;
			int j_ = p.second;
			if (i  == i_|| i == j_)
			{
				float dif = compare_results[j].first; 
				if (dif < smallest_dif)
				{
					smallest_dif = dif;
					if (i == i_)
						smallest_index = j_;
					else
						smallest_index = i_;
				}
			}
		}
		if (smallest_index != -1)
		{
			resemblance_pairs.push_back({ descs[i].indices[0] ,descs[smallest_index].indices[0] });
		}
	}*/

	file.close();

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
	m->calculated_symmetry_pairs = resemblance_pairs;
	int correct_count = 0;
	/*for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		int ground_truth = m->ground_truth_symmetry_pairs[index1];
		if (index2 == ground_truth)
		{
			correct_count++;
		}
	}

	std::cout << " total == " << correct_count << " " << m->calculated_symmetry_pairs.size() << std::endl;*/
	
	//create descriptors and match thme
	//m->update_raylib_mesh();
	return descs;
}





//using the midpoints
void NLateralMapping_get_best_pairs(TrilateralMesh* m, std::vector<NLateralDescriptor>& descs,
	std::vector<std::pair<unsigned int, unsigned int>>& resemblance_pairs , unsigned int midpoint, unsigned int midpoint_2,
	unsigned int histogram_size) 
{
	std::vector<std::pair<unsigned int, unsigned int>> new_pairs;
	for (size_t i = 0; i < descs.size(); i++)
	{
		unsigned int index = descs[i].indices[0];
		std::vector<std::pair<unsigned int, unsigned int>> pairs_i;
		std::vector<std::pair<float, unsigned int>> hist_values;
		for (size_t j = 0; j < resemblance_pairs.size(); j++)
		{
			if (resemblance_pairs[j].first == index || resemblance_pairs[j].second == index)
			{
				pairs_i.push_back(resemblance_pairs[j]);
			}
		}
		for (size_t j = 0; j < pairs_i.size(); j++)
		{
			unsigned int other_index;
			if (resemblance_pairs[j].first == index)
			{
				other_index = resemblance_pairs[j].second; 
			}
			else
			{
				other_index = resemblance_pairs[j].first;
			}
			std::vector<unsigned int> indices;
			indices.push_back(index);
			indices.push_back(midpoint_2);
			indices.push_back(midpoint);
			NLateralDescriptor desc_i = NLateral_generate_descriptor(m, indices);
			desc_i.create_histogram_area(m, histogram_size,0);
			desc_i.create_histogram_HKS(m, histogram_size,0);
			desc_i.area_histogram[0].normalize(1);
			desc_i.hks_histogram[0].normalize(1);

			indices.clear();
			indices.push_back(other_index);
			indices.push_back(midpoint_2);
			indices.push_back(midpoint);
			NLateralDescriptor desc_i_other = NLateral_generate_descriptor(m, indices);
			desc_i_other.create_histogram_area(m, histogram_size,0);
			desc_i_other.create_histogram_HKS(m, histogram_size,0);
			desc_i_other.area_histogram[0].normalize(1);
			desc_i_other.hks_histogram[0].normalize(1);
			
			
			float dif_area = Histogram_ChiSquareDistance(desc_i.area_histogram[0], desc_i_other.area_histogram[0]);
			float dif_hks = Histogram_ChiSquareDistance(desc_i.hks_histogram[0], desc_i_other.hks_histogram[0]);
			float dif = std::sqrtf((dif_area * dif_area) + (dif_hks * dif_hks));
			hist_values.push_back(std::make_pair(dif,other_index));
		}
		std::sort(hist_values.begin(), hist_values.end());
		if (hist_values.size() > 0)
		{
			new_pairs.push_back(std::make_pair(index, hist_values[0].second));
		}
	}



	resemblance_pairs = new_pairs;
}


std::vector<NLateralDescriptor> NLateralMapping_generate_via_midpoints(TrilateralMesh* m , std::vector<unsigned int>& agd_point_indices, float sweep_distance, float min_geo_tau
,float fuziness , float distance_to_mid_param , float hks_dif_param , float closeness_param )
{
	std::vector<NLateralDescriptor> descs;
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	float biggest_dijkstra = -INFINITY;
	if (agd_point_indices.size() == 0)
	{
		agd_point_indices = Geodesic_avg_dijkstra_modified(m, sweep_distance, 1, false, biggest_dijkstra);
		for (size_t i = 0; i < 3; i++)
		{
			agd_point_indices = Geodesic_min_dijkstra(m, agd_point_indices, sweep_distance, min_geo_tau, false);
		}
	}

	//point_indices = NLateral_sweepDistance(m, point_indices, sweep_distance);
	for (size_t i = 0; i < agd_point_indices.size(); i++)
	{
		int index = agd_point_indices[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 255;
		m->raylib_mesh.colors[index * 4 + 2] = 0;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}
	int hist_size = 5;
	//descs = NLateral_generate_closest_points(m, agd_point_indices,  N, hist_size);
	unsigned int mid_point_index = -1;
	unsigned int mid_point_index_2 = -1;
	Geodesic_mid_point_w_AGD(m, mid_point_index, mid_point_index_2);

	//generate trilateral descriptors from midpoints
	descs = NLateral_generate_with_midpoints(m, agd_point_indices, mid_point_index, mid_point_index_2 , fuziness,biggest_dijkstra);

	std::vector<std::vector<float>> hist_diffs;
	std::ofstream file("../../Trilateral/Mesh/descriptor.txt");
	for (size_t i = 0; i < descs.size(); i++)
	{
		std::vector<float> hist_i;
		for (size_t j = 0; j < 3; j++)
		{
			descs[i].area_histogram[j].normalize(1);
			descs[i].hks_histogram[j].normalize(1);
		}
		for (size_t j = 0; j < descs.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			for (size_t k = 0; k < 3; k++)
			{
				descs[j].area_histogram[k].normalize(1);
				descs[j].hks_histogram[k].normalize(1);
			}
			float dif_area_arr[3] = { 0,0,0 };
			float dif_hks_arr[3] = { 0,0,0 };
			for (size_t k = 0; k < 3; k++)
			{
				dif_area_arr[k] = Histogram_ChiSquareDistance(descs[i].area_histogram[k], descs[j].area_histogram[k]);
				dif_hks_arr[k] = Histogram_ChiSquareDistance(descs[i].hks_histogram[k], descs[j].hks_histogram[k]);
			}
			glm::vec3 dif_area_arr_vec = glm::vec3(dif_area_arr[0], dif_area_arr[1], dif_area_arr[2]);
			glm::vec3 dif_hks_arr_vec = glm::vec3(dif_hks_arr[0], dif_hks_arr[1], dif_hks_arr[2]);
			float dif_area = glm::length(dif_area_arr_vec);
			float dif_hks = glm::length(dif_hks_arr_vec);


			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descs[i].indices[0]] - m->normalized_heat_kernel_signature[descs[j].indices[0]]);
			//is_hks = NLateral_compare_HKS(m , descs[i] , descs[j] , hks_dif_param);
			bool is_hks = hks_dif < hks_dif_param;

			bool is_points_close_to_midpoint = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], mid_point_index, distance_to_mid_param, file);
			bool is_points_far_from_each_other = Nlateral_compare_closeness(m, descs[i], descs[j], mid_point_index, closeness_param, file);
			//bool is_ratio = NLateral_compare_path_ratio(m, descs[i], descs[j], ratio_dif_param, file);
			//file << " is depth " << is_depth << std::endl;
			//file << " histogram diff " << hist_diffs[i][j] << std::endl;
			if (is_hks && is_points_close_to_midpoint  && is_points_far_from_each_other) /* && is_hks  && is_gaussian && && is_points_close && is_area_dif*/
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				compare_results.push_back(std::make_pair(std::sqrtf((dif_area * dif_area) + (dif_hks * dif_hks)), std::make_pair(i, j)));

			}


		}
		hist_diffs.push_back(hist_i);
	}
	std::vector<bool> used(descs.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ descs[i].indices[0], descs[j].indices[0] });
			used[i] = true;  // Mark these objects as used
			//used[j] = true;  // Mark these objects as used
		}
	}
	m->calculated_symmetry_pairs = resemblance_pairs;
	return descs;
}