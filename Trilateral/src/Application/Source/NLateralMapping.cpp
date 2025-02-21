#include "../Include/NLateralMapping.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/Sampling.h"
#include "../Include/MetricCalculations.h"
#include "../Include/VarianceMinimizingTransportPlan.h"
#include "../Include/Geodesic.h"
#include "../Include/ShapeDiameter.h"
#include "../Include/CurvatureGeneration.h"

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
,float fuziness , float distance_to_mid_param , float hks_dif_param , float closeness_param ,float sdf_param ,  int hist_no , int min_agd_param , float& biggest_dijkstra, std::vector<unsigned int>&
original_agd_vertices , float voronoi_param)
{
	std::vector<NLateralDescriptor> descs;
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	std::vector<unsigned int> secondary_curve; 
	if (agd_point_indices.size() == 0)
	{
		agd_point_indices = Geodesic_avg_dijkstra_modified(m, sweep_distance, 1, false, biggest_dijkstra);
		original_agd_vertices = agd_point_indices; 
		for (size_t i = 0; i < min_agd_param; i++)
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
	//descs = NLateral_generate_closest_points(m, agd_point_indices,  N, hist_size);
	unsigned int mid_point_index = -1;
	unsigned int mid_point_index_2 = -1;
	float best_distance; 
	Geodesic_mid_point_w_AGD(m, mid_point_index, mid_point_index_2, best_distance);
	//secondary_curve = Geodesic_generate_secondary_curve_w_midpoints(m,mid_point_index , mid_point_index_2 );

	//generate trilateral descriptors from midpoints
	std::ofstream file("../../Trilateral/Mesh/descriptor.txt");
	descs = NLateral_generate_with_midpoints(m, agd_point_indices, mid_point_index, mid_point_index_2 , fuziness,biggest_dijkstra, hist_no);
	m->sdf = computeSDF(m, 30, 15);
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

			file << " i " << i << " j " << j << std::endl;

			float dif = VarianceMin_compare(m, descs[i], descs[j], true, hist_no, 1);
			file << " dif " << dif << std::endl;

			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descs[i].indices[0]] - m->normalized_heat_kernel_signature[descs[j].indices[0]]);
			bool is_hks = hks_dif < hks_dif_param;

			bool is_points_close_to_midpoint = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], mid_point_index, distance_to_mid_param, file);
			bool is_points_close_to_midpoint_2 = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], mid_point_index_2, distance_to_mid_param, file);
			bool is_points_close_to_midpoint_reverse = NLateral_compare_distance_to_midpoint_reverse(m, descs[i], descs[j],mid_point_index, mid_point_index_2, distance_to_mid_param, file);
			bool is_points_far_from_each_other = Nlateral_compare_closeness(m, descs[i], descs[j], mid_point_index, closeness_param, file);
			bool is_sdf = NLateral_compare_SDF(m, descs[i], descs[j], m->sdf, sdf_param , file);
			//bool is_voronoi = NLateral_compare_voronoi(m, descs[i], descs[j], mid_point_index, 0.15, file);
			//bool is_midpoint_conv = NLateral_compare_divergence(m, descs[i], descs[j], mid_point_index, mid_point_index_2, file);
			//std::vector<unsigned int> path_i_j = conv_int_to_unsigned(Geodesic_between_two_points(*m, descs[i].indices[0], descs[j].indices[0]));
			unsigned int no_of_hit = 0;
			/*bool is_hit = Geodesic_path_intersection(m, secondary_curve, path_i_j, no_of_hit);
			bool is_same_agd = Nlateral_compare_closest_AGD(m, descs[i], descs[j], original_agd_vertices, file);
			if (no_of_hit % 2 == 0)
			{
				is_hit = false; 
			}
			file <<  " no of hit to midpath " << no_of_hit << " is hit " << is_hit << std::endl;
			*/
			//bool is_ratio = NLateral_compare_path_ratio(m, descs[i], descs[j], ratio_dif_param, file);
			//file << " is depth " << is_depth << std::endl;
			//file << " histogram diff " << hist_diffs[i][j] << std::endl;

			if (is_hks  && is_sdf ) // && is_close //&& is_voronoi 
			//&& is_points_far_from_each_other &&  !is_hit  && is_same_agd*/)  
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				//compare_results.push_back(std::make_pair(std::sqrtf((dif_area * dif_area) + (dif_hks * dif_hks)), std::make_pair(i, j)));
				compare_results.push_back(std::make_pair(dif, std::make_pair(i, j)));
			}


		}
		hist_diffs.push_back(hist_i);
	}
	bool first = false; 
	std::vector<bool> used(descs.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ descs[i].indices[0], descs[j].indices[0] });
			used[i] = true;  // Mark these objects as used
			
			if (!first)
			{
				first = !first; 
				std::cout << " the entry with least dist " << i << " " << j << std::endl; 
			}
			//used[j] = true;  // Mark these objects as used
		}
	}
	m->calculated_symmetry_pairs = resemblance_pairs;
	
	std::vector<unsigned int> voronoi_point_set;
	std::vector<unsigned int> points_that_every_curve_has(m->vertices.size(), 0);
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		// part 1 - voronoi boundary
		std::vector<float> distances_1 = Geodesic_dijkstra(*m, resemblance_pairs[i].first);
		std::vector<float> distances_2 = Geodesic_dijkstra(*m, resemblance_pairs[i].second);
		for (size_t j = 0; j < m->vertices.size(); j++)
		{
			float dif = std::abs(distances_1[j] - distances_2[j]);
			if (dif < 1e-2)
			{
				//voronoi_point_set.push_back(j);
				points_that_every_curve_has[j] += 1;
			}
		}

	}
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (points_that_every_curve_has[i] == resemblance_pairs.size())
		{
			voronoi_point_set.push_back(i);
		}
	}
	m->color_points(voronoi_point_set, YELLOW);
	
	
	return descs;
}
// try to add a new one with the highest quality
// quality here being the 
void NLateralMapping_get_new_matchings(TrilateralMesh* m, Curvature& c, std::vector<NLateralDescriptor>& descs, float distance_to_mid_param, float hks_dif_param)
{
	std::ofstream file;
	float best_quality = -INFINITY; 
	std::pair<unsigned int , unsigned int> best_index;

	std::vector<unsigned int> points_far = NLateralMapping_get_unmathced_areas(m , c , 0.1,false); 
	std::vector<NLateralDescriptor> descs_inside;
	for (size_t i = 0; i < points_far.size(); i++)
	{
		for (size_t j = 0; j < descs.size(); j++)
		{
			if (points_far[i] == descs[j].indices[0])
			{
				descs_inside.push_back(descs[j]);
			}
		}
	}
	
	for (size_t i = 0; i < descs_inside.size(); i++)
	{
		for (size_t j = 0; j < descs_inside.size(); j++)
		{
			if (i == j)
			{
				continue; 
			}
			bool is_already_sym = false; 
			for (size_t k = 0; k < m->calculated_symmetry_pairs.size(); k++)
			{
				if (descs_inside[i].indices[0] == m->calculated_symmetry_pairs[k].first)
				{
					is_already_sym = true; 
					break; 
				}
			}
			if (is_already_sym)
			{
				continue;
			}
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descs_inside[i].indices[0]] - m->normalized_heat_kernel_signature[descs_inside[j].indices[0]]);
			bool is_hks = hks_dif < hks_dif_param;
			//bool is_hit = Curvature_curve_intersection_with_two_points(m, c, descs[i].indices[0], descs[j].indices[0], 0.85);
			bool is_points_close_to_midpoint = NLateral_compare_distance_to_midpoint(m, descs_inside[i], descs_inside[j], descs_inside[0].indices[1], distance_to_mid_param, file);
			bool is_points_close_to_midpoint_2 = NLateral_compare_distance_to_midpoint(m, descs_inside[i], descs_inside[j], descs_inside[0].indices[2], distance_to_mid_param, file);
			//bool is_far_away_from_others = NLateral_check_far_away(m, descs_inside[i], descs_inside[j], 0.05);
			
			if ( !( is_points_close_to_midpoint && is_points_close_to_midpoint_2 ))
			{
				continue; 
			}
			unsigned int midpoint = Geodesic_get_midpoint_from_path(m, descs_inside[i].indices[0], descs_inside[j].indices[0]);
			std::vector<float> distances = Geodesic_dijkstra(*m, midpoint);
			float quality_k = 0;
			for (size_t k = 0; k < m->calculated_symmetry_pairs.size(); k++)
			{
				float dist_from_first = distances[m->calculated_symmetry_pairs[k].first];
				float dist_from_second = distances[m->calculated_symmetry_pairs[k].second];
				
				quality_k += std::min(dist_from_first , dist_from_second) / std::max(dist_from_first, dist_from_second);
			}
			
			if (best_quality < quality_k)
			{
				best_quality = quality_k;
				best_index = std::make_pair(i,j);
			}
		}
	}

	m->calculated_symmetry_pairs.push_back(std::make_pair(best_index.first, best_index.second));
	CurvePoints cp;
	cp.mid_point = Geodesic_get_midpoint_from_path(m, descs_inside[best_index.first].indices[0], descs_inside[best_index.second].indices[0]);
	cp.correspondence = std::make_pair(descs_inside[best_index.first].indices[0], descs_inside[best_index.second].indices[0]);
	c.curve_points.push_back(cp);
	build_curvature(m, c);
}


void NLateralMapping_get_new_matchings(TrilateralMesh* m, Curvature& c , std::vector<NLateralDescriptor>& descs, float distance_to_mid_param ,float sdf_param ,
float hks_dif_param , int hist_no )
{
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int>>> compre_results ;
	std::ofstream file; 
	for (size_t i = 0; i < descs.size(); i++)
	{
		for (size_t j = 0; j < descs.size(); j++)
		{
			bool is_continue = true; 
			if (i == j)
			{
				continue; 
			}
			//also check if already together
			for (size_t k = 0; k < m->calculated_symmetry_pairs.size(); k++)
			{
				unsigned int desc1_i0 = m->calculated_symmetry_pairs[k].first; 
				if (desc1_i0 == descs[i].indices[0])
				{
					is_continue = false;
					break; 
				}
			}
			if (!is_continue)
			{
				continue; 
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

			file << " i " << i << " j " << j << std::endl;

			float dif = VarianceMin_compare(m, descs[i], descs[j], true, hist_no, 1);

			unsigned int no_of_hit = 0;
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[descs[i].indices[0]] - m->normalized_heat_kernel_signature[descs[j].indices[0]]);
			bool is_hks = hks_dif < hks_dif_param;
			bool is_hit = Curvature_curve_intersection_with_two_points(m, c, descs[i].indices[0], descs[j].indices[0], 0.85);
			bool is_points_close_to_midpoint = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], descs[0].indices[1], distance_to_mid_param, file);
			bool is_points_close_to_midpoint_2 = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], descs[0].indices[2], distance_to_mid_param, file);
			//bool is_new_curve_makes_sense = Curvature_compare_with_new_addition(m , c ,  descs[i], descs[j] );
			bool is_sdf = NLateral_compare_SDF(m, descs[i], descs[j], m->sdf, sdf_param, file);
			//bool is_voronoi = NLateral_compare_voronoi(m, descs[i], descs[j], descs[0].indices[1], 0.25, file);

			if (is_hks && is_points_close_to_midpoint && is_points_close_to_midpoint_2 && is_sdf )
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				//compare_results.push_back(std::make_pair(std::sqrtf((dif_area * dif_area) + (dif_hks * dif_hks)), std::make_pair(i, j)));
				compre_results.push_back(std::make_pair(dif, std::make_pair(i, j)));
			}
		}
	}
	std::vector<bool> used(descs.size(), false);
	std::vector<std::pair<unsigned int, unsigned int> > resemblance_pairs;
	std::sort(compre_results.begin(), compre_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compre_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			resemblance_pairs.push_back({ descs[i].indices[0], descs[j].indices[0] });
			break; // di it once
			used[i] = true;  // Mark these objects as used
			//used[j] = true;  // Mark these objects as used
		}
	}

	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int no_of_hits = 0;
		bool is_hit = Curvature_curve_intersection_with_two_points(m, c, m->calculated_symmetry_pairs[i].first
			, m->calculated_symmetry_pairs[i].second, no_of_hits);
		std::cout << " no of hits " << no_of_hits << std::endl; 
	}

	//append
	for (size_t i = 0; i < resemblance_pairs.size(); i++)
	{
		m->calculated_symmetry_pairs.push_back(resemblance_pairs[i]);
	}
}

std::vector<unsigned int> NLateralMapping_get_unmathced_areas(TrilateralMesh* m, Curvature& c, float param , bool is_color)
{
	std::vector<unsigned int> colored_points;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		std::vector<float> distances = Geodesic_dijkstra(*m, i);
		auto max_element_auto = std::max_element(distances.begin(), distances.end());
		float max_elem = *max_element_auto;
		bool is_far = true;
		for (size_t j = 0; j < m->calculated_symmetry_pairs.size(); j++)
		{
			float dist = param * max_elem;
			int index_1 = m->calculated_symmetry_pairs[j].first;
			int index_2 = m->calculated_symmetry_pairs[j].second;
			if (distances[index_1] < dist || distances[index_2] < dist)
			{
				is_far = false; 
				break; 
			}
		}
		if (is_far)
		{
			colored_points.push_back(i);
		}
	}
	if (is_color)
	{
		m->color_points(colored_points, RED);
	}
	return colored_points; 
}