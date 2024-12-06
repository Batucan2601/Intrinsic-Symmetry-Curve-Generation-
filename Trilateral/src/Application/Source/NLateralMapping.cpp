#include "../Include/NLateralMapping.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/Sampling.h"
#include "../Include/MetricCalculations.h"
#include "../Include/VarianceMinimizingTransportPlan.h"
#include "../Include/Geodesic.h"
#include "../Include/ShapeDiameter.h"
std::pair<std::vector<NLateralDescriptor>,std::vector<NLateralDescriptor>> NlateralMap_point_matching_copy_symmetric_points(TrilateralMesh* m, Skeleton& skeleton,Plane& plane,
	int dvorak_enpoint_no, float sweep_distance, float hks_dif_param, float curv_param, float norm_angle_param, float skel_dist_param, float n_ring_param,
	float area_dif_param,float skel_point_dist_param, int N)
{

	int size = m->vertices.size();
	int mesh_mid_point_index = -1;
	std::vector<NLateralDescriptor> desc_neg;
	std::vector<NLateralDescriptor> desc_pos;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<std::pair<int, float>> hks_pairs;
	std::vector<DvorakPairs> dvorak_pairs;
	std::vector<unsigned int> fps_indices;
	std::vector<unsigned int> fps_sym_indices;
	std::vector<unsigned int> points_pos;
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;
	//hks_pairs = HKS_extraction_significant_points(m, dvorak_enpoint_no);
	//sweep 
	//hks_pairs = HKS_sweep_distance(m, hks_pairs, sweep_distance);
	//dvorak_pairs = HKS_to_dvorak_pairs(m, hks_pairs);
	

	Metric_set_gaussian(m, dvorak_enpoint_no, sweep_distance);
	Metric_set_N(N);
	if (plane.isNull())
	{
		plane = generate_dominant_symmetry_plane(glm::vec3(0, 0, 0), *m); //dont use this
	}
	//get one side
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (get_point_status_from_plane(&plane, &(m->vertices[i])) > 0   )
		{
			points_pos.push_back(i);
		}
	}
	fps_indices = furthest_point_sampling_on_partial_points(m, dvorak_enpoint_no, points_pos);
	
	

	for (size_t i = 0; i < dvorak_enpoint_no; i++)
	{
		int index = fps_indices[i];
		int sym = m->ground_truth_symmetry_pairs[index];
		m->raylib_mesh.colors[sym * 4] = 0;
		m->raylib_mesh.colors[sym * 4 + 1] = 255;
		m->raylib_mesh.colors[sym * 4 + 2] = 0; 
		m->raylib_mesh.colors[sym * 4 + 3] = 255;
		fps_sym_indices.push_back(sym);
	}
	SkeletonTree skelTree = skeleton_generate_skeleton_tree(m, skeleton);
	int hist_size = 10; 
	desc_pos = NLateral_generate_closest_points(m, skeleton, fps_indices, skelTree, N, 500 , hist_size);
	desc_neg = NLateral_generate_closest_points(m, skeleton, fps_sym_indices, skelTree, N, 500, hist_size);


	for (size_t i = 0; i < fps_indices.size(); i++)
	{
		int index = fps_indices[i];
		DvorakPairs p;
		p.p_index = index;
		p.gaussian_curv = gaussian_curvature(m, index);
		dvorak_pairs.push_back(p);
	}
	for (size_t i = 0; i < fps_sym_indices.size(); i++)
	{
		int index = fps_sym_indices[i];
		DvorakPairs p;
		p.p_index = index;
		p.gaussian_curv = gaussian_curvature(m, index);
		dvorak_pairs.push_back(p);
	}

	//std::vector<std::vector<float>> hist_diffs=  VarianceMin_compare_all(m, desc_pos, desc_neg,true,10,1 );
	std::vector<std::vector<float>> hist_diffs;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		std::vector<float> hist_i;
		desc_pos[i].histogram.normalize(1);
		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			desc_neg[j].histogram.normalize(1);
			float dif = Histogram_ChiSquareDistance( desc_pos[i].histogram, desc_neg[j].histogram);

			hist_i.push_back(dif);
		}
		hist_diffs.push_back(hist_i);
	}

	// maximum n ring 
	float maximum_n_ring = -INFINITY;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			float dif = std::fabs(desc_pos[i].n_ring_area - desc_neg[j].n_ring_area);
			if (maximum_n_ring < dif)
			{
				maximum_n_ring = dif;
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
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		float area_i = desc_pos[i].area;
		for (size_t j = i; j < desc_neg.size(); j++)
		{
			float area_j = desc_neg[j].area;
			float dif = std::abs(area_i - area_j);
			if (dif > maximum_area_dif)
			{
				maximum_area_dif = dif;
			}

		}
	}

	//maximum skel point dif
	float maximum_skel_point_dist = -INFINITY; 
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		float skel_i = desc_pos[i].skel_point_dist;
		for (size_t j = i; j < desc_neg.size(); j++)
		{
			float skel_j = desc_neg[j].skel_point_dist;
			float dif = std::abs(skel_i - skel_j);
			if (dif > maximum_area_dif)
			{
				maximum_skel_point_dist = dif;
			}

		}
	}

	//maximum gaussian curvature
	float maximum_gaussian_curve = -INFINITY;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		float gaussian_i = dvorak_pairs[i].gaussian_curv;
		for (size_t j = i; j < desc_neg.size(); j++)
		{
			float gaussian_j = dvorak_pairs[j + desc_pos.size()].gaussian_curv;
			float div = std::abs(gaussian_i / gaussian_j);
			if (div > maximum_gaussian_curve)
			{
				maximum_gaussian_curve = div;
			}

		}
	}
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		float smallest = INFINITY;
		int index = -1;
		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			int fps_index_i;
			int fps_index_j;
			float hks_dif = std::abs(m->normalized_heat_kernel_signature[desc_pos[i].indices[0]] - m->normalized_heat_kernel_signature[desc_neg[j].indices[0]]);
			bool is_hks = hks_dif < hks_dif_param;
			bool is_skel_dist_far = 1;// NLateral_compare_skeldist_mid(m, desc_pos[i], desc_pos[j], skel_point_dist_param, maximum_skel_dist);
			float n_ring = std::abs(desc_pos[i].n_ring_area - desc_neg[j].n_ring_area);
			bool is_n_ring_close = (n_ring / maximum_n_ring) < n_ring_param;
			float area_dif = std::abs(desc_pos[i].area - desc_neg[j].area);
			bool is_area_dif = area_dif / maximum_area_dif < area_dif_param;
			float skel_point_dist_dif = std::abs(desc_pos[i].skel_point_dist- desc_neg[j].skel_point_dist);
			bool is_skel_point_dist = skel_point_dist_dif / maximum_skel_point_dist < skel_point_dist_param; 
			
			float gaussian_curve= std::abs(dvorak_pairs[i].gaussian_curv / dvorak_pairs[j + desc_pos.size()].gaussian_curv);
			bool is_gaussian = gaussian_curve / maximum_gaussian_curve < curv_param;
			std::cout << "hks " <<  is_hks << std::endl;
			std::cout << "skel dist " << is_skel_dist_far <<  std::endl;
			std::cout << "n ring " <<  is_n_ring_close << std::endl;
			std::cout << "area dif " <<  is_area_dif << std::endl;
			std::cout << "gaussian dif " <<  is_area_dif << std::endl;
			if (  is_hks && is_n_ring_close && is_area_dif && is_skel_dist_far && is_skel_point_dist /* && is_gaussian*/)
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = hist_diffs[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}

			/*std::cout << " i " << i << " j " << j << " " << hist_diffs[i][j] << std::endl;
			std::pair<float, std::pair<unsigned int, unsigned int>> res;
			res.first = hist_diffs[i][j];
			res.second = std::make_pair(i, j);
			compare_results.push_back(res);*/
			
		}
	}
	std::vector<bool> used(desc_pos.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] /* && !used[j]*/) {
			resemblance_pairs.push_back({ desc_pos[i].indices[0], desc_neg[j].indices[0] });
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

	m->calculated_symmetry_pairs = resemblance_pairs;
	
	int correct_count = 0;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = m->calculated_symmetry_pairs[i].first;
		int index2 = m->calculated_symmetry_pairs[i].second;
		int ground_truth = m->ground_truth_symmetry_pairs[index1];
		if (index2 == ground_truth)
		{
			correct_count++;
		}
	}

	std::cout << " total == " << correct_count << " " << m->calculated_symmetry_pairs.size() << std::endl;
	//create descriptors and match thme
	m->update_raylib_mesh();
	return std::make_pair(desc_pos,desc_neg);
}


std::vector<NLateralDescriptor> NlateralMap_point_matching_w_average_geodesic(TrilateralMesh* m, Skeleton& skeleton, 
	int dvorak_enpoint_no, float sweep_distance, float hks_dif_param, float curv_param, float norm_angle_param, float skel_dist_param, float ratio_dif_param,
	float area_dif_param, float skel_point_dist_param,float paths_dif_param,float min_geo_tau,int avg_n_ring, int skel_depth_param,float tri_hist_param,
	float distance_to_mid_param , float sdf_param , int N)
{
	int size = m->vertices.size();
	int mesh_mid_point_index = -1;
	std::vector<NLateralDescriptor> descs;
	glm::vec3 mesh_mid_point;
	// 1 - get dvork significant points
	std::vector<std::pair<int, float>> hks_pairs;
	std::vector<DvorakPairs> dvorak_pairs;
	std::vector<unsigned int> point_indices;
	std::vector<unsigned int> points_pos;
	std::vector<std::pair<float, std::pair<unsigned int, unsigned int> >> compare_results;
	std::vector<std::pair<unsigned int, unsigned int>> resemblance_pairs;

	Metric_set_gaussian(m, dvorak_enpoint_no, sweep_distance);
	Metric_set_N(N);

	point_indices = Geodesic_avg_dijkstra_modified(m, dvorak_enpoint_no, sweep_distance, avg_n_ring, true);
	for (size_t i = 0; i < 3; i++)
	{
		point_indices = Geodesic_min_dijkstra(m, dvorak_enpoint_no, point_indices, sweep_distance, min_geo_tau, true);
	}

	//point_indices = NLateral_sweepDistance(m, point_indices, sweep_distance);
	for (size_t i = 0; i < point_indices.size(); i++)
	{
		int index = point_indices[i];
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 255;
		m->raylib_mesh.colors[index * 4 + 2] = 0;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}
	SkeletonTree skelTree = skeleton_generate_skeleton_tree(m, skeleton);
	int hist_size = 10;
	descs = NLateral_generate_closest_points(m, skeleton, point_indices, skelTree, N, 500, hist_size);


	for (size_t i = 0; i < point_indices.size(); i++)
	{
		int index = point_indices[i];
		DvorakPairs p;
		p.p_index = index;
		p.gaussian_curv = gaussian_curvature(m, index);
		dvorak_pairs.push_back(p);
	}


	//std::vector<std::vector<float>> hist_diffs=  VarianceMin_compare_all(m, descs,true,10,1 );
	std::vector<std::vector<float>> hist_diffs;
	for (size_t i = 0; i < descs.size(); i++)
	{
		std::vector<float> hist_i;
		descs[i].histogram.normalize(1);
		for (size_t j = 0; j < descs.size(); j++)
		{
			descs[j].histogram.normalize(1);
			float dif = Histogram_ChiSquareDistance(descs[i].histogram, descs[j].histogram);
			//float dif = Histogram_ChiSquareDistance(descs[i].histogram, descs[j].histogram);

			hist_i.push_back(dif);
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
	//maximum skel distance
	float maximum_skel_dist= -INFINITY;
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
	float maximum_gaussian_curve = -INFINITY;
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
	}
	//maximum sdf
	float maximum_sdf = ShapeDiameter_calculate_simple_max_dif(m, point_indices);
	unsigned int mid_point_index = NLateral_get_closest_index_to_midpoint(m, point_indices);
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
			//is_hks = NLateral_compare_HKS(m,descs[i], descs[j], hks_dif_param,file);
			is_hks = hks_dif < hks_dif_param;
			file << " is hks " << is_hks <<  std::endl; 
			//bool is_hks = hks_dif < hks_dif_param;
			bool is_skel_dist_far = NLateral_compare_skeldist_mid(m, descs[i], descs[j], skel_dist_param, maximum_skel_dist,file);
			file << " is skel dist  " <<  is_skel_dist_far << std::endl;

			//float n_ring = std::abs(descs[i].n_ring_area - descs[j].n_ring_area);
			//bool is_n_ring_close = (n_ring / maximum_n_ring) < n_ring_param;
			//float area_dif = std::abs(descs[i].area - descs[j].area);
			//bool is_area_dif = area_dif / maximum_area_dif < area_dif_param;
			bool is_area_dif = NLateral_compare_area(m,descs[i],descs[j],maximum_area_dif , area_dif_param,file);
			file << " is area dif " << is_area_dif <<  std::endl;

			float skel_point_dist_dif = std::abs(descs[i].skel_point_dist - descs[j].skel_point_dist);
			file << " skel point dist " << std::endl;
			file << skel_point_dist_dif << std::endl;
			bool is_skel_point_dist = skel_point_dist_dif / maximum_skel_point_dist < skel_point_dist_param;
			file << " is skel point dist " << is_skel_point_dist;
 
			float ratio_dif = std::fabs(descs[i].paths_ratio - descs[j].paths_ratio);
			bool is_ratio_dif = ratio_dif < ratio_dif_param;
			file << "  desc path ratio dif " << ratio_dif << std::endl;
			file << "  desc path ratio  " << is_ratio_dif << std::endl;

 			float gaussian_curve = std::abs(dvorak_pairs[i].gaussian_curv / dvorak_pairs[j].gaussian_curv);
			bool is_gaussian = gaussian_curve / maximum_gaussian_curve < curv_param;
			//file << " area dif " << area_dif << " " << is_area_dif << std::endl;
			bool is_depth = NLatera_compare_depth(m, descs[i], descs[j], skel_depth_param, file);
			bool is_endpoint = Nlateral_check_endpoint(m, skeleton, descs[i], descs[j]);
			bool is_nlateral_dist_midpoint = NLateral_compare_distance_to_midpoint(m, descs[i], descs[j], mid_point_index,
			distance_to_mid_param, file);
			bool is_sdf = NLateral_compare_SDF(m, descs[i], descs[j], maximum_sdf, sdf_param, file);
			
			//file << " is depth " << is_depth << std::endl;
			if (is_hks && is_sdf && is_skel_dist_far /* && is_depth */ && is_endpoint && is_nlateral_dist_midpoint && is_gaussian
			&& is_sdf)
			{
				std::pair<float, std::pair<unsigned int, unsigned int>> res;
				res.first = hist_diffs[i][j];
				res.second = std::make_pair(i, j);
				compare_results.push_back(res);
			}

			/*std::cout << " i " << i << " j " << j << " " << hist_diffs[i][j] << std::endl;
			std::pair<float, std::pair<unsigned int, unsigned int>> res;
			res.first = hist_diffs[i][j];
			res.second = std::make_pair(i, j);
			compare_results.push_back(res);*/

		}
	}
	file.close();
	std::vector<bool> used(descs.size(), false);
	std::sort(compare_results.begin(), compare_results.end());
	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compare_results) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] /* && !used[j]*/) {
			resemblance_pairs.push_back({ descs[i].indices[0], descs[j].indices[0] });
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
	m->update_raylib_mesh();
	return descs;
}