#include "../Include/NLateralMapping.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/Sampling.h"
#include "../Include/MetricCalculations.h"
#include "../Include/VarianceMinimizingTransportPlan.h"
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
			float skel_dist = std::abs(desc_pos[i].skel_dist_mid - desc_neg[j].skel_dist_mid);
			bool is_skel_dist_far = (skel_dist / maximum_skel_dist) < skel_dist_param;
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