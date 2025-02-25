#include "../Include/NLateralDescriptor.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/Skeleton.h"
#include "../Include/CurvatureGeneration.h"

void NLateralMapping_get_best_pairs(TrilateralMesh* m, std::vector<NLateralDescriptor>& descs,
	std::vector<std::pair<unsigned int, unsigned int>>& resemblance_pairs, unsigned int midpoint,
	unsigned int midpoint_2, unsigned int histogram_size);

std::vector<NLateralDescriptor> NLateralMapping_generate_via_voronoi_midpoints(TrilateralMesh* m, std::vector<unsigned int>& agd_point_indices, float sweep_distance, float min_geo_tau
	, float fuziness, float distance_to_mid_param, float hks_dif_param, float closeness_param, float sdf_param, int hist_no, int min_agd_param, float& biggest_dijkstra, std::vector<unsigned int>&
	original_agd_vertices, float voronoi_param);
void NLateralMapping_get_new_matchings(TrilateralMesh* m, Curvature& c, std::vector<NLateralDescriptor>& descs, float distance_to_mid_param, float sdf_param,
	float hks_dif_param, int hist_no);

void NLateralMapping_get_new_matchings(TrilateralMesh* m, Curvature& c, std::vector<NLateralDescriptor>& descs, float distance_to_mid_param, float hks_dif_param);

std::vector<unsigned int> NLateralMapping_get_unmathced_areas(TrilateralMesh* m, Curvature& c, float param, bool is_color);

bool NLateral_check_voronoi_area(TrilateralMesh* m, unsigned int p1, unsigned int p2, float voronoi_param, float area_dif_param);
