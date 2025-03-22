#pragma once
#include "TrilateralMesh.h"
#include "NLateralDescriptor.h"
struct Voronoi
{
	Voronoi();
	Voronoi(TrilateralMesh* m, unsigned int p1, unsigned int p2, float param);
	TrilateralMesh* m;
	std::vector<unsigned int> indices;
	float param;
	int no_of_partition = 0;
	float distance_to_index(unsigned int index); //mesh index
	std::vector<int> status; // if on right 0 if on voronoi -1 other part1 
	void generate_voronoi_parts();
	void color();
	void connect_boundary();
	int p1, p2;

	std::vector<unsigned int> get_closest_path(unsigned int index1, unsigned int index2);
};

float NLateral_generate_descriptors_with_random_voronoi_points(
	TrilateralMesh* m, unsigned int p1, unsigned int p2, float voronoi_param, float fuziness, int hist_no
	, int no_of_point_param);

Voronoi Voronoi_get_closest_voronoi(TrilateralMesh* m, float voronoi_param);

Voronoi Voronoi_destroy_wrong_matches_and_recalculate(TrilateralMesh* m, float voronoi_param, Voronoi& v);
Voronoi Voronoi_check_pair_closeness_and_recalculate(TrilateralMesh* m, float voronoi_param, float dist_to_voronoi_param, Voronoi& v);
Voronoi Voronoi_algorithm_in_action(TrilateralMesh* m, float voronoi_param, float hks_param, float sdf_param, float dist_param, int hist_no, float fuziness, int no_of_points, std::vector<unsigned int>& agd_points);
Voronoi Voronoi_algorithm_in_action(TrilateralMesh* m,Voronoi& v,  float voronoi_param, float hks_param, float sdf_param, float dist_param, int hist_no, float fuziness, int no_of_points, std::vector<unsigned int>& agd_points);
void Voronoi_add_point(TrilateralMesh* m, Voronoi& voronoi, int selected_agd_index, float voronoi_param, float hks_param, float sdf_param , float dist_param, unsigned int midpoint, int hist_no, float fuziness
	, std::vector<unsigned int>& agd_points);

Voronoi Voronoi_add_points_and_recalculate(TrilateralMesh* m, Voronoi& v, float voronoi_param, float hks_param,float sdf_param, float dist_param, int hist_no, float fuziness, int no_of_points
	, std::vector<unsigned int>& agd_points);

void Voronoi_color_every_pairs_voronoi(TrilateralMesh* m, float voronoi_param);
Voronoi Voronoi_show_voronoi(TrilateralMesh* m, unsigned int pair_no, float voronoi_param);

void Voronoi_prune_voronoi(TrilateralMesh* m, Voronoi& voronoi, float voronoi_param);


void Voronoi_get_two_pair_and_generate_a_path(TrilateralMesh* m, Voronoi& voronoi, int voronoi_p1 ,  float fuzziness, int voronoi_division_no);