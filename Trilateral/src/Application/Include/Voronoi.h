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
	float distance_to_index(unsigned int index); //mesh index
	std::vector<int> status; // if on right 0 if on voronoi -1 other part1 
	void generate_voronoi_parts();
	void color();

};

float NLateral_generate_descriptors_with_random_voronoi_points(
	TrilateralMesh* m, unsigned int p1, unsigned int p2, float voronoi_param, float fuziness, int hist_no
	, int no_of_point_param);

Voronoi Voronoi_get_closest_voronoi(TrilateralMesh* m, float voronoi_param);

Voronoi Voronoi_destroy_wrong_matches_and_recalculate(TrilateralMesh* m, float voronoi_param, Voronoi& v);
Voronoi Voronoi_check_pair_closeness_and_recalculate(TrilateralMesh* m, float voronoi_param, float dist_to_voronoi_param, Voronoi& v);
