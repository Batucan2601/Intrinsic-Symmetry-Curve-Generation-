#pragma once 
#include "MeshFactory.h"
#include "CoreTypeDefs.h"

typedef struct
{
	int p_index;
	float gaussian_curv;
}DvorakPairs;
Plane dvorak_generate_plane(MeshFactory& mesh_fac , int selected_index );
std::vector<DvorakPairs> dvorak_extraction_of_significant_points(Mesh* m , int P );
std::vector<DvorakPairs> dvorak_extraction_of_significant_points(Mesh* m, std::vector<unsigned int>& indices);
//criterion
bool dvorak_curvature_similarity_criterion(std::vector<DvorakPairs>& best_pairs , float S ,int index1 , int index2  );
bool dvorak_normal_angle_criterion(Mesh* m, std::vector<DvorakPairs>& best_pairs, int index1, int index2, float c_min_norm);
std::vector<std::pair<int, int>> dvorak_chose_criterion(Mesh* m, std::vector<DvorakPairs>& significant_points, float S, float c_min_norm);
std::vector<Plane> dvorak_generate_candidate_planes(Mesh* m, std::vector<std::pair<int, int>>& best_pairs);
void dvorak_show_signifcant_points(Mesh* m , int P );
float gaussian_curvature(Mesh* mesh, int index);