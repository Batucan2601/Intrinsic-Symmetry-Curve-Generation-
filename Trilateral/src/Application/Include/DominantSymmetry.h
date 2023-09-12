#pragma once 
#include "MeshFactory.h"
#include "CoreTypeDefs.h"
Plane generate_dominant_symmetry_plane(int seletected_mesh , MeshFactory & mesh_fac );
void generate_two_separate_mesh_using_dominant_symmetry_plane(Plane plane, Mesh* mesh_to_be_separated, Mesh* m1, Mesh* m2, std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2);
void match_two_meshes_with_fps(Mesh* m1, Mesh* m2 ,std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2 , int no_of_samples);