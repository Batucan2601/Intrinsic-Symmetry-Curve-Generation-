#pragma once 
#include "MeshFactory.h"
#include "CoreTypeDefs.h"
Plane generate_dominant_symmetry_plane(int seletected_mesh , MeshFactory & mesh_fac , float sym_plane_iter);
Plane generate_dominant_symmetry_plane(Mesh mesh , float sym_plane_iter);
Plane generate_dominant_symmetry_plane(const glm::vec3& plane_point , Mesh mesh);

void dom_sym_generate_two_separate_mesh_using_dominant_symmetry_plane(Plane plane, Mesh* mesh_to_be_separated, Mesh* m1, Mesh* m2, std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2);
void match_two_meshes_with_fps(Mesh* selected_mesh, Mesh* m1, Mesh* m2, std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2, int no_of_samples);
std::vector<TrilateralDescriptor> match_two_meshes_with_fps(Mesh* selected_mesh, Plane* plane, int no_of_samples);

void dom_sym_save_plane(Plane& plane, Mesh* m);
bool dom_sym_read_plane(MeshFactory& mesh_fac, int selected_mesh, Plane& plane);
