#pragma once 
#include "MeshFactory.h"
#include "TrilateralDescriptor.h"
#include "CoreTypeDefs.h"
Plane generate_dominant_symmetry_plane(int seletected_mesh , MeshFactory & mesh_fac , float sym_plane_iter);
Plane generate_dominant_symmetry_plane(TrilateralMesh* mesh , float sym_plane_iter);
Plane generate_dominant_symmetry_plane(const glm::vec3& plane_point , TrilateralMesh mesh);

void dom_sym_generate_two_separate_mesh_using_dominant_symmetry_plane(Plane plane, TrilateralMesh* mesh_to_be_separated, TrilateralMesh* m1, TrilateralMesh* m2, std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2);
void match_two_meshes_with_fps(TrilateralMesh* selected_mesh, TrilateralMesh* m1, TrilateralMesh* m2, std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2, int no_of_samples);
std::vector<TrilateralDescriptor> match_two_meshes_with_fps(TrilateralMesh* selected_mesh, Plane* plane, int no_of_samples);

void dom_sym_write_plane(TrilateralMesh* m, Plane& plane, std::string path);
bool dom_sym_read_plane(TrilateralMesh* m , Plane& plane , std::string path);
void dom_sym_draw_plane();
