#pragma once 
#include "../Include/TrilateralMesh.h"
#include "../Include/CoreTypeDefs.h"
#include <vector>


void KIDS_read_meshes();
void KIDS_dom_sym_generate_or_read_planes(float convergence_ratio);
void KIDS_generate_gaussians(int no_of_gaussian, float sweep_distance);
void KIDS_select_mesh(TrilateralMesh& m, int meshNo);
void KIDS_endpoint_matching_w_gaussian(int gaussian_endpoint, float convergence_ratio);
void KIDS_endpoint_matching_w_NLateral(int gaussian_end_point_no, float sweep_distance, int N );
void KIDS_save_descriptors(int index, float convergence_ratio);