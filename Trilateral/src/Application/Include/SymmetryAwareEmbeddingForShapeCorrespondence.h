#pragma once

#include "../Include/MeshFactory.h"
#include <GL/glew.h>
#include <utility>
#include <algorithm>
#include <queue>
#include <stack>
#include <algorithm> 
#include <math.h>
#include <eigen/Eigen/Dense>
#include "Sampling.h"
#include "CoreTypeDefs.h"
#include <src/Application/Include/CoreTypeDefs.h>
// https://geometry.stanford.edu/papers/yyg-saefsc-16/yyg-saefsc-16.pdf

float generate_symmetry_score(Mesh mesh, Plane* p1 , bool is_simplify_active, std::vector<bool>* is_vertex_exist_vector  );
Plane generate_isomap_embedding(Mesh* mesh, bool simplify_mesh , float simplification_percentage );
Eigen::MatrixXd ComputeClassicalMds(const Eigen::MatrixXd& D, const unsigned target_dim);
Plane generate_symmetry_plane_dividing_classical_MDS(Mesh* mesh);
Mesh compute_landmark_MDS(Mesh* mesh ,  const unsigned target_dim , const int no_of_landmarks = 100 );
void trilateral_symmetry_with_landmark_MDS_with_plane(Mesh* mesh ,  const unsigned target_dim , const int no_of_landmarks = 100 );
Plane trilateral_symmetry_with_landmark_MDS_with_plane(Mesh* mesh, const unsigned target_dim, const int no_of_landmarks, const int no_of_trilateral_points , float& error_percentage);
void create_trilateral_sym_w_landmarl_with_planes(std::vector<Mesh> mesh, const unsigned target_dim, const int no_of_landmarks, const int no_of_trilateral_points , std::string filename);
