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

float generate_symmetry_score(Mesh mesh, Plane* p1);
Plane generate_isomap_embedding(Mesh* mesh, bool simplify_mesh , float simplification_percentage );
Eigen::MatrixXd ComputeClassicalMds(const Eigen::MatrixXd& D, const unsigned target_dim);
Plane generate_symmetry_plane_dividing_classical_MDS(Mesh* mesh);
Eigen::MatrixXd compute_landmark_MDS(Mesh* mesh ,  const unsigned target_dim , const int no_of_landmarks = 100 );