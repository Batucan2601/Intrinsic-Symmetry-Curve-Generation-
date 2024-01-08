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
