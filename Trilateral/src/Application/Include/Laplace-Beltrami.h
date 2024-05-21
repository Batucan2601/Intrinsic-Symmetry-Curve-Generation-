#pragma once 
#include "eigen/Eigen/dense"
#include "MeshFactory.h"
#include "glm/glm.hpp"
#include "Skeleton.h"
#include "../Include/NLateralDescriptor.h"
//#include "spectra/Spectra/SymEigsSolver.h"

// chatgpt
Eigen::MatrixXd embed_mesh_to_2d(Mesh& mesh); 
Eigen::MatrixXd embed_mesh_endpoints_to_2d(Mesh& mesh, Skeleton& skeleton, NLateralParameters& nLateralParameters);

Eigen::MatrixXd generate_L(Mesh * mesh  );
Eigen::MatrixXd generate_W_cotangent_laplacian(Mesh* mesh);
Eigen::MatrixXd generate_A(Mesh* mesh);