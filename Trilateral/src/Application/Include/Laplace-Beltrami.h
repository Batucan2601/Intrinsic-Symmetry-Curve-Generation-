#pragma once 
#include "eigen/Eigen/dense"
#include "eigen/Eigen/sparse"
#include "MeshFactory.h"
#include "glm/glm.hpp"
#include "Skeleton.h"
#include "../Include/NLateralDescriptor.h"
//#include "spectra/Spectra/SymEigsSolver.h"

// chatgpt
Eigen::MatrixXd embed_mesh_to_2d(TrilateralMesh& mesh); 
Eigen::MatrixXd embed_mesh_endpoints_to_2d(TrilateralMesh& mesh, Skeleton& skeleton, NLateralParameters& nLateralParameters);

Eigen::MatrixXd generate_L(TrilateralMesh * mesh  );
Eigen::MatrixXd generate_W_cotangent_laplacian(TrilateralMesh* mesh);
Eigen::MatrixXd generate_A(TrilateralMesh* mesh);

Eigen::SparseMatrix<double> cotangent_laplacian(const TrilateralMesh& mesh);

bool check_if_matrix_symmetric(const Eigen::SparseMatrix<double>& mat);
Eigen::SparseMatrix<double> normalize_laplacian(const Eigen::SparseMatrix<double>& L);
Eigen::SparseMatrix<double> regularize_matrix(const Eigen::SparseMatrix<double>& L, double epsilon);