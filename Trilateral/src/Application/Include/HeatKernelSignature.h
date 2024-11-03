#pragma once 
#include "../Include/Laplace-Beltrami.h"
#include "../Include/MeshFactory.h"
#include "../Include/TrilateralDescriptor.h"

void HKS_extract_kernel_signature(TrilateralMesh* m);

void HKS_read_kernel_signature(TrilateralMesh* m);

std::vector<std::vector<double>> HKS_compute_kernel(TrilateralMesh* m, std::pair<Eigen::VectorXd, Eigen::MatrixXd>& eigen_pairs, const std::vector<double>& timeSteps
, int time_step_no);

void HKS_hks_on_descriptor(TrilateralMesh* m, TrilateralDescriptor& desc);