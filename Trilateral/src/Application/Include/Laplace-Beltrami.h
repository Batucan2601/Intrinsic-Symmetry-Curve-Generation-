#include "eigen/Eigen/dense"
#include "MeshFactory.h"
#include "glm/glm.hpp"
//#include "spectra/Spectra/SymEigsSolver.h"
Eigen::MatrixXd generate_L(Mesh * mesh  );
Eigen::MatrixXd generate_W_cotangent_laplacian(Mesh* mesh);
Eigen::MatrixXd generate_A(Mesh* mesh);