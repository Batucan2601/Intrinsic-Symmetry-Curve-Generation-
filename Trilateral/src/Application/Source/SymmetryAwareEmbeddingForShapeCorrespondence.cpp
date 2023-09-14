#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
#include "../Include/CoreTypeDefs.h"

void generate_isomap_embedding(Mesh* mesh)
{
	int M = 3 ; //dimensions 

	int N = mesh->vertices.size();

	Eigen::MatrixXf delta(N, N);
	// generate geodesic distances matrix ( delta)
	for (size_t i = 0; i < mesh->vertices.size(); i++)
	{
		std::vector<float> distances = compute_geodesic_distances_fibonacci_heap_distances(*mesh, i);
		for (size_t j = 0; j < mesh->vertices.size(); j++)
		{
			delta(i, j) = distances[j];
		}
	}

	//generate J
	Eigen::MatrixXf J(N, N);
	J.fill(-1.0f / (float)N);

	//diagonal
	for (size_t i = 0; i < N; i++)
	{
		J(i, i) = 1 - (1.0f / (float)N);
	}

	//generate K
	Eigen::MatrixXf K(N, N);
	// K = -1/2 J delta^2 J
	K = -1 / 2 * J * (N * N) * J; 

	//eigendecompose isomap

	//solves the eignevector
	Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
	eigensolver.compute(K);
	Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
	Eigen::MatrixXf eigen_vectors = eigensolver.eigenvectors().real();

	// diagonalize and take sqrt
	Eigen::MatrixXf eigen_values_diagonal(M , M );
	eigen_values_diagonal.fill(0);
	for (size_t i = 0; i < M; i++)
	{
		eigen_values_diagonal(i, i) = sqrt(eigen_values(i));
	}

	// take eigen value pairs
	// stated 1-2 1-3 1-4  2-1 3-1 4-1
	//generate 6 plane 
	Plane p_1_2;
	Plane p_1_3;
	Plane p_1_4;

	Plane p_2_1;
	Plane p_3_1;
	Plane p_4_1;


}
