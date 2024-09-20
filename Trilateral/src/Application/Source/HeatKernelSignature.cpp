#include "../Include/HeatKernelSignature.h"
#include "eigen/Eigen/EigenValues"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <numeric>

int no_of_samples = 10 ;
int no_of_eigen_val = 3;

void HKS_extract_kernel_signature(Mesh* m)
{
	Eigen::SparseMatrix<double>  cotangent_mat_sparse = cotangent_laplacian(*m);
	cotangent_mat_sparse = normalize_laplacian(cotangent_mat_sparse);

	check_if_matrix_symmetric(cotangent_mat_sparse);

	Spectra::SparseSymMatProd<double> op(cotangent_mat_sparse);

	// Create a SymEigsSolver object using the matrix operation object
	// Compute 3 eigenvalues, with a shift towards the smallest values
	//Spectra::SymEigsSolver< double, Spectra::SmallestAlge, Spectra::SparseSymMatProd<double> > eigs(&op, 3, 6);
	Eigen::Index num_of_eigenvalues = std::min(Eigen::Index(3), cotangent_mat_sparse.rows());
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> solver(op , num_of_eigenvalues,10);
	for (int k = 0; k < cotangent_mat_sparse.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(cotangent_mat_sparse, k); it; ++it) {
			if (std::isnan(it.value()) || std::isinf(it.value())) {
				std::cerr << "Invalid matrix entry at (" << it.row() << ", " << it.col() << "): " << it.value() << std::endl;
				return;
			}
		}
	}
	int max_iterations = 1000;  // Increase maximum number of iterations
	double tolerance = 1e-8;    // Set a smaller tolerance for convergence

	// Initialize and compute the eigenvalues
	solver.init();
	int nconv = solver.compute(Spectra::SortRule::SmallestAlge ,2000);
	Eigen::VectorXd eigenvalues;
	if (solver.info() == Spectra::CompInfo::Successful)
		eigenvalues = solver.eigenvalues();

	float minT = (abs(4 * log(10) / eigenvalues(eigenvalues.size() - 1)));
	float maxT = (abs(4 * log(10) / eigenvalues(1)));

	const double logTMin = log(minT);
	const double logTMax = log(maxT);

	const double tIncrement = (logTMax - logTMin) / (no_of_samples- 1);

	double sumRingAreas = std::accumulate(m->areas.begin(), m->areas.end(), 0.0);
	double sqrtSumRingAreas = sqrt(sumRingAreas);

	const size_t verSize = m->vertices.size();

	std::vector<std::vector<double >> point_descriptor(m->vertices.size());
	for (size_t v = 0; v < verSize; v++)
	{
		if (v < verSize - 1)
		{
			std::cout << "%" << (100 * v) / verSize << " completed for calculating hks descriptor" << "\r";
		}
		else
		{
			std::cout << "%" << 100 << " completed for calculating hks descriptor" << "\n";
		}

		//Initialize descriptor of the vertex
		std::vector<double> descriptor(no_of_samples);

		double currT = logTMin;
		for (size_t t = 0; t < no_of_samples; t++)
		{
			const double logScaleCurrT = exp(currT);
			double heatTrace = 0.0;
			double descVal = 0.0;
			for (size_t eig = 0; eig < no_of_eigen_val; eig++)
			{
				const double expInvEigVal = exp(sumRingAreas * abs(eigenvalues(eig)) * logScaleCurrT * (-1.0));
				const double EigFuncVal = sqrtSumRingAreas * eigenvalues(v, eig);
				const double sqEigFuncVal = EigFuncVal * EigFuncVal;
				descVal += (expInvEigVal * sqEigFuncVal);
				heatTrace += expInvEigVal;
			}
			descriptor[t] = descVal / heatTrace;
			currT += tIncrement;
		}
		point_descriptor[v] = descriptor;
	}


}