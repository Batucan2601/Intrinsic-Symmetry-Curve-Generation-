#include "../Include/HeatKernelSignature.h"
#include "eigen/Eigen/EigenValues"
#include <numeric>

int no_of_samples = 10 ;
int no_of_eigen_val = 5;

void HKS_extract_kernel_signature(Mesh* m)
{
	Eigen::SparseMatrix<double>  cotangent_mat_sparse = cotangent_laplacian(*m);

	Eigen::MatrixXd cotangent_mat = Eigen::MatrixXd(cotangent_mat_sparse);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(cotangent_mat);
	if (es.info() != Eigen::Success) {
		std::cerr << "Eigendecomposition failed!" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Extract the eigenvectors corresponding to the smallest non-zero eigenvalues
	Eigen::VectorXd eigenvalues = es.eigenvalues();
	Eigen::MatrixXd eigenvectors = es.eigenvectors();

	float minT = 1;
	float maxT = 10;
	if (1)
	{
		minT = (abs(4 * log(10) / eigenvalues(eigenvalues.size() - 1)));
		maxT = (abs(4 * log(10) / eigenvalues(1)));
	}

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