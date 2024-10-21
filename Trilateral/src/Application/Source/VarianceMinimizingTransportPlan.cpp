#include "../Include/VarianceMinimizingTransportPlan.h"

static Eigen::MatrixXd generate_cost_function(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2);
static std::vector<float> get_voronoi_areas(TrilateralMesh* m, TrilateralDescriptor& desc1);
Eigen::MatrixXd sinkhornOptimalTransport(const Eigen::MatrixXd& costMatrix, const Eigen::VectorXd& sourceWeights, const Eigen::VectorXd& targetWeights, double epsilon, int maxIter, double tolerance);
static std::vector<float> get_voronoi_areas(TrilateralMesh* m, TrilateralDescriptor& desc1)
{
	std::vector<float> weights; 
	float sum = 0;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		int index = desc1.visited_indices[i];
		float voronoi_area = m->areas[index] / 3.0f;
		weights.push_back(voronoi_area);
		sum += voronoi_area;
	}

	//normalize
	for (size_t i = 0; i < weights.size(); i++)
	{
		weights[i] = weights[i] / sum; 
	}
	return weights;
}

static Eigen::MatrixXd generate_cost_function(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2)
{
	std::vector<float> desc1_weights = get_voronoi_areas(m, desc1);
	std::vector<float> desc2_weights = get_voronoi_areas(m, desc2);
	int size_X = desc1_weights.size();
	int size_Y = desc2_weights.size();

	//generate |X| x |Y| matri
	Eigen::MatrixXd cost(size_X, size_Y);
	for (size_t i = 0; i < size_X; i++)
	{
		for (size_t j = 0; j < size_Y; j++)
		{
			cost(i, j) = std::fabs(desc1_weights[i] - desc2_weights[j]);
		}
	}
	return cost; 
}
// for now include only the INSIDE of descriptors 
float VarianceMin_compare(TrilateralMesh* m,TrilateralDescriptor& desc1, TrilateralDescriptor& desc2)
{
	//create and fill the weight 
	Eigen::VectorXd desc1_weights =stdVectorToEigenVectorXd ( get_voronoi_areas(m,desc1) );
	Eigen::VectorXd desc2_weights =stdVectorToEigenVectorXd( get_voronoi_areas(m,desc2) );
	Eigen::MatrixXd cost_matrix =  generate_cost_function(m, desc1, desc2);
	Eigen::MatrixXd transport_plan =  sinkhornOptimalTransport(cost_matrix, desc1_weights, desc2_weights, 1e-3, 1000, 1e-3);
	float  totalCost = (transport_plan.cwiseProduct(cost_matrix)).sum();
	return totalCost;
}

Eigen::MatrixXd sinkhornOptimalTransport(const Eigen::MatrixXd& costMatrix,
	const Eigen::VectorXd& sourceWeights,
	const Eigen::VectorXd& targetWeights,
	double epsilon, int maxIter, double tolerance) {
	int n = sourceWeights.size();
	int m = targetWeights.size();
	Eigen::MatrixXd K = (-costMatrix / epsilon).array().exp().matrix();
	Eigen::VectorXd u = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd v = Eigen::VectorXd::Ones(m);

	// Iteratively update u and v to find the optimal transport plan
	for (int iter = 0; iter < maxIter; ++iter) {
		Eigen::VectorXd uPrev = u;
		u = sourceWeights.array() / (K * v).array();
		v = targetWeights.array() / (K.transpose() * u).array();

		// Check for convergence
		if ((u - uPrev).norm() < tolerance) {
			break;
		}
	}

	// Compute the final transport plan
	Eigen::MatrixXd transportPlan = u.asDiagonal() * K * v.asDiagonal();;
	return transportPlan;
}


std::vector<std::vector<float>> VarianceMin_compare_all(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& desc_neg)
{
	std::vector<std::vector<float>> comparisons;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		std::vector<float> comparison_i; 
		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			float res = VarianceMin_compare(m, desc_pos[i], desc_neg[j]);
			comparison_i.push_back(res);
		}
		comparisons.push_back(comparison_i);
	}
	return comparisons;
}