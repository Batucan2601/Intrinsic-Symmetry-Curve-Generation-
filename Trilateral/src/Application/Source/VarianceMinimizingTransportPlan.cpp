#include "../Include/VarianceMinimizingTransportPlan.h"
#include "../Include/Geodesic.h"
#include "../Include/HistogramFunctions.h"
#include "../Include/TrilateralMap.h"

static Eigen::MatrixXd generate_cost_function(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2);
static std::vector<float> get_voronoi_areas(TrilateralMesh* m, TrilateralDescriptor& desc1);
Eigen::MatrixXd sinkhornOptimalTransport(const Eigen::MatrixXd& costMatrix, const Eigen::VectorXd& sourceWeights, const Eigen::VectorXd& targetWeights, double epsilon, int maxIter, double tolerance);
float VarianceMin_compare_w_CDF(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2, Histogram& cdf_desc1, Histogram& cdf_desc2);
static Eigen::MatrixXd generate_cost_function_w_PDF_CDF(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2,
	Histogram& h1, Histogram& h2);

static std::vector<float> get_voronoi_areas(TrilateralMesh* m, TrilateralDescriptor& desc1);
static std::vector<float> get_voronoi_areas_w_paths(TrilateralMesh* m, TrilateralDescriptor& desc1);
static std::vector<float> get_voronoi_areas_divided(TrilateralMesh* m, TrilateralDescriptor& desc1 , int division_no);
static std::vector<float> get_n_ring_areas_divided(TrilateralMesh* m, TrilateralDescriptor& desc1, int division_no, int N_ring);

static std::vector<float> get_n_ring_areas_divided(TrilateralMesh* m, TrilateralDescriptor& desc1, int division_no, int N_ring)
{
	std::vector<float>  distances_p1 = Geodesic_dijkstra(*m, desc1.p1);
	std::vector<float> voronoi_areas_divided(division_no, 0);
	float longest = -INFINITY;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		float dist = distances_p1[desc1.visited_indices[i]];
		if (dist > longest)
		{
			longest = dist;
		}
	}
	for (size_t i = 0; i < desc1.path_1_2.size(); i++)
	{
		float dist = distances_p1[desc1.path_1_2[i]];
		if (dist > longest)
		{
			longest = dist;
		}
	}
	for (size_t i = 0; i < desc1.path_1_3.size(); i++)
	{
		float dist = distances_p1[desc1.path_1_3[i]];
		if (dist > longest)
		{
			longest = dist;
		}
	}
	for (size_t i = 0; i < desc1.path_2_3.size(); i++)
	{
		float dist = distances_p1[desc1.path_2_3[i]];
		if (dist > longest)
		{
			longest = dist;
		}
	}

	float step = longest / division_no;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		int index = desc1.visited_indices[i];
		float dist = distances_p1[index];
		int hist_index = dist / step;
		float n_ring_area= get_N_ring_area(m, index, N_ring) ;
		if (hist_index == division_no)
		{
			hist_index--;
		}

		voronoi_areas_divided[hist_index] += n_ring_area;// *m->normalized_heat_kernel_signature[index];
	}

	//get paths point areas
	for (size_t i = 0; i < desc1.path_1_2.size(); i++)
	{
		int index = desc1.path_1_2[i];
		float dist = distances_p1[index];
		int hist_index = dist / step;
		float area = m->areas[index];
		if (hist_index == division_no)
		{
			hist_index--;
		}
		voronoi_areas_divided[hist_index] += area / 3.0;
	}
	//get paths point areas
	for (size_t i = 0; i < desc1.path_1_3.size(); i++)
	{
		int index = desc1.path_1_3[i];
		float dist = distances_p1[index];
		int hist_index = dist / step;
		float area = m->areas[index];
		if (hist_index == division_no)
		{
			hist_index--;
		}
		voronoi_areas_divided[hist_index] += area / 3.0;// *m->normalized_heat_kernel_signature[index];
	}
	//get paths point areas
	for (size_t i = 0; i < desc1.path_2_3.size(); i++)
	{
		int index = desc1.path_2_3[i];
		float dist = distances_p1[index];
		int hist_index = dist / step;
		float area = m->areas[index];
		if (hist_index == division_no)
		{
			hist_index--;
		}
		voronoi_areas_divided[hist_index] += area / 3.0;// *m->normalized_heat_kernel_signature[index];
	}
	return voronoi_areas_divided;
}
static std::vector<float> get_voronoi_areas_divided(TrilateralMesh* m, TrilateralDescriptor& desc1, int division_no)
{
	std::vector<float>  distances_p1 = Geodesic_dijkstra(*m, desc1.p1);
	std::vector<float> voronoi_areas_divided(division_no, 0);
	float longest = -INFINITY;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		float dist = distances_p1[desc1.visited_indices[i]];
		if (dist > longest)
		{
			longest = dist; 
		}
	}
	float step = longest / division_no;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		int index = desc1.visited_indices[i];
		float dist = distances_p1[index];
		int hist_index = dist / step;
		float voronoi_area = m->areas[index] / 3.0f;
		if (hist_index == division_no)
		{
			hist_index--;
		}

		voronoi_areas_divided[hist_index] += voronoi_area;
	}
	return voronoi_areas_divided;
}
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
	/*for (size_t i = 0; i < weights.size(); i++)
	{
		weights[i] = weights[i] / sum; 
	}*/
	return weights;
}
static std::vector<float> get_voronoi_areas_w_paths(TrilateralMesh* m, TrilateralDescriptor& desc1)
{

	std::vector<float> weights;
	float sum = 0;
	for (size_t i = 0; i < desc1.path_1_2.size(); i++)
	{
		int index = desc1.path_1_2[i];
		float voronoi_area = m->areas[index] / 3.0f;
		weights.push_back(voronoi_area);
		sum += voronoi_area;
	}
	for (size_t i = 0; i < desc1.path_1_3.size(); i++)
	{
		int index = desc1.path_1_3[i];
		float voronoi_area = m->areas[index] / 3.0f;
		weights.push_back(voronoi_area);
		sum += voronoi_area;
	}
	for (size_t i = 0; i < desc1.path_2_3.size(); i++)
	{
		int index = desc1.path_2_3[i];
		float voronoi_area = m->areas[index] / 3.0f;
		weights.push_back(voronoi_area);
		sum += voronoi_area;
	}
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		int index = desc1.visited_indices[i];
		float voronoi_area = m->areas[index] / 3.0f;
		weights.push_back(voronoi_area);
		sum += voronoi_area;
	}

	//normalize
	/*for (size_t i = 0; i < weights.size(); i++)
	{
		weights[i] = weights[i] / sum;
	}*/
	return weights;
}

static Eigen::MatrixXd generate_cost_function(TrilateralMesh* m, TrilateralDescriptor&desc1, TrilateralDescriptor& desc2)
{
	//std::vector<float> desc1_weights = get_voronoi_areas_w_paths(m, desc1);
	//std::vector<float> desc2_weights = get_voronoi_areas_w_paths(m, desc2);
	Eigen::VectorXd desc1_weights = desc1.weight;
	Eigen::VectorXd desc2_weights = desc2.weight;

	int size_X = desc1_weights.size();
	int size_Y = desc2_weights.size();

	float desc1_hks = m->normalized_heat_kernel_signature[desc1.p1];
	float desc2_hks = m->normalized_heat_kernel_signature[desc2.p1];
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

static Eigen::MatrixXd generate_cost_function_w_PDF_CDF(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2,
	Histogram& h1, Histogram& h2)
{
	std::vector<float> desc1_weights = get_voronoi_areas(m, desc1);
	std::vector<float> desc2_weights = get_voronoi_areas(m, desc2);
	int size_X = desc1_weights.size();
	int size_Y = desc2_weights.size();

	// 1 - we need to know each point's distance from p1
	std::vector<float> dist_array_1 = Geodesic_dijkstra(*m, desc1.p1);
	float desc_1_max_size = -INFINITY;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		int index = desc1.visited_indices[i];
		float dist = dist_array_1[index];
		if (dist > desc_1_max_size)
		{
			desc_1_max_size = dist;
		}
	}
	std::vector<float> dist_array_2 = Geodesic_dijkstra(*m, desc2.p1);
	float desc_2_max_size = -INFINITY;
	for (size_t i = 0; i < desc2.visited_indices.size(); i++)
	{
		int index = desc2.visited_indices[i];
		float dist = dist_array_2[index];
		if (dist > desc_2_max_size)
		{
			desc_2_max_size = dist;
		}
	}
	//generate |X| x |Y| matri
	Eigen::MatrixXd cost(size_X, size_Y);
	for (size_t i = 0; i < size_X; i++)
	{
		int hist_no_i = dist_array_1[desc1.visited_indices[i]] / desc_1_max_size * h1.size();
		if (h1.size() == hist_no_i)
		{
			hist_no_i = h1.size() - 1;
		}
		float hist_i = h1[hist_no_i];
		for (size_t j = 0; j < size_Y; j++)
		{
			int hist_no_j = dist_array_2[desc2.visited_indices[j]] / desc_2_max_size * h2.size();
			if (h2.size() == hist_no_j)
			{
				hist_no_j = h2.size() - 1;
			}
			float hist_j = h2[hist_no_j];
			cost(i, j) = std::fabs(desc1_weights[i] - desc2_weights[j]) * std::fabs(1 -  (hist_i - hist_j));
			if (std::isnan(cost(i, j)))
			{
				std::cout << " nan " << std::endl; 
			}
		}
	}

	//lets normalize
	//double minCoeff = cost.minCoeff();
	//double maxCoeff = cost.maxCoeff();
	//cost = (cost.array() - minCoeff) / (maxCoeff - minCoeff);


	return cost;
}
// for now include only the INSIDE of descriptors 
float VarianceMin_compare(TrilateralMesh* m,TrilateralDescriptor desc1, TrilateralDescriptor desc2, bool is_normalize, int division_no ,int N_ring_no )
{
	
	if (is_normalize)
	{
		desc1.weight = desc1.weight  / desc1.weight.sum();
		desc2.weight = desc2.weight  / desc2.weight.sum();
	}
	Eigen::MatrixXd cost_matrix =  generate_cost_function(m, desc1, desc2);
	Eigen::MatrixXd transport_plan =  sinkhornOptimalTransport(cost_matrix, desc1.weight, desc2.weight, 1e-1, 10000, 1e-8);
	float  totalCost = (transport_plan.cwiseProduct(cost_matrix)).sum();
	
	std::cout << "====================================================================" << std::endl;
	std::cout << "transport plan " << transport_plan << std::endl;
	std::cout << "====================================================================" << std::endl;
	std::cout << " cost matrix " << cost_matrix << std::endl;
	std::cout << "====================================================================" << std::endl;
	std::cout << "weight 1 " << desc1.weight << std::endl;
	std::cout << "====================================================================" << std::endl;
	std::cout << "weight 2 " << desc2.weight << std::endl;
	std::cout << "====================================================================" << std::endl;
	std::cout << "total cost" << totalCost  << std::endl;
	std::cout << "====================================================================" << std::endl;
	return totalCost;
}

double logSumExp(const Eigen::VectorXd& vec) {

	double maxVal = vec.maxCoeff();
	return maxVal + std::log((vec.array() - maxVal).exp().sum());
}
double kl_divergence(const Eigen::VectorXd& P, const Eigen::VectorXd& Q) {
	// Check that P and Q have the same size
	assert(P.size() == Q.size() && "Vectors P and Q must be the same size");

	// Check that both vectors are probability distributions (sum to 1)
	assert(std::abs(P.sum() - 1.0) < 1e-9 && "Vector P must sum to 1");
	assert(std::abs(Q.sum() - 1.0) < 1e-9 && "Vector Q must sum to 1");

	// Compute KL divergence
	double kl_div = 0.0;
	for (int i = 0; i < P.size(); ++i) {
		if (P(i) > 0 && Q(i) > 0) {  // Avoid division by zero and log(0)
			kl_div += P(i) * std::log(P(i) / Q(i));
		}
	}
	return kl_div;
}
// Sinkhorn algorithm to solve regularized optimal transport problem
Eigen::MatrixXd sinkhornOptimalTransport(const Eigen::MatrixXd& costMatrix,
	const Eigen::VectorXd& sourceWeights,
	const Eigen::VectorXd& targetWeights,
	double epsilon,  int maxIter, double tolerance) {
	int X = sourceWeights.size();
	int Y = targetWeights.size();
	Eigen::MatrixXd H = (-costMatrix / epsilon).array().exp().matrix();
	Eigen::VectorXd v = Eigen::VectorXd::Ones(X);
	Eigen::VectorXd w = Eigen::VectorXd::Ones(Y);
	Eigen::VectorXd v_prev = v;
	double kl_prev = 0;
	std::cout << "new " << std::endl;
	for (int iter = 0; iter < maxIter; ++iter) {

		v =  H * (targetWeights.array() * w.array()).matrix();
		v = Eigen::VectorXd::Ones(X).array() * v.array() ;

		w = H.transpose() * (sourceWeights.array() * v.array()).matrix();
		w = Eigen::VectorXd::Ones(Y).array() * w.array();

		double kl_cur = kl_divergence(v / v.sum(), v_prev / v_prev.sum());
		if (iter != 0)
		{
			double diff = std::abs(kl_cur - kl_prev);
			std::cout << "divergence " << kl_cur <<  " difference in divergence " <<  diff << std::endl;
			
			if (diff < tolerance)
			{
				break;
			}
		}
		v_prev = v;
		kl_prev = kl_cur; 
	}

	Eigen::MatrixXd transportPlan = sourceWeights.asDiagonal().toDenseMatrix() * v.asDiagonal().toDenseMatrix() * H
	* w.asDiagonal().toDenseMatrix() * targetWeights.asDiagonal().toDenseMatrix();
	return transportPlan;
}

// for now include only the INSIDE of descriptors 
float VarianceMin_compare_w_CDF(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2 , Histogram& cdf_desc1 , Histogram& cdf_desc2)
{
	//create and fill the weight 
	Eigen::VectorXd desc1_weights = stdVectorToEigenVectorXd(get_voronoi_areas(m, desc1));
	if (desc1_weights.array().isNaN().any())
	{
		std::cout << "The matrix contains NaN values." << std::endl;
	}
	Eigen::VectorXd desc2_weights = stdVectorToEigenVectorXd(get_voronoi_areas(m, desc2));
	if (desc2_weights.array().isNaN().any()) 
	{
		std::cout << "The matrix contains NaN values." << std::endl;
	}

	Eigen::MatrixXd cost_matrix = generate_cost_function_w_PDF_CDF(m, desc1, desc2 , cdf_desc1 , cdf_desc2 );
	//Eigen::MatrixXd cost_matrix = generate_cost_function(m, desc1, desc2  );
	if (cost_matrix.array().isNaN().any()) {
		std::cout << "The matrix contains NaN values." << std::endl;
	}
	Eigen::MatrixXd transport_plan = sinkhornOptimalTransport(cost_matrix, desc1_weights, desc2_weights, 1e-3, 1000, 1e-3);
	if (transport_plan.array().isNaN().any()) {
		std::cout << "The matrix contains NaN values." << std::endl;
	}
	float  totalCost = (transport_plan.cwiseProduct(cost_matrix)).sum();
	return totalCost;
}
// for now include only the INSIDE of descriptors 
float VarianceMin_compare_w_PDF(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2, Histogram& cdf_desc1, Histogram& cdf_desc2)
{
	//create and fill the weight 
	Eigen::VectorXd desc1_weights = stdVectorToEigenVectorXd(get_voronoi_areas(m, desc1));
	Eigen::VectorXd desc2_weights = stdVectorToEigenVectorXd(get_voronoi_areas(m, desc2));
	Eigen::MatrixXd cost_matrix = generate_cost_function_w_PDF_CDF(m, desc1, desc2, cdf_desc1, cdf_desc2);
	Eigen::MatrixXd transport_plan = sinkhornOptimalTransport(cost_matrix, desc1_weights, desc2_weights, 5 * 1e-4, 300, 5 * 1e-4);
	float  totalCost = (transport_plan.cwiseProduct(cost_matrix)).sum();
	return totalCost;
}


std::vector<std::vector<float>> VarianceMin_compare_all(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& desc_neg,
bool is_normalize, int division_no ,int N_ring_no )
{
	std::vector<std::vector<float>> comparisons;
	for (size_t i = 0; i < desc_pos.size() ; i++)
	{
		Eigen::VectorXd desc1_weights = stdVectorToEigenVectorXd(get_n_ring_areas_divided(m, desc_pos[i], division_no, N_ring_no));
		desc_pos[i].weight = desc1_weights;
	}
	for (size_t i = 0; i < desc_neg.size(); i++)
	{
		Eigen::VectorXd desc1_weights = stdVectorToEigenVectorXd(get_n_ring_areas_divided(m, desc_neg[i], division_no, N_ring_no));
		desc_neg[i].weight = desc1_weights;
	}
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		std::vector<float> comparison_i; 
 		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			std::cout << " source = " << i << " target == " << desc_pos.size() + j << std::endl; 
			float res = VarianceMin_compare(m, desc_pos[i], desc_neg[j], is_normalize , division_no, N_ring_no);
			comparison_i.push_back(res);
		}
		comparisons.push_back(comparison_i);
	}

	return comparisons;
}
std::vector<std::vector<float>> VarianceMin_compare_all_w_CDF(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& desc_neg
, int hist_div_no)
{
	std::vector<Histogram> cdf_pos;
	std::vector<Histogram> cdf_neg;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		Histogram cdf = TrilateralDescriptor_generate_cdf_of_areas(m, desc_pos[i], hist_div_no);
		cdf_pos.push_back(cdf);
	}
	for (size_t i = 0; i < desc_neg.size(); i++)
	{
		Histogram cdf = TrilateralDescriptor_generate_cdf_of_areas(m, desc_neg[i], hist_div_no);
		cdf_neg.push_back(cdf);
	}

	std::vector<std::vector<float>> comparisons;
	for (size_t i = 0; i < desc_pos.size(); i++)
	{
		std::vector<float> comparison_i;
		for (size_t j = 0; j < desc_neg.size(); j++)
		{
			float res = VarianceMin_compare_w_CDF(m, desc_pos[i], desc_neg[j] , cdf_pos[i] , cdf_neg[j]);
			comparison_i.push_back(res);
		}
		comparisons.push_back(comparison_i);
	}
	return comparisons;
}
