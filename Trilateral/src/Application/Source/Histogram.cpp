#include "../Include/Histogram.h"
#include "../Include/Geodesic.h"
#include "../Include/ROI.h"
#define _USE_MATH_DEFINES
#include "math.h"

void Histogram::init(int N)
{
	this->histogram.resize(N, 0);
}
//normalize to N  ( usually 1 ) 
void Histogram::normalize(float N)
{
	//normalize histogram.
	float histogram_sum = 0;
	for (size_t i = 0; i < this->histogram.size(); i++)
	{
		histogram_sum += this->histogram[i];
	}
	for (size_t i = 0; i < this->histogram.size(); i++)
	{
		this->histogram[i] /= histogram_sum;
		this->histogram[i] *= N;
	}

}
int Histogram::size()
{
	return  this->histogram.size();
}
float& Histogram::operator[](int index)
{
	return this->histogram[index];
}

float Histogram_L2Norm(const Histogram& h1, const Histogram& h2)
{
	if (h1.histogram.size() != h1.histogram.size())
	{
		return -1; 
	}
	Eigen::VectorXd  vector_h1 =  stdVectorToEigenVectorXd(h1.histogram);
	Eigen::VectorXd  vector_h2 =  stdVectorToEigenVectorXd(h2.histogram);

	return (vector_h1 - vector_h2).norm();
}

float Histogram_ChiSquareDistance(const Histogram& h1, const Histogram& h2)
{
	if (h1.histogram.size() != h2.histogram.size()) {
		std::cerr << "Error: Vectors must be of the same length." << std::endl;
		return -1;
	}
	float chiSquareDist = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		float numerator = std::pow(h1.histogram[i] - h2.histogram[i], 2);
		float denominator = h1.histogram[i] + h2.histogram[i];
		if (denominator != 0) {
			chiSquareDist += numerator / denominator;
		}
	}

	return chiSquareDist;
}

float Histogram_klDivergence(const Histogram& h1, const Histogram& h2 ) {
	// Check if both vectors are of the same length

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		sumP += h1.histogram[i];
		sumQ += h2.histogram[i];
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return -1;
	}

	double klDiv = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		if (h1.histogram[i] > 0.0 && h2.histogram[i] > 0.0) {
			klDiv += h1.histogram[i] * log(h1.histogram[i] / h2.histogram[i]);
		}
	}

	return klDiv;
}

float Histogram_earthMoversDistance(const Histogram& h1 , const Histogram& h2) {
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1;
	}

	// Check if the vectors are valid probability distributions
	double h1_sum = 0.0, h2_sum = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		h1_sum += h1.histogram[i];
		h2_sum += h2.histogram[i];
	}
	if (fabs(h1_sum - 1.0) > 1e-9 || fabs(h2_sum- 1.0) > 1e-9) {
		return -1;
	}

	// Compute the Earth Mover's Distance
	double cumulativeDifference = 0.0;
	double totalDistance = 0.0;

	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		cumulativeDifference += h1.histogram[i] - h2.histogram[i];
		totalDistance += fabs(cumulativeDifference);
	}

	return totalDistance;
}

float Histogram_bhattacharyyaDistance(const Histogram& h1, const Histogram& h2) 
{
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1; 
	}

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		sumP += h1.histogram[i];
		sumQ += h2.histogram[i];
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return  -1; 
	}

	// Calculate the Bhattacharyya coefficient
	double bc = 0.0;  // Bhattacharyya coefficient
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		if (h1.histogram[i] > 0.0 && h2.histogram[i] > 0.0) {
			bc += sqrt(h1.histogram[i] * h2.histogram[i]);
		}
	}

	// Return the Bhattacharyya distance
	return -log(bc);
}

float Histogram_jensenShannonDivergence(const Histogram& h1 , const Histogram& h2) 
{
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1; 
	}

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		sumP += h1.histogram[i];
		sumQ += h2.histogram[i];
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return -1; 
	}

	// Create the midpoint vector M = (P + Q) / 2
	std::vector<float> M(h1.histogram.size(), 0.0);
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		M[i] = 0.5 * (h1.histogram[i] + h2.histogram[i]);
	}

	// Compute the Jensen-Shannon Divergence
	Histogram M_h = Histogram_Vector_to_Hist(M);
	float jsd = 0.5 * Histogram_klDivergence(h1, M_h) + 0.5 * Histogram_klDivergence(h2, M_h);

	return jsd;
}

float Histogram_kolmogorovSmirnovTest(const Histogram& h1, const Histogram& h2) 
{
	if (h1.histogram.empty() || h2.histogram.empty()) {
		return -1;
	}

	// Combine both samples into a single sorted array
	std::vector<float> combined = h1.histogram;
	combined.insert(combined.end(),h2.histogram.begin(), h2.histogram.end());
	sort(combined.begin(), combined.end());

	// Calculate the CDF for each sample
	std::vector<double> cdf1(h1.histogram.size(), 0.0);
	std::vector<double> cdf2(h2.histogram.size(), 0.0);

	int n1 = h1.histogram.size();
	int n2 = h2.histogram.size();

	// Create step function counts for each sample
	int count1 = 0, count2 = 0;
	double dMax = 0.0;  // Maximum KS statistic

	// Calculate the cumulative step functions and the maximum deviation (KS statistic)
	for (const double& value : combined) {
		while (count1 < n1 && h1.histogram[count1] <= value) {
			count1++;
		}
		while (count2 < n2 && h2.histogram[count2] <= value) {
			count2++;
		}

		// Calculate the CDFs for each sample at the current point
		double cdf1Val = (double)count1 / n1;
		double cdf2Val = (double)count2 / n2;

		// Update the maximum deviation
		dMax = std::max(dMax, fabs(cdf1Val - cdf2Val));
	}

	return dMax;
}









Histogram2D::Histogram2D(int rows, int cols)
{
	this->histogram.resize(rows, cols);
	this->histogram.setZero();
}

void Histogram2D::normalize(float num)
{
	float sum = 0;
	for (size_t i = 0; i < this->histogram.rows(); i++)
	{
		for (size_t j = 0; j < this->histogram.cols(); j++)
		{
			sum += this->histogram(i, j);
		}
	}

	for (size_t i = 0; i < this->histogram.rows(); i++)
	{
		for (size_t j = 0; j < this->histogram.cols(); j++)
		{
			this->histogram(i, j) = this->histogram(i, j)/ sum  * num ;
		}
	}
}

std::pair<int, int> Histogram2D::size()
{
	return std::pair<int, int>(this->histogram.rows()  , this->histogram.cols());
}


float Histogram2D_earthMoversDistance(const Histogram2D& h1, const Histogram2D& h2) 
{
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1;
	}

	// Check if the vectors are valid probability distributions
	double h1_sum = 0.0, h2_sum = 0.0;
	for (size_t i = 0; i < h1.histogram.rows(); ++i) 
	{
		for (size_t j = 0; i < h1.histogram.cols(); ++i)
		{
			h2_sum += h2.histogram(i,j);
			h1_sum += h1.histogram(i,j);
		}
	}
	if (fabs(h1_sum - 1.0) > 1e-9 || fabs(h2_sum - 1.0) > 1e-9) {
		return -1;
	}
	// Compute the Earth Mover's Distance
	double cumulativeDifference = 0.0;
	double totalDistance = 0.0;

	for (size_t i = 0; i < h1.histogram.rows(); ++i) 
	{
		for (size_t j = 0; j < h1.histogram.cols(); j++)
		{
			cumulativeDifference += h1.histogram(i, j) - h2.histogram(i,j);
			totalDistance += fabs(cumulativeDifference);
		}
	}

	return totalDistance;
}

float Histogram2D_bhattacharyyaDistance(const Histogram2D& h1, const Histogram2D& h2)
{
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1;
	}

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.rows(); i++) 
	{
		for (size_t j = 0; j < h1.histogram.cols(); j++)
		{
			sumP += h1.histogram(i,j);
			sumQ += h2.histogram(i,j);
		}
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return  -1;
	}

	// Calculate the Bhattacharyya coefficient
	double bc = 0.0;  // Bhattacharyya coefficient
	for (size_t i = 0; i < h1.histogram.rows(); ++i) {
		for (size_t j = 0; j < h1.histogram.cols(); ++j) {
			if (h1.histogram(i, j) > 0.0 && h2.histogram(i,j)> 0.0) {
				bc += sqrt(h1.histogram(i,j) * h2.histogram(i, j));
			}
		}
	}

	// Return the Bhattacharyya distance
	return -log(bc);
}


float Histogram2D_ChiSquareDistance(const Histogram2D& h1, const Histogram2D& h2)
{
	if (h1.histogram.size() != h2.histogram.size()) {
		std::cerr << "Error: Vectors must be of the same length." << std::endl;
		return -1;
	}
	float chiSquareDist = 0.0;
	for (size_t i = 0; i < h1.histogram.rows(); ++i) 
	{
		for (size_t j = 0; j < h1.histogram.cols(); ++j)
		{
			float numerator = std::pow(h1.histogram(i,j) - h2.histogram(i,j), 2);
			float denominator = h1.histogram(i, j) + h2.histogram(i, j);
			if (denominator != 0) {
				chiSquareDist += numerator / denominator;
			}
		}
		
	}

	return chiSquareDist;
}

float Histogram2D_L2Norm(const Histogram2D& h1, const Histogram2D& h2)
{
	float total_dif = 0;
	for (size_t i = 0; i < h1.histogram.rows(); ++i)
	{
		for (size_t j = 0; j < h1.histogram.cols(); ++j)
		{
			float dif = abs(h1.histogram(i, j) - h2.histogram(i, j));
			total_dif += dif;
		}

	}
	return total_dif;
}


float Histogram_L2Norm_DifferentSize(Histogram& h1, Histogram& h2)
{
	Histogram smaller;
	Histogram bigger;

	int big_size = std::max(h1.histogram.size(), h2.histogram.size());
	if (h1.histogram.size() < h2.histogram.size())
	{
		smaller = h1;
		bigger = h2; 
	}
	else
	{
		smaller = h2;
		bigger = h1;
	}

	// just pad smaller with 0's
	int size_dif = bigger.size() - smaller.size();
	for (size_t i = 0; i < size_dif; i++)
	{
		smaller.histogram.push_back(0.0f);
	}

	return Histogram_L2Norm(smaller, bigger);
}

std::vector<float> Histogram_to_Vector(Histogram& h1)
{
	std::vector<float> vec;
	vec = h1.histogram;
	return vec;
}

Histogram Histogram_Vector_to_Hist(std::vector<float>& h1)
{
	Histogram h; 
	h.init(h1.size());
	h.histogram = h1; 
	return h;
}