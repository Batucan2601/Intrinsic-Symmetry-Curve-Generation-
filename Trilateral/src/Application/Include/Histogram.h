#pragma once
#include  "../Include/TrilateralMesh.h"
#include <eigen/Eigen/Dense>

struct Histogram
{
	std::vector<float> histogram;
	void normalize(float num);
	void init(int N);
	int size();
	float& operator[](int index);
};
struct Histogram2D
{
	Eigen::MatrixXf histogram;
	Histogram2D(int rows, int cols );
	void normalize(float num);
	std::pair<int,int> size();
	float& at( int row , int col);
};
float Histogram_L2Norm(const Histogram& h1 , const Histogram& h2 );
float Histogram_ChiSquareDistance(const Histogram& h1 , const Histogram& h2 );
float Histogram_klDivergence(const Histogram& h1, const Histogram& h2);
float Histogram_jensenShannonDivergence(const Histogram& h1, const Histogram& h2);
float Histogram_bhattacharyyaDistance(const Histogram& h1, const Histogram& h2);
float Histogram_kolmogorovSmirnovTest(const Histogram& h1, const Histogram& h2);

float Histogram2D_earthMoversDistance(const Histogram2D& h1, const Histogram2D& h2);
float Histogram2D_bhattacharyyaDistance(const Histogram2D& h1, const Histogram2D& h2);
float Histogram2D_ChiSquareDistance(const Histogram2D& h1, const Histogram2D& h2);
float Histogram2D_L2Norm(const Histogram2D& h1, const Histogram2D& h2);

float Histogram_L2Norm_DifferentSize(Histogram& h1, Histogram& h2);
std::vector<float> Histogram_to_Vector(Histogram& h1);
Histogram Histogram_Vector_to_Hist(std::vector<float>& h1);