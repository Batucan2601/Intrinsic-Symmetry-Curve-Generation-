#pragma once
#include  "../Include/TrilateralMesh.h"

struct Histogram
{
	std::vector<float> histogram;
	Histogram(int size);
	Histogram(std::vector<float> vec);
	Histogram();
	void normalize(float num);
	int size();
	float& operator[](int index);
};

float Histogram_L2Norm(const Histogram& h1 , const Histogram& h2 );
float Histogram_ChiSquareDistance(const Histogram& h1 , const Histogram& h2 );
float Histogram_klDivergence(const Histogram& h1, const Histogram& h2);
float Histogram_jensenShannonDivergence(const Histogram& h1, const Histogram& h2);
float Histogram_bhattacharyyaDistance(const Histogram& h1, const Histogram& h2);
float Histogram_kolmogorovSmirnovTest(const Histogram& h1, const Histogram& h2);

Histogram  Histogram_triangle_area(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited);
