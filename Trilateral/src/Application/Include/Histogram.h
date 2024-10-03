#pragma once
#include  "../Include/TrilateralMesh.h"

struct Histogram
{
	std::vector<float> histogram;
	Histogram(int size);
	Histogram();
	void normalize(float num);
	int size();
	float& operator[](int index);
};


//Histogram Histogram_generate_

Histogram  Histogram_triangle_area(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited);
