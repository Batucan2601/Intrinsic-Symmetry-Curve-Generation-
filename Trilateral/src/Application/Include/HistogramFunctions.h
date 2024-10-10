#pragma once
#include "../Include/CoreTypeDefs.h"
#include "../Include/TrilateralDescriptor.h"

Histogram histogram_roi_area_detailed(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited);
Histogram  Histogram_triangle_area(TrilateralMesh* m, TrilateralDescriptor& desc, int division_no, std::vector<int>& global_is_visited);