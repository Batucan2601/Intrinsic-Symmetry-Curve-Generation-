#pragma once
#include  "../Include/TrilateralMap.h"
#include  "../Include/CoreTypeDefs.h"
#include  "../Include/TrilateralMesh.h"

//region of interest
std::vector<int> ROI_trilateral(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no, bool is_color);

