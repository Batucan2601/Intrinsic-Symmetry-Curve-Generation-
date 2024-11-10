#pragma once
#include  "../Include/TrilateralMap.h"
#include  "../Include/CoreTypeDefs.h"
#include  "../Include/TrilateralMesh.h"
#include  "../Include/TrilateralDescriptor.h"

//region of interest
std::vector<unsigned int> breadth_first_search(TrilateralMesh* m, int point_index, std::vector<int> is_visited);
void ROI_trilateral(TrilateralMesh* m, TrilateralDescriptor&desc, int division_no, bool is_color);

