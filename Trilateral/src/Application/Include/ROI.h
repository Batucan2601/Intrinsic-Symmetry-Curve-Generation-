#pragma once
#include  "../Include/TrilateralMap.h"
#include  "../Include/CoreTypeDefs.h"
#include  "../Include/TrilateralMesh.h"

//region of interest
std::vector<int> ROI_trilateral(TrilateralMesh* m, TrilateralDescriptor&desc, int division_no, bool is_color);

