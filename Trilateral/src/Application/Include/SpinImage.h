#pragma once
#include "Histogram.h"
#include "CoreTypeDefs.h"


Histogram2D SpinImage_generate_spin_image(TrilateralMesh* m, int reference_point_index,
	std::vector<int>& vertices_in_tri_area, int spin_image_width);