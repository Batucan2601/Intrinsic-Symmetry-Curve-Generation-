#include "CoreTypeDefs.h"
#include "TrilateralDescriptor.h"
float VarianceMin_compare(TrilateralMesh* m, TrilateralDescriptor& desc1, TrilateralDescriptor& desc2);

std::vector<std::vector<float>> VarianceMin_compare_all(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& des_neg);