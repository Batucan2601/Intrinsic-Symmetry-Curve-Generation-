#include "CoreTypeDefs.h"
#include "TrilateralDescriptor.h"
#include "NLateralDescriptor.h"
float VarianceMin_compare(TrilateralMesh* m, TrilateralDescriptor desc1, TrilateralDescriptor desc2, bool is_normalize, int division_no ,int N_ring_no);

std::vector<std::vector<float>> VarianceMin_compare_all(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_pos, std::vector<TrilateralDescriptor>& des_neg,
bool is_normalize , int division_no, int N_ring_no);
std::vector<std::vector<float>> VarianceMin_compare_all_w_CDF(TrilateralMesh* m, std::vector<TrilateralDescriptor>& desc_pos,
	std::vector<TrilateralDescriptor>& desc_neg, int hist_div_no);
std::vector<std::vector<float>> VarianceMin_compare_all(TrilateralMesh* m, std::vector<NLateralDescriptor>& desc_pos, std::vector<NLateralDescriptor>& desc_neg,
	bool is_normalize, int division_no, int N_ring_no);
std::vector<std::vector<float>> VarianceMin_compare_all(TrilateralMesh* m, std::vector<NLateralDescriptor>& descriptors,
	bool is_normalize, int division_no, int N_ring_no);