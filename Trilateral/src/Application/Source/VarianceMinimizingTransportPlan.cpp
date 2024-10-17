#include "../Include/VarianceMinimizingTransportPlan.h"

static std::vector<float> get_voronoi_areas(TrilateralMesh* m, TrilateralDescriptor& desc1)
{
	std::vector<float> weights; 
	float sum = 0;
	for (size_t i = 0; i < desc1.visited_indices.size(); i++)
	{
		int index = desc1.visited_indices[i];
		float voronoi_area = m->areas[index] / 3.0f;
		weights.push_back(voronoi_area);
		sum += voronoi_area;
	}

	//normalize
	for (size_t i = 0; i < weights.size(); i++)
	{
		weights[i] = weights[i] / sum; 
	}
	return weights;
}
// for now include only the INSIDE of descriptors 
void VarianceMin_generate_map(TrilateralMesh* m,TrilateralDescriptor& desc1, TrilateralDescriptor& desc2)
{
	//create and fill the weight 
	std::vector<float> desc1_weights = get_voronoi_areas(m,desc1);
	std::vector<float> desc2_weights = get_voronoi_areas(m,desc2);
	
	{

	}

}