#include "../Include/FuzzyGeodesic.h"
static float gaussian(float d_x_p , float d_x_q , float d_p_q, float sigma);

static float gaussian(float d_x_p, float d_x_q, float d_p_q , float sigma)
{
	float res = (d_x_p + d_x_q - d_p_q) / sigma;

	res = std::exp(res);

	return res; 
}

FuzzyGeodesicList FuzzyGeodesic_calculateFuzzyGedoesic(Mesh* m, int startIndex, int endIndex , float sigma)
{
	FuzzyGeodesicList fuzzyList; 
	// 1 - calculate dijkstra 
	std::vector<int> path =draw_with_fib_heap_implementation(*m, startIndex, endIndex);

	glm::vec3 p = m->vertices[startIndex];
	glm::vec3 q = m->vertices[endIndex];;
	fuzzyList.startIndex = startIndex;
	fuzzyList.endIndex = endIndex;

	// calcualte total distance
	float total_geodesic_dist = 0;
	std::vector<float> distances(path.size());
	distances[0] = 0;
 	for (size_t i = 1; i < path.size()-1; i++)
	{
		float dist = glm::distance(m->vertices[path[i]], m->vertices[path[i - 1]]);
		total_geodesic_dist += dist; 
		distances[i] += dist + distances[i-1];
	}

	for (size_t i = 0; i < path.size(); i++)
	{
		float d_x_p = distances[i];
		float d_x_q = total_geodesic_dist - distances[i];
		float gaussian_fuzziness = gaussian(d_x_p, d_x_q , total_geodesic_dist, sigma );
		fuzzyList.fuzzyGeodesicList[i].fuzziness = gaussian_fuzziness;
		fuzzyList.fuzzyGeodesicList[i].pointIndex = path[i];
	}

	return fuzzyList;
}

void FuzzyGeodesic_colorMesh(Mesh* m, const FuzzyGeodesicList& fuzzyList)
{

}