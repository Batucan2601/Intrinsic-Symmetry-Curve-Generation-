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

// use fuzziness as a direct distance
// for each point
// get the area covered by fuzziness around the point 
void FuzzyGeodesic_FuzzyArea(Mesh* m, const FuzzyGeodesicList& fuzzyList, bool color )
{
	int N = m->vertices.size();
	for (size_t i = 0; i < fuzzyList.fuzzyGeodesicList.size(); i++)
	{
		std::vector<int> isVertexInside(N, false);
		std::vector<float> distances =  compute_geodesic_distances_fibonacci_heap_distances(*m, fuzzyList.fuzzyGeodesicList[i].pointIndex);
		for (size_t j = 0; j < N; j++)
		{
			if (distances[j] > fuzzyList.fuzzyGeodesicList[i].fuzziness)
			{
				isVertexInside[j] = true; 
			}
		}

		float area = 0;
		//get all triangles
		for (size_t i = 0; i < m->triangles.size(); i+=3)
		{
			int index1 = m->triangles[i];
			int index2 = m->triangles[i + 1];
			int index3 = m->triangles[i + 2];
			if (isVertexInside[index1] && isVertexInside[index2] && isVertexInside[index3])
			{
				area = area + compute_triangle_area(m->vertices[index1], m->vertices[index2], m->vertices[index3]);
				if (color)
				{
					m->colors[index1].g = 255; 
				}
			}
		}

	}

}