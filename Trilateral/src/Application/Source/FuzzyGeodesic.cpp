#include "../Include/FuzzyGeodesic.h"
#include "../Include/CoreTypeDefs.h"
#include "../Include/Geodesic.h"


static float fuzzy_gaussian(float d_x_p, float d_x_q, float d_p_q , float sigma)
{
	float res = (d_x_p + d_x_q - d_p_q) / sigma;

	res = std::exp(res);

	return res; 
}

FuzzyGeodesicList FuzzyGeodesic_calculateFuzzyGedoesic(TrilateralMesh* m, int startIndex, int endIndex , float fuzziness_sigma)
{
	FuzzyGeodesicList fuzzyList; 
	// 1 - calculate dijkstra 
	std::vector<int> path =Geodesic_between_two_points(*m, startIndex, endIndex);

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
		FuzzyGeodesic f_geo;
		std::vector<float> distances_for_x  = Geodesic_dijkstra(*m ,path[i]);
		float d_x_p = distances_for_x[startIndex];
		float d_x_q = distances_for_x[endIndex];
		f_geo.fuzziness = fuzzy_gaussian(d_x_p, d_x_q , total_geodesic_dist, fuzziness_sigma);
		f_geo.pointIndex = path[i];
		fuzzyList.fuzzyGeodesicList.push_back(f_geo);
	}

	return fuzzyList;
}

// use fuzziness as a direct distance
// for each point
// get the area covered by fuzziness around the point 
float FuzzyGeodesic_FuzzyArea(TrilateralMesh* m, const FuzzyGeodesicList& fuzzyList, bool color )
{
	int N = m->vertices.size();
	std::vector<int> isPath(N, false);
	


	float area = 0;

	for (size_t i = 0; i < fuzzyList.fuzzyGeodesicList.size(); i++)
	{
		std::vector<int> isVertexInside(N, false);
		std::vector<float> distances =  Geodesic_dijkstra(*m, fuzzyList.fuzzyGeodesicList[i].pointIndex);
		for (size_t j = 0; j < N; j++)
		{
			if (distances[j] < fuzzyList.fuzzyGeodesicList[i].fuzziness)
			{
				isVertexInside[j] = true; 
			}
		}

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
					m->colors[index1].r = 0;
					m->colors[index1].g = 255; 
					m->colors[index1].b = 0; 
				}

			}
		}

	}

	//for the path itself 
	if (color)
	{
		for (size_t i = 0; i < fuzzyList.fuzzyGeodesicList.size(); i++)
		{
			m->colors[fuzzyList.fuzzyGeodesicList[i].pointIndex].r = 255; 
			m->colors[fuzzyList.fuzzyGeodesicList[i].pointIndex].g = 0; 
			m->colors[fuzzyList.fuzzyGeodesicList[i].pointIndex].b = 0; 
		}
	}

	return area;
}