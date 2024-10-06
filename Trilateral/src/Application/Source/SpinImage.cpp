#include "../Include/SpinImage.h"


Histogram2D SpinImage_generate_spin_image(TrilateralMesh* m, int reference_point_index,
	std::vector<int>& vertices_in_tri_area, int spin_image_width)
{

	int N = vertices_in_tri_area.size(); //number of vertices in area
	Histogram2D histogram(spin_image_width, spin_image_width);
	// 1 - select a reference point
	glm::vec3 reference_point = m->vertices[reference_point_index];
	glm::vec3 n = m->normals[reference_point_index];

	// 2.1 - find our biggest beta
	float biggest_beta = -INFINITY;
	for (size_t i = 0; i < N; i++)
	{
		glm::vec3 dif = glm::vec3(reference_point- m->vertices[vertices_in_tri_area[i]]);
		float dot = glm::dot(n, dif);
		if (biggest_beta < dot)
		{
			dot = biggest_beta;
		}
	}
	// 2.2 - find our smallest beta
	float smallest_beta = INFINITY;
	for (size_t i = 0; i < N; i++)
	{
		glm::vec3 dif = glm::vec3(reference_point - m->vertices[vertices_in_tri_area[i]]);
		float dot = glm::dot(n, dif);
		if (smallest_beta > dot)
		{
			dot = smallest_beta;
		}
	}
	// 2.3 - find our biggest alpha 
	float biggest_alpha = -INFINITY;
	for (size_t i = 0; i < N; i++)
	{
		glm::vec3 dif = glm::vec3(reference_point - m->vertices[vertices_in_tri_area[i]]);
		float beta = glm::dot(n, dif);
		float len = glm::length(dif);
		float alpha = sqrtf(pow(len, 2) - pow(beta, 2));
		if (alpha > biggest_alpha)
		{
			biggest_alpha = alpha;
		}
	}
	for (size_t i = 0; i < N; i++)
	{
		glm::vec3 p = reference_point; 
		glm::vec3 x = m->vertices[vertices_in_tri_area[i]];

		float beta = glm::dot(n, x - p);
		float len = glm::length(x - p);
		float alpha = sqrtf( pow(len, 2) - pow(beta, 2));

		//normalzie alpha
		int index_a = (alpha / biggest_alpha) * spin_image_width;
		int index_b = ( (beta - smallest_beta) / (biggest_beta - smallest_beta)) * spin_image_width;

		if (index_a == spin_image_width)
		{
			index_a--;
		}
		if (index_b == spin_image_width)
		{
			index_b--;
		}
		histogram.histogram(index_a, index_b) = histogram.histogram(index_a, index_b) + 1 ;
	}
	//normalize 
	histogram.normalize(1);
	return histogram;
}