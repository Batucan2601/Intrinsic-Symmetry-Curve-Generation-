#define _USE_MATH_DEFINES 
#include "../Include/DvorakEstimatingApprox.h"
#include <cmath>
static float gaussian_curvature(Mesh* mesh , int index )
{
	float area = mesh->areas[index];
	float A_i = area / 3; // voronoi

	float gaussian_curvature = 0;
	float phi = 0; 
	float area = 0;
	for (size_t i = 0; i < mesh->adjacenies[index].size(); i++)
	{
		glm::vec3 p1 = mesh->vertices[index]; //itself 
		glm::vec3 p2 = mesh->vertices[mesh->adjacenies[index][i].first]; // second point
		for (size_t j = 0; j < mesh->adjacenies[index].size(); j++)
		{
			for (size_t k = 0; k < mesh->adjacenies[mesh->adjacenies[index][i].first].size(); k++)
			{
				int index_j = mesh->adjacenies[index][j].first;
				int index_k = mesh->adjacenies[mesh->adjacenies[index][i].first][k].first;
				if (index_j == index_k)
				{
					//calculate the cosine
					glm::vec3 p3 = mesh->vertices[index_j];
					float cosine = glm::dot(p3-p1 , p2-p1) / (glm::length(p2 - p1) * glm::length(p3 - p1));
					phi += std::acosf(cosine);
					area += compute_triangle_area(p1, p2, p3) / 3;
				}
			}
		}

	}
	gaussian_curvature = 1.0 / area * (2 * M_PI );
	return gaussian_curvature;
}
Plane dvorak_generate_plane(MeshFactory& mesh_fac, int selected_index)
{
	Mesh* m = &mesh_fac.mesh_vec[selected_index];

}