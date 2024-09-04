#include "../Include/DvorakEstimatingApprox.h"

static float gaussian_curvature(Mesh* mesh , int index )
{
	float area = mesh->areas[index];
	float A_i = area / 3; // voronoi

	std::vector<float> phi; 
	for (size_t i = 0; i < mesh->adjacenies[index].size(); i++)
	{
		float phi = 0;
		glm::vec3 p1 = mesh->vertices[index]; //itself 
		glm::vec3 p2 = mesh->vertices[mesh->adjacenies[index][i].first]; // second point
		for (size_t j = 0; j < mesh->adjacenies[index].size(); j++)
		{

		}
	}
}
Plane dvorak_generate_plane(MeshFactory& mesh_fac, int selected_index)
{
	Mesh* m = &mesh_fac.mesh_vec[selected_index];

}