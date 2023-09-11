#include "../Include/DominantSymmetry.h"
#include "../Include/glm/glm.hpp"
#include "../../External/Include/eigen/Eigen/Core"
#include "../../External/Include/eigen/Eigen/Eigenvalues" 
using Eigen::MatrixXd;
Plane generate_dominant_symmetry_plane(int seletected_mesh, MeshFactory& mesh_fac) 
{
	Mesh mesh = mesh_fac.mesh_vec[seletected_mesh];
	// generate PCA weights are same and 1 for now 
	float s = mesh.vertices.size();
	glm::vec3 m(0.0f,0.0f,0.0f); 
	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		m += mesh.vertices[i];
	}
	m = m / s; 

	MatrixXd Co(3, 3);

	Co(0, 0) = 0;
	Co(0, 1) = 0;
	Co(0, 2) = 0;
	Co(1, 0) = 0;
	Co(1, 1) = 0;
	Co(1, 2) = 0;
	Co(2, 0) = 0;
	Co(2, 1) = 0;
	Co(2, 2) = 0;
	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		glm::vec3 pi_m;
		pi_m = mesh.vertices[i] - m;
		Eigen::VectorXd pi(3);
		pi(0) = pi_m.x;
		pi(1) = pi_m.y;
		pi(2) = pi_m.z;

		Eigen::MatrixXd  Co_i = pi * pi.transpose();
		Co = Co + Co_i;
	}
	Co = Co / s;

	//// get the eigenvectors 
	Eigen::EigenSolver<MatrixXd> es;
	es.compute(Co/* computeEigenvectors = */ );
	Eigen::MatrixXd eigen_vecs = es.eigenvectors().real();
	Eigen::VectorXd eigen_values = es.eigenvalues().real();

	double biggest_value = -INFINITY;
	int biggest_index = -1; 
	//get the best eigen value
	for (size_t i = 0; i < eigen_values.rows(); i++)
	{
		if (biggest_value < (float)eigen_values(i ))
		{
			biggest_value = (float)eigen_values(i);
			biggest_index = i; 
		}
	}
	std::cout << " eigne values "
		<< std::endl << es.eigenvalues().col(0).real() << std::endl;
	// generate the 3 planes
	Plane planes[3];
	planes[0].normal = glm::vec3(eigen_vecs.col(0).real()(0), eigen_vecs.col(0).real()(1), eigen_vecs.col(0).real()(2));
	planes[0].point = m;
	planes[1].normal = glm::vec3(eigen_vecs.col(1).real()(0), eigen_vecs.col(1).real()(1), eigen_vecs.col(1).real()(2));
	planes[1].point = m;
	planes[2].normal = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	planes[2].point = m;
	
	float minimum_distances[3]; //for each plane 
	for (size_t p = 0; p <3 ; p++) //for all planes 
	{
		float minimum_of_di = INFINITY;
		for (size_t i = 0; i < mesh.vertices.size(); i++)
		{
			// generate sir
			glm::vec3 s_ir = symmetry_point_from_plane( &planes[p] , &mesh.vertices[i] );
			float s_ir_status = get_point_status_from_plane(&planes[p], &s_ir);
			for (size_t j = 0; j < mesh.vertices.size(); j++)
			{
				if (s_ir_status * get_point_status_from_plane(&planes[p], &mesh.vertices[j]) >= 0 )
				{
					float dist = glm::distance(s_ir, mesh.vertices[j]);
					if (dist < minimum_of_di)
					{
						minimum_of_di = dist;
					}
				}
			}
		}
		minimum_distances[p] = minimum_of_di; 
	}
	
	int smallest_dist_index = 0; 
	float smallest_dist = INFINITY; 
	for (size_t i = 0; i < 3; i++)
	{
		if (smallest_dist > minimum_distances[i])
		{
			smallest_dist = minimum_distances[i];
			smallest_dist_index = i; 
		}
	}

	Mesh plane_mesh = generate_mesh_from_plane(&planes[smallest_dist_index], &m);
	mesh_fac.add_mesh(plane_mesh);
	return planes[smallest_dist_index];
	
}
/* We separate the mesh into two in order to get symmetry sets
* m1 represents the points where you get the + sign when you plug the vertices in the plane
* m2 represents the minus sign 
*/
void generate_two_separate_mesh_using_dominant_symmetry_plane(Plane plane, Mesh* mesh_to_be_separated, Mesh* m1, Mesh* m2)
{

	std::vector<int> indices_for_m1(mesh_to_be_separated->vertices.size(), -1); //input old index, gives new index 
	std::vector<int> indices_for_m2(mesh_to_be_separated->vertices.size(), -1); //input old index, gives new index 
	// separate vertices
	for (size_t i = 0; i < mesh_to_be_separated->vertices.size(); i++)
	{
		if (get_point_status_from_plane( &plane , &mesh_to_be_separated->vertices[i]) > 0 )
		{
			indices_for_m1[i] = m1->vertices.size();

			m1->vertices.push_back(mesh_to_be_separated->vertices[i]);
		}
		else
		{
			indices_for_m2[i] = m2->vertices.size();

			m2->vertices.push_back(mesh_to_be_separated->vertices[i]);
		}
	}

	//separate adjacencies
	for (size_t i = 0; i < mesh_to_be_separated->adjacenies.size(); i++)
	{
		float sign_i = get_point_status_from_plane(&plane, &mesh_to_be_separated->vertices[i]);
		std::vector<std::pair<int, float>> adjacencies_m1; 
		std::vector<std::pair<int, float>> adjacencies_m2; 
		for (size_t j = 0; j < mesh_to_be_separated->adjacenies[i].size(); j++)
		{
			float sign_j = get_point_status_from_plane(&plane, &mesh_to_be_separated->vertices[j]);

			if (sign_i  * sign_j > 0 )
			{
				if (sign_i > 0)
				{
					adjacencies_m1.push_back(std::pair <int,float>(indices_for_m1[j], glm::distance(mesh_to_be_separated->vertices[i], mesh_to_be_separated->vertices[j])));
				}
				if (sign_i < 0)
				{
					adjacencies_m2.push_back(std::pair <int, float>(indices_for_m2[j], glm::distance(mesh_to_be_separated->vertices[i], mesh_to_be_separated->vertices[j])));
				}
			}
		}
		m1->adjacenies.push_back(adjacencies_m1);
		m2->adjacenies.push_back(adjacencies_m2);
	}
	//separate triangles
	for (size_t i = 0; i < mesh_to_be_separated->triangles.size(); i += 3 )
	{
		float sign_i = get_point_status_from_plane(&plane, &mesh_to_be_separated->vertices[mesh_to_be_separated->triangles[i]]);
		float sign_i_1 = get_point_status_from_plane(&plane, &mesh_to_be_separated->vertices[mesh_to_be_separated->triangles[i+1]]);
		float sign_i_2 = get_point_status_from_plane(&plane, &mesh_to_be_separated->vertices[mesh_to_be_separated->triangles[i+2]]);

		if (sign_i > 0 && sign_i_1 > 0 && sign_i_2 > 0)
		{
			m1->triangles.push_back(indices_for_m1[i]);
			m1->triangles.push_back(indices_for_m1[i+1]);
			m1->triangles.push_back(indices_for_m1[i+2]);
		}
		else if (sign_i < 0 && sign_i_1 < 0 && sign_i_2 < 0)
		{
			m2->triangles.push_back(indices_for_m2[i]);
			m2->triangles.push_back(indices_for_m2[i + 1]);
			m2->triangles.push_back(indices_for_m2[i + 2]);
		}
	}

}
