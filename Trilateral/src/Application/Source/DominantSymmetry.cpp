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
	glm::vec3 C(0.0f, 0.0f, 0.0f);

	MatrixXd Co(3, 3);
	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		glm::vec3 pi_m;
		pi_m = mesh.vertices[i] - m;
		Co(0, 0) = pi_m.x * pi_m.x;
		Co(0, 1) = pi_m.x * pi_m.y;
		Co(0, 2) = pi_m.x * pi_m.z;
		
		Co(1, 0) = pi_m.y * pi_m.x;
		Co(1, 1) = pi_m.y * pi_m.y;
		Co(1, 2) = pi_m.y * pi_m.z;
		
		Co(2, 0) = pi_m.z * pi_m.x;
		Co(2, 1) = pi_m.z * pi_m.y;
		Co(2, 2) = pi_m.z * pi_m.z;
	}
	//// get the eigenvectors 
	Eigen::EigenSolver<MatrixXd> es;
	es.compute(Co, /* computeEigenvectors = */ true);
	Eigen::MatrixXcd eigen_vecs = es.eigenvectors();
	std::cout << "The first eigenvector of the 3x3 matrix of ones is:"
		<< std::endl << es.eigenvectors().col(0).real()(0)  << std::endl;
	// 3 x 3 3 eigencolumns 
	Plane plane1;
	plane1.normal = glm::vec3(eigen_vecs.col(0).real()(0), eigen_vecs.col(0).real()(1), eigen_vecs.col(0).real()(2));
	plane1.point = m;
	Plane plane2;
	plane2.normal = glm::vec3(eigen_vecs.col(1).real()(0), eigen_vecs.col(1).real()(1), eigen_vecs.col(1).real()(2));
	plane2.point = m;
	Plane plane3;
	plane2.normal = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	plane3.point = m;
	//generate planes using eigen_vecs and point m
	Mesh plane_mesh1 = generate_mesh_from_plane(&plane1, &m);
	Mesh plane_mesh2 = generate_mesh_from_plane(&plane2, &m);
	Mesh plane_mesh3 = generate_mesh_from_plane(&plane3, &m);

	//easy way for now !!! get the plane which cuts y 
	if (std::abs(plane1.normal.z) >= std::abs(plane2.normal.z) && std::abs(plane1.normal.z) >= std::abs(plane3.normal.z) )
	{
		mesh_fac.add_mesh(plane_mesh1);
		return plane1;
	}
	else if (std::abs(plane2.normal.z) >= std::abs(plane3.normal.z) && std::abs(plane2.normal.z) >= std::abs(plane1.normal.z) )
	{
		mesh_fac.add_mesh(plane_mesh2);
		return plane2;
	}
	else if(std::abs(plane3.normal.z) >= std::abs(plane1.normal.z) && std::abs(plane3.normal.z) >= std::abs(plane2.normal.z) )
	{
		mesh_fac.add_mesh(plane_mesh3);
		return plane3;
	}
}