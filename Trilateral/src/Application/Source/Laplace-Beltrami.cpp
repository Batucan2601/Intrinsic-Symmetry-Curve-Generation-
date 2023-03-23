#include "../Include/Laplace-Beltrami.h"
using Eigen::MatrixXd;
using Eigen::VectorXcd;

MatrixXd generate_L(Mesh* mesh)
{
	// method =>  L = A_-1 (D - W )
	MatrixXd L(mesh->vertices.size(), mesh->vertices.size());

	// 1 - degree matrix
	MatrixXd D(mesh->vertices.size(), mesh->vertices.size());
	for (size_t i = 0; i < mesh->vertices.size(); i++)
	{
		D(i, i) = mesh->adjacenies[i].size();
	}
	// 2 - Weight matrix
	MatrixXd W = generate_W_cotangent_laplacian(mesh);

	// 3 - diagonal matrix and it's inverse 
	MatrixXd A = generate_A(mesh);
	MatrixXd A_inv = A.inverse();
	L = A_inv * (D - W);
	VectorXcd L_eigenValues = L.eigenvalues();

	std::ostringstream oss;
	std::cout << oss.str() << std::endl; 
	//
}
MatrixXd generate_W_cotangent_laplacian(Mesh* mesh)
{
	MatrixXd W(mesh->vertices.size(), mesh->vertices.size());
	for (size_t i = 0; i < mesh->triangles.size(); i+= 3 )
	{
		//for all of the 3 edges
		glm::vec3 p1 = mesh->vertices[mesh->triangles[i]];
		glm::vec3 p2 = mesh->vertices[mesh->triangles[i+1]];
		glm::vec3 p3 = mesh->vertices[mesh->triangles[i+2]];

		// W_p1_p2
		float cos_p3 = glm::dot(p1 - p3, p2 - p3) / (glm::length(p1 - p3) * glm::length(p2 - p3));
		float sin_p3 = std::sqrt(1 - (cos_p3 * cos_p3));
		float cot_p3 = cos_p3 / sin_p3;
		W(mesh->triangles[i], mesh->triangles[i + 1]) += cot_p3; 
		// W_p1_p3
		float cos_p2 = glm::dot(p1 - p2, p3 - p2) / (glm::length(p1 - p2) * glm::length(p3 - p2));
		float sin_p2 = std::sqrt(1 - (cos_p2 * cos_p2));
		float cot_p2 = cos_p2 / sin_p2;
		W(mesh->triangles[i], mesh->triangles[i + 2]) += cot_p2;
		// W_p2_p3
		float cos_p1 = glm::dot(p2 - p1, p3 - p1) / (glm::length(p2 - p1) * glm::length(p3 - p1));
		float sin_p1 = std::sqrt(1 - (cos_p1 * cos_p1));
		float cot_p1 = cos_p1 / sin_p1;
		W(mesh->triangles[i + 1], mesh->triangles[i + 2]) += cot_p1;

	}
	return W;
}
Eigen::MatrixXd generate_A(Mesh* mesh)
{
	MatrixXd A(mesh->vertices.size(), mesh->vertices.size());
	for (size_t i = 0; i < mesh->triangles.size(); i += 3)
	{
		//for all of the 3 edges
		glm::vec3 p1 = mesh->vertices[mesh->triangles[i]];
		glm::vec3 p2 = mesh->vertices[mesh->triangles[i + 1]];
		glm::vec3 p3 = mesh->vertices[mesh->triangles[i + 2]];

		float area = compute_triangle_area(p1, p2, p3);

		A(mesh->triangles[i]) += area / 3;
		A(mesh->triangles[i + 1]) += area / 3;
		A(mesh->triangles[i + 2]) += area / 3;

	}
	return A; 
}