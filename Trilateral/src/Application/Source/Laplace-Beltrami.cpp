#include "../Include/Laplace-Beltrami.h"
#include "eigen/Eigen/EigenValues"
#include "eigen/Eigen/Sparse"
#include <vector>
#include "../Include/NLateralDescriptor.h"
#include "../Include/CoreTypeDefs.h"
#include <windows.h>
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::Vector3d;

// Function to compute the cotangent of an angle
static double cotangent(const glm::vec3& v1, const glm::vec3& v2) {
	return glm::dot(v1,v2) / glm::length(glm::cross(v1,v2));
}

// Function to compute the area of the Voronoi region around a vertex
static double voronoi_area(const Mesh& mesh, int vertex) {
	double area = 0.0;
	for (int i = 0; i < mesh.triangles.size(); i += 3 ) {
		if ( mesh.triangles[i] == vertex || mesh.triangles[i + 1] == vertex || mesh.triangles[i + 2] == vertex) {
			glm::vec3 v0 = mesh.vertices[mesh.triangles[i]];
			glm::vec3 v1 = mesh.vertices[mesh.triangles[i + 1]];
			glm::vec3 v2 = mesh.vertices[mesh.triangles[i + 2]];
			double a = glm::length((v1 - v0));
			double b = glm::length((v2 - v1));
			double c = glm::length((v0 - v2));
			double s = (a + b + c) / 2.0;
			area += sqrt(s * (s - a) * (s - b) * (s - c)) / 3.0;
		}
	}
	return area;
}

// Function to assemble the cotangent Laplacian matrix
Eigen::SparseMatrix<double> cotangent_laplacian(const Mesh& mesh)
{
	int n = mesh.vertices.size();
	Eigen::SparseMatrix<double> L(n, n);
	std::vector<Eigen::Triplet<double>> triplets;

	for (int index = 0; index < mesh.triangles.size(); index += 3) {
		int i = mesh.triangles[index];
		int j = mesh.triangles[index + 1];
		int k = mesh.triangles[index + 2];

		glm::vec3 vi = mesh.vertices[i];
		glm::vec3 vj = mesh.vertices[j];
		glm::vec3 vk = mesh.vertices[k];

		double cot_alpha = cotangent(vj - vi, vk - vi);
		double cot_beta = cotangent(vi - vj, vk - vj);
		double cot_gamma = cotangent(vi - vk, vj - vk);

		// Symmetrically add off-diagonal entries
		triplets.emplace_back(i, j, -cot_alpha / 2.0 );
		triplets.emplace_back(j, i, -cot_alpha / 2.0 );

		triplets.emplace_back(j, k, -cot_beta / 2.0);
		triplets.emplace_back(k, j, -cot_beta / 2.0);

		triplets.emplace_back(k, i, -cot_gamma / 2.0);
		triplets.emplace_back(i, k, -cot_gamma / 2.0);

		// Diagonal entries (sum of cotangent weights for the current vertex)
		triplets.emplace_back(i, i, (cot_alpha + cot_gamma) / 2.0);
		triplets.emplace_back(j, j, (cot_alpha + cot_beta) / 2.0);
		triplets.emplace_back(k, k, (cot_beta + cot_gamma) / 2.0);
	}

	L.setFromTriplets(triplets.begin(), triplets.end());
	L.makeCompressed();
	return L;
}

Eigen::SparseMatrix<double> normalize_laplacian(const Eigen::SparseMatrix<double>& L)
{
	Eigen::VectorXd D = L.diagonal();
	Eigen::SparseMatrix<double> D_inv_sqrt(L.rows(), L.cols());

	std::vector<Eigen::Triplet<double>> triplets;
	for (int i = 0; i < D.size(); ++i) {
		if (D(i) > 1e-9) {  // Avoid division by zero for isolated vertices
			triplets.emplace_back(i, i, 1.0 / std::sqrt(D(i)));
		}
	}

	D_inv_sqrt.setFromTriplets(triplets.begin(), triplets.end());

	// L_normalized = D_inv_sqrt * L * D_inv_sqrt
	return D_inv_sqrt * L * D_inv_sqrt;
}
Eigen::SparseMatrix<double> regularize_matrix(const Eigen::SparseMatrix<double>& L, double epsilon)
{
	Eigen::SparseMatrix<double> I(L.rows(), L.cols());
	I.setIdentity();
	return L + epsilon * I;
}

bool check_if_matrix_symmetric(const Eigen::SparseMatrix<double>& mat)
{
	if (mat.rows() != mat.cols()) {
		std::cerr << "Matrix is not square!" << std::endl;
		return false;
	}

	for (int k = 0; k < mat.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
			if (std::abs(it.value() - mat.coeff(it.col(), it.row())) > 1e-9) {
				std::cerr << "Matrix is not symmetric at (" << it.row() << ", " << it.col() << ")" << std::endl;
				return false;
			}
		}
	}
	std::cout << "Matrix is symmetric." << std::endl;
	return true;
}
// Function to compute the Laplace-Beltrami operator
static MatrixXd laplace_beltrami(const Mesh& mesh) {
	MatrixXd L = cotangent_laplacian(mesh);
	int n = mesh.vertices.size();
	MatrixXd A = MatrixXd::Zero(n, n);

	for (int i = 0; i < n; ++i) {
		A(i, i) = voronoi_area(mesh, i);
	}

	return A.inverse() * L;
}

// Function to embed the original vertices to a 2D domain
MatrixXd embed_mesh_to_2d(Mesh& mesh) {

	// Compute Laplace-Beltrami operator
	MatrixXd L = laplace_beltrami(mesh);

	// Compute eigendecomposition
	Eigen::SelfAdjointEigenSolver<MatrixXd> es(L);
	if (es.info() != Eigen::Success) {
		std::cerr << "Eigendecomposition failed!" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Extract the eigenvectors corresponding to the smallest non-zero eigenvalues
	Eigen::VectorXd eigenvalues = es.eigenvalues();
	MatrixXd eigenvectors = es.eigenvectors();

	// Skip the first eigenvector (corresponding to eigenvalue 0)
	// The second and third smallest eigenvalues correspond to the embedding
	MatrixXd embedding(mesh.vertices.size(), 2);
	embedding.col(0) = eigenvectors.col(1); // Second smallest eigenvalue
	embedding.col(1) = eigenvectors.col(2); // Third smallest eigenvalue


	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		mesh.vertices[i].z = 0.0f;
		mesh.vertices[i].x = embedding(0,i);
		mesh.vertices[i].y = embedding(1,i);
	}
	
	return embedding;
}

Eigen::MatrixXd embed_mesh_endpoints_to_2d(Mesh& mesh, Skeleton& skeleton, NLateralParameters& nLateralParameters )
{
	std::vector<unsigned int> mesh_indices; //for skeleton end points
	std::vector<NLateralDescriptor> n_lateral_list;
	std::map<int,int> pair_index_to_mesh_index;  
	int size_of_endpoints = 0;
	skeleton_calculate_closest_mesh_points(skeleton ,&mesh , mesh_indices);
	n_lateral_list = get_N_lateral_descriptor_using_closest_pairs(&mesh, mesh_indices, nLateralParameters);
	Mesh mesh_endpoints = mesh;
	mesh_endpoints.colors.clear();
	mesh_endpoints.vertices.clear();
	mesh_endpoints.triangles.clear();

	size_of_endpoints = mesh_indices.size();

	//make it into a mesh
	for (size_t i = 0; i < size_of_endpoints; i++)
	{
		mesh_endpoints.vertices.push_back(mesh.vertices[mesh_indices[i]]);
		mesh_endpoints.colors.push_back(glm::vec3(0.0f , 255.0f , 0.0f));

		pair_index_to_mesh_index[mesh_indices[i]] =  i ;
	}
	//generate Degree matrix
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(size_of_endpoints, size_of_endpoints);
	for (size_t i = 0; i < size_of_endpoints; i++)
	{
		D(i, i) = n_lateral_list[i].point_indices.size() - 1;
	}
	//generate adjacency matrix 
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size_of_endpoints, size_of_endpoints);
	for (size_t i = 0; i < size_of_endpoints; i++)
	{
		for (size_t j = 1; j < n_lateral_list[i].point_indices.size(); j++)
		{
			A(i, pair_index_to_mesh_index[n_lateral_list[i].point_indices[j]]) = 1;
		}
	}
	//Laplacian matrix
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(size_of_endpoints, size_of_endpoints);
	L = D - A;
	
	// Compute eigendecomposition
	Eigen::SelfAdjointEigenSolver<MatrixXd> es(L);
	if (es.info() != Eigen::Success) {
		std::cerr << "Eigendecomposition failed!" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Extract the eigenvectors corresponding to the smallest non-zero eigenvalues
	Eigen::VectorXd eigenvalues = es.eigenvalues();
	MatrixXd eigenvectors = es.eigenvectors();

	// Skip the first eigenvector (corresponding to eigenvalue 0)
	// The second and third smallest eigenvalues correspond to the embedding
	MatrixXd embedding(size_of_endpoints, 3);
	embedding.col(0) = eigenvectors.col(1); // Second smallest eigenvalue
	embedding.col(1) = eigenvectors.col(2); // Third smallest eigenvalue

	for (size_t i = 0; i < n_lateral_list.size(); i++)
	{
		for (size_t j = 1; j  < n_lateral_list[i].point_indices.size(); j += 2 )
		{
			mesh_endpoints.triangles.push_back(pair_index_to_mesh_index[n_lateral_list[i].point_indices[j]]);
			mesh_endpoints.triangles.push_back(pair_index_to_mesh_index[n_lateral_list[i].point_indices[0]]);
			mesh_endpoints.triangles.push_back(pair_index_to_mesh_index[n_lateral_list[i].point_indices[j + 1]]);
		}
	}

	//embed 
	for (size_t i = 0; i < mesh_endpoints.vertices.size(); i++)
	{
		mesh_endpoints.vertices[i].z = 0.0f;
		mesh_endpoints.vertices[i].x = embedding(i, 0) * 20;
		mesh_endpoints.vertices[i].y = embedding(i, 1) * 20;
	}
	mesh = mesh_endpoints;

	return embedding; 
}


// non auto generated.
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
	// for now !! 
	return A;
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