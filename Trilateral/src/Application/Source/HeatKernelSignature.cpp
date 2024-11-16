#include "../Include/HeatKernelSignature.h"
#include "eigen/Eigen/EigenValues"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <numeric>
#include <fstream>
#include "../Include/Geodesic.h"

int no_of_samples = 10 ;
int no_of_eigen_val = 3;

void HKS_extract_kernel_signature(TrilateralMesh* m)
{
	Eigen::SparseMatrix<double>  cotangent_mat_sparse = cotangent_laplacian(*m);
	cotangent_mat_sparse = normalize_laplacian(cotangent_mat_sparse);

	check_if_matrix_symmetric(cotangent_mat_sparse);

	Spectra::SparseSymMatProd<double> op(cotangent_mat_sparse);

	// Create a SymEigsSolver object using the matrix operation object
	// Compute 3 eigenvalues, with a shift towards the smallest values
	//Spectra::SymEigsSolver< double, Spectra::SmallestAlge, Spectra::SparseSymMatProd<double> > eigs(&op, 3, 6);
	Eigen::Index num_of_eigenvalues = std::min(Eigen::Index(3), cotangent_mat_sparse.rows());
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> solver(op , num_of_eigenvalues,10);
	for (int k = 0; k < cotangent_mat_sparse.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(cotangent_mat_sparse, k); it; ++it) {
			if (std::isnan(it.value()) || std::isinf(it.value())) {
				std::cerr << "Invalid matrix entry at (" << it.row() << ", " << it.col() << "): " << it.value() << std::endl;
				return;
			}
		}
	}
	int max_iterations = 1000;  // Increase maximum number of iterations
	double tolerance = 1e-8;    // Set a smaller tolerance for convergence

	// Initialize and compute the eigenvalues
	solver.init();
	int nconv = solver.compute(Spectra::SortRule::SmallestAlge ,2000);
	Eigen::VectorXd eigenvalues;
	if (solver.info() == Spectra::CompInfo::Successful)
		eigenvalues = solver.eigenvalues();

	float minT = (abs(4 * log(10) / eigenvalues(eigenvalues.size() - 1)));
	float maxT = (abs(4 * log(10) / eigenvalues(1)));

	const double logTMin = log(minT);
	const double logTMax = log(maxT);

	const double tIncrement = (logTMax - logTMin) / (no_of_samples- 1);

	double sumRingAreas = std::accumulate(m->areas.begin(), m->areas.end(), 0.0);
	double sqrtSumRingAreas = sqrt(sumRingAreas);

	const size_t verSize = m->vertices.size();

	std::vector<std::vector<double >> point_descriptor(m->vertices.size());
	for (size_t v = 0; v < verSize; v++)
	{
		if (v < verSize - 1)
		{
			std::cout << "%" << (100 * v) / verSize << " completed for calculating hks descriptor" << "\r";
		}
		else
		{
			std::cout << "%" << 100 << " completed for calculating hks descriptor" << "\n";
		}

		//Initialize descriptor of the vertex
		std::vector<double> descriptor(no_of_samples);

		double currT = logTMin;
		for (size_t t = 0; t < no_of_samples; t++)
		{
			const double logScaleCurrT = exp(currT);
			double heatTrace = 0.0;
			double descVal = 0.0;
			for (size_t eig = 0; eig < no_of_eigen_val; eig++)
			{
				const double expInvEigVal = exp(sumRingAreas * abs(eigenvalues(eig)) * logScaleCurrT * (-1.0));
				const double EigFuncVal = sqrtSumRingAreas * eigenvalues(v, eig);
				const double sqEigFuncVal = EigFuncVal * EigFuncVal;
				descVal += (expInvEigVal * sqEigFuncVal);
				heatTrace += expInvEigVal;
			}
			descriptor[t] = descVal / heatTrace;
			currT += tIncrement;
		}
		point_descriptor[v] = descriptor;
	}


}

void HKS_read_kernel_signature(TrilateralMesh* m)
{
	//go read the hks file 
	std::string path = "../../Trilateral/Mesh/off/HKS/";
	std::string file_name = m->file_name;
	std::string to_insert = ".hks";
	file_name.insert(file_name.size() - 4, to_insert);
	path = path + file_name;
	std::ifstream file(path );  // Open the file for reading

	// Check if the file was successfully opened
	if (!file.is_open()) {
		std::cout << path << std::endl; 
		std::cerr << "Could not open the file!" << std::endl;
		return;
	}

	std::string line;
	// Extract integers from the stream
	int line_count = 0;
	// Read the file line by line
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		if (line_count > 1 &&  line_count < m->vertices.size() + 2)
		{
			float number;
			std::vector<float> numbers;
			// Extract integers from the stream
			while (iss >> number) {
				numbers.push_back(number);
			}
			m->normalized_heat_kernel_signature.push_back(numbers[3]); 
		}
		line_count += 1;
	}
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (m->normalized_heat_kernel_signature[i]< 1e-10)
		{
			m->normalized_heat_kernel_signature[i] += 1e-3; // dont make them 0 
		}
	}
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		m->raylib_mesh.colors[i * 4] = (unsigned char ) (m->normalized_heat_kernel_signature[i] * 255); 
		m->raylib_mesh.colors[i * 4 + 1] = 0; 
		m->raylib_mesh.colors[i * 4 + 2] = 0; 
		m->raylib_mesh.colors[i * 4 + 3] = 255; 
	}
	m->update_raylib_mesh();
	file.close();  // Close the file after reading
	return;
}

std::vector<std::vector<double>> HKS_compute_kernel(TrilateralMesh* m, std::pair<Eigen::VectorXd, Eigen::MatrixXd>& eigen_pairs , const std::vector<double>& timeSteps
, int time_step_no)
{
	int nPoints = m->vertices.size();
	int nTimeSteps = timeSteps.size();
	Eigen::MatrixXd hks = Eigen::MatrixXd::Zero(nPoints, nTimeSteps);
	double biggest = -1;

	eigen_pairs.first = eigen_pairs.first / eigen_pairs.first.maxCoeff();

	for (int t = 0; t < nTimeSteps; ++t) {
		double time = timeSteps[t];
		for (int i = 0; i < nPoints; ++i) {
			double hksValue = 0.0;
			for (int k = 0; k < eigen_pairs.first.size(); ++k) {
				double eigenValue = eigen_pairs.first(k);
				//std::cout << " 1 - " << eigenValue << std::endl;
				double expTerm = exp(-eigenValue * time);
				//std::cout << " 2 - " << expTerm << std::endl;
				hksValue += expTerm * pow(eigen_pairs.second(i, k), 2);
				//std::cout << " 3 - " << hksValue << std::endl;

			}
			hks(i,t) = hksValue;
			//std::cout << " 4 - " << hks(i, t) << std::endl;
			if (hks(i,t) > biggest)
			{
				biggest = hks(i,t);

			}
		}
	}
	if (hks.array().isNaN().any())
	{
		hks(0, 0) = 1;

		return std::vector<std::vector<double>>();
	}
	Eigen::VectorXd hks_total(nPoints);
	hks_total = hks_total.array() * 0;

	for (size_t i = 0; i < nPoints; i++)
	{
		hks_total(i) = hks_total(i) + hks.row(i).sum();
	}
	hks_total = (hks_total.array() - hks_total.minCoeff()) / (hks_total.maxCoeff() - hks_total.minCoeff());
	for (size_t i = 0; i < hks.cols(); i++)
	{
		hks.col(i) = (hks.col(i).array() - hks.col(i).minCoeff()) / (hks.col(i).maxCoeff()- hks.col(i).minCoeff());

	}
	//std::cout << hks.col(timeSteps.size() - 1);
	for (size_t i = 0; i < nPoints; i++)
	{
		m->raylib_mesh.colors[i * 4] = hks(i, time_step_no) * 255;
		//m->raylib_mesh.colors[i * 4] = hks_total(i) * 255;
		m->raylib_mesh.colors[i * 4 + 1] = 0;
		m->raylib_mesh.colors[i * 4 + 2] = 0;
		m->raylib_mesh.colors[i * 4 + 3] = 255;

	}
	m->update_raylib_mesh();
	return std::vector<std::vector<double>>();
}



void HKS_hks_on_descriptor(TrilateralMesh* m, TrilateralDescriptor& desc)
{
	//  1 - generate a sparse matrix with vertices inside
	TrilateralDescriptor_generate_mesh_inside(m, desc);
	Eigen::SparseMatrix<double> L = laplace_beltrami(&desc.m_inside);

	std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigen_pairs = laplace_beltrami_eigendecompose(L, 3);
	std::vector<double> time_steps = { 0.1, 0.2, 0.5, 1, 2 };
	
	int nPoints = desc.m_inside.vertices.size();
	int nTimeSteps = time_steps.size();
	std::vector<std::vector<double>> hks(nPoints, std::vector<double>(nTimeSteps, 0.0));
	double biggest = -1;

	for (int t = 0; t < nTimeSteps; ++t) {
		double time = time_steps[t];
		for (int i = 0; i < nPoints; ++i) {
			double hksValue = 0.0;
			for (int k = 0; k < eigen_pairs.first.size(); ++k) {
				double eigenValue = eigen_pairs.first(eigen_pairs.first.size() - k - 1);
				double expTerm = exp(-eigenValue * time);
				hksValue += expTerm * pow(eigen_pairs.second(i, eigen_pairs.first.size() - k - 1), 2);
			}
			hks[i][t] = hksValue;
			if (hks[i][t] > biggest)
			{
				biggest = hks[i][t];
			}
		}
	}

	// color
	for (size_t i = 0; i < nPoints; i++)
	{
		m->raylib_mesh.colors[desc.m_inside.mesh_to_desc_map[i] * 4] = (unsigned char)(hks[i][0] / biggest * 255);
		m->raylib_mesh.colors[desc.m_inside.mesh_to_desc_map[i] * 4 + 1] = 0;
		m->raylib_mesh.colors[desc.m_inside.mesh_to_desc_map[i] * 4 + 2] = 0;
		m->raylib_mesh.colors[desc.m_inside.mesh_to_desc_map[i] * 4 + 3] = 255;
	}
	m->update_raylib_mesh();
}

std::vector<std::pair<int, float>> HKS_sweep_distance(TrilateralMesh* m ,std::vector<std::pair<int,float>> pair  ,float distance)
{

	std::vector<std::pair<int, float>> extracted_pairs;
	while (pair.size() > 0)
	{
		std::pair<int, float> p = pair[0];
		extracted_pairs.push_back(p);
		pair.erase(pair.begin(), pair.begin() + 1);
		int index = p.first;
		std::vector<float> distances = Geodesic_dijkstra(*m, index);
		std::vector<unsigned int> close_vertices;
		for (size_t i = 0; i < pair.size(); i++)
		{
			if (distances[pair[i].first] < distance)
			{
				close_vertices.push_back(i);
			}
		}

		std::sort(close_vertices.begin(), close_vertices.end(), std::greater<unsigned int>());

		// Remove elements from the vector starting from the highest index
		for (int index : close_vertices) {
			if (index >= 0 && index < pair.size()) {
				pair.erase(pair.begin() + index);
			}
			else {
				std::cerr << "Invalid index: " << index << std::endl;
			}
		}

	}
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		m->raylib_mesh.colors[i * 4] = 0;
		m->raylib_mesh.colors[i * 4 + 1] = 0;
		m->raylib_mesh.colors[i * 4 + 2] = 0;
		m->raylib_mesh.colors[i * 4 + 3] = 255;
	}
	for (size_t i = 0; i < extracted_pairs.size(); i++)
	{
		int index = extracted_pairs[i].first;
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 0;
		m->raylib_mesh.colors[index * 4 + 2] = 255;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}
	m->update_raylib_mesh();

	return extracted_pairs;
}
std::vector<std::pair<int, float>> HKS_extraction_significant_points(TrilateralMesh* m, int P)
{
	int N = m->vertices.size();
	std::vector<std::pair<int, float>> indexedArr;

	for (int i = 0; i < N; ++i) {
		indexedArr.emplace_back(i, m->normalized_heat_kernel_signature[i]);
	}
	// Sort the vector by the second element (value)
	std::sort(indexedArr.begin(), indexedArr.end(), [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
		return a.second > b.second;
		});

	std::vector<std::pair<int, float>> arr;
	for (size_t i = 0; i < P; i++)
	{
		arr.push_back(indexedArr[i]);
	}

	for (size_t i = 0; i < P; i++)
	{
		int index = arr[i].first;
		m->raylib_mesh.colors[index * 4] = 0;
		m->raylib_mesh.colors[index * 4 + 1] = 0;
		m->raylib_mesh.colors[index * 4 + 2] = 255;
		m->raylib_mesh.colors[index * 4 + 3] = 255;
	}
	m->update_raylib_mesh();

	return arr;
}


std::vector<DvorakPairs> HKS_to_dvorak_pairs(TrilateralMesh* m,std::vector<std::pair<int, float>>& pairs)
{
	std::vector<DvorakPairs> dvorak_pairs; 
	for (size_t i = 0; i < pairs.size(); i++)
	{
		DvorakPairs dvorak;
		int index = pairs[i].first;
		dvorak.p_index = index ;
		float gaussian = gaussian_curvature(m, index);
		
		dvorak.gaussian_curv = gaussian;
		dvorak.p_index = index;
		dvorak_pairs.push_back(dvorak);
	}
	return dvorak_pairs;
}