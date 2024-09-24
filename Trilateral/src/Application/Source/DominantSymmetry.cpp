#include "../Include/DominantSymmetry.h"
#include "../Include/glm/glm.hpp"
#include "../../External/Include/eigen/Eigen/Core"
#include "../../External/Include/eigen/Eigen/Eigenvalues" 
#include "../Include/Sampling.h"
#include "../Include/TrilateralMap.h"
#include "../Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
using Eigen::MatrixXd;

static Plane rotate_plane(Plane plane, float rotation_degree);

Plane generate_dominant_symmetry_plane(int seletected_mesh, MeshFactory& mesh_fac , int sym_iter_no) 
{
	Mesh mesh = mesh_fac.mesh_vec[seletected_mesh];
	
	Plane plane = generate_dominant_symmetry_plane(mesh , sym_iter_no);

	float s = mesh.vertices.size();
	glm::vec3 m(0.0f, 0.0f, 0.0f);
	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		m += mesh.vertices[i];
	}

	m = m / s;
	Mesh plane_mesh = generate_mesh_from_plane(&plane, &m);
	mesh_fac.add_mesh(plane_mesh);

	return plane;
	
}
Plane generate_dominant_symmetry_plane(Mesh mesh , int sym_iter_no)
{
	// generate PCA weights are same and 1 for now 
	float s = mesh.vertices.size();
	glm::vec3 m(0.0f, 0.0f, 0.0f);
	for (size_t i = 0; i < mesh.vertices.size(); i++)
	{
		m += mesh.vertices[i];
	}
	m = m / s;

	// surfel weights 
	std::vector<float> surfel_weights(mesh.vertices.size() , 0);
	surfel_weights = mesh_point_surfel_normalized(&mesh);


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

		Eigen::MatrixXd  Co_i = surfel_weights[i] * pi * pi.transpose();
		Co = Co + Co_i;
	}
	Co = Co / s;

	//// get the eigenvectors 
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Co);
	Eigen::MatrixXd eigen_vecs = es.eigenvectors().real();
	Eigen::VectorXd eigen_values = es.eigenvalues().real();

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


	float eigen_val1 = es.eigenvalues().col(0).real()(0);
	float eigen_val2 = es.eigenvalues().col(0).real()(1);
	float eigen_val3 = es.eigenvalues().col(0).real()(2);
	int eigen_second_best_index = -1;

	//get second best 
	if ( (eigen_val1 > eigen_val2 && eigen_val1 < eigen_val3) || 
		(eigen_val1 > eigen_val3 && eigen_val1 < eigen_val2))
	{
		eigen_second_best_index  = 0;
	}
	if ( (eigen_val2 > eigen_val1 && eigen_val2 < eigen_val3) || 
		(eigen_val2 > eigen_val3 && eigen_val2 < eigen_val1))
	{
		eigen_second_best_index = 1;

	}
	if ( (eigen_val3 > eigen_val1 && eigen_val3 < eigen_val2) ||
	(eigen_val3 > eigen_val3 && eigen_val3 < eigen_val1) )
	{
		eigen_second_best_index = 2;

	}

	Plane best_plane = planes[eigen_second_best_index];
	best_plane.normal = glm::normalize(best_plane.normal);
	int sym_no = 0; 
	//for( int p = 0; p < sym_iter_no; p++ )
	while( 1 )
	{
		std::vector<float> di_vec;
		for (size_t i = 0; i < mesh.vertices.size(); i++)
		{
			glm::vec3 s_ir = symmetry_point_from_plane(&best_plane, &mesh.vertices[i]);
			float s_ir_status = get_point_status_from_plane(&best_plane, &s_ir);
			float d_i= INFINITY;
			int minimum_index = -1;
			for (size_t j = 0; j < mesh.vertices.size(); j++)
			{
				if (i == j)
				{
					continue;
				}
				float s_j = get_point_status_from_plane(&best_plane, &mesh.vertices[j]);
				if (s_j * s_ir_status < 0) //different side 
				{
					continue;
				}
				float distance = glm::distance(s_ir, mesh.vertices[j]);
				if (distance < d_i)
				{
					d_i  = distance;
				}
			}
			di_vec.push_back(d_i);
		}
		std::vector<float> di_vec_sorted = di_vec;

		std::sort(di_vec_sorted.begin(), di_vec_sorted.end());
		float di_median = di_vec_sorted[di_vec_sorted.size() / 2];
		float c = 1.5;
		float sigma = c * di_median;

		//redo surfels
		for (size_t i = 0; i < mesh.vertices.size(); i++)
		{
			float d_i = di_vec[i];
			if (d_i > sigma)
			{
				surfel_weights[i] = (2 * sigma * sigma) /  pow((sigma*sigma + d_i * d_i ),2);
			}
			else
			{
				surfel_weights[i] = 0;
			}
		}

		//redo PCA
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

			Eigen::MatrixXd  Co_i = surfel_weights[i] * pi * pi.transpose();
			Co = Co + Co_i;
		}
		Co = Co / s;

		//// get the eigenvectors 
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Co);
		Eigen::MatrixXd eigen_vecs = es.eigenvectors().real();
		Eigen::VectorXd eigen_values = es.eigenvalues().real();

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


		float eigen_val1 = es.eigenvalues().col(0).real()(0);
		float eigen_val2 = es.eigenvalues().col(0).real()(1);
		float eigen_val3 = es.eigenvalues().col(0).real()(2);
		int eigen_second_best_index = -1;

		//get second best 
		if ((eigen_val1 > eigen_val2 && eigen_val1 < eigen_val3) ||
			(eigen_val1 > eigen_val3 && eigen_val1 < eigen_val2))
		{
			eigen_second_best_index = 0;
		}
		if ((eigen_val2 > eigen_val1 && eigen_val2 < eigen_val3) ||
			(eigen_val2 > eigen_val3 && eigen_val2 < eigen_val1))
		{
			eigen_second_best_index = 1;

		}
		if ((eigen_val3 > eigen_val1 && eigen_val3 < eigen_val2) ||
			(eigen_val3 > eigen_val3 && eigen_val3 < eigen_val1))
		{
			eigen_second_best_index = 2;

		}
		//check convergence
		//get their formula
		float A_A = best_plane.normal.x;
		float A_B = best_plane.normal.y;
		float A_C = best_plane.normal.z;
		float A_D = - ( best_plane.normal.x * best_plane.point.x + best_plane.normal.y * best_plane.point.y + best_plane.normal.z * best_plane.point.z);

		float B_A = planes[eigen_second_best_index].normal.x;
		float B_B = planes[eigen_second_best_index].normal.y;
		float B_C = planes[eigen_second_best_index].normal.z;
		float B_D = -(planes[eigen_second_best_index].normal.x * planes[eigen_second_best_index].point.x + 
		planes[eigen_second_best_index].normal.y * planes[eigen_second_best_index].point.y + 
		planes[eigen_second_best_index].normal.z * planes[eigen_second_best_index].point.z);
		
		// convergence formula
		float convergence = sqrt(pow((A_A - B_A), 2) + pow((A_B - B_B), 2) + pow((A_C - B_C), 2) + pow((A_D - B_D), 2) );
		sym_no++;
		if (convergence < 0.1)
		{
			break; 
		}

		best_plane = planes[eigen_second_best_index];
		best_plane.normal = glm::normalize(best_plane.normal);
}
	


	return best_plane;
}
Plane generate_dominant_symmetry_plane(const glm::vec3& plane_point, Mesh mesh)
{

	// generate PCA weights are same and 1 for now 
	float s = mesh.vertices.size();
	glm::vec3 m(0.0f, 0.0f, 0.0f);
	int N = mesh.vertices.size();
	m = plane_point;

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
	for (size_t i = 0; i < N; i++)
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
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Co);
	Eigen::MatrixXd eigen_vecs = es.eigenvectors().real();
	Eigen::VectorXd eigen_values = es.eigenvalues().real();

	double biggest_value = -INFINITY;
	int biggest_index = -1;
	//get the best eigen value
	for (size_t i = 0; i < eigen_values.rows(); i++)
	{
		if (biggest_value < (float)eigen_values(i))
		{
			biggest_value = (float)eigen_values(i);
			biggest_index = i;
		}
	}

	
	// generate the 3 planes
	Plane planes[3];
	planes[0].normal = glm::vec3(eigen_vecs.col(0).real()(0), eigen_vecs.col(0).real()(1), eigen_vecs.col(0).real()(2));
	planes[0].point = m;
	planes[1].normal = glm::vec3(eigen_vecs.col(1).real()(0), eigen_vecs.col(1).real()(1), eigen_vecs.col(1).real()(2));
	planes[1].point = m;
	planes[2].normal = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	planes[2].point = m;

	std::vector<float> distances[3];
	float minimum_distances[3] = { 0 , 0 ,0 }; //for each plane 
	for (size_t p = 0; p < 3; p++) //for all planes 
	{
		for (size_t i = 0; i < mesh.vertices.size(); i++)
		{
			float minimum_of_di = INFINITY;
			// generate sir
			glm::vec3 s_ir = symmetry_point_from_plane(&planes[p], &mesh.vertices[i]);
			float s_ir_status = get_point_status_from_plane(&planes[p], &s_ir);
			for (size_t j = 0; j < mesh.vertices.size(); j++)
			{
				if (s_ir_status * get_point_status_from_plane(&planes[p], &mesh.vertices[j]) >= 0)
				{
					float dist = glm::distance(s_ir, mesh.vertices[j]);
					if (dist < minimum_of_di)
					{
						minimum_of_di = dist;
					}
				}
			}
			distances[p].push_back(minimum_of_di); //for median calculation
			minimum_distances[p] += minimum_of_di;
		}
		//minimum_distances[p] = minimum_of_di;
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

//#pragma region rotation around 3 degree
//	float symmetry_correspondence_score = -INFINITY;
//	int symmetry_correspondence_index = 0;
//	Plane rot_planes[9];
//	for (int p = 0; p < 3; p++)
//	{
//		// itself
//		// 1 - rotate +3 degree
//		// 2 - rotate -3 degree 
//		float degree = -3;
//		for (int i = 0; i < 3; i++)
//		{
//			Plane rotated_plane = rotate_plane(planes[p], degree);
//			rot_planes[p * 3 + i] = rotated_plane;
//			float symm_score = generate_symmetry_score( mesh ,&rotated_plane);
//			if (symm_score > symmetry_correspondence_score)
//			{
//				symmetry_correspondence_score = symm_score;
//				symmetry_correspondence_index = p * 3 + i;
//			}
//			degree += 3; 
//		}
//	}
//	// translate in both normal and normal's opposite direction 
//	Plane best_plane = rot_planes[symmetry_correspondence_index];
//	float step_size = 0.01f;
//	float prev_sym_score = symmetry_correspondence_score;
//	Plane temp_best_plane_normal_dir = best_plane;
//	float sym_score_normal_dir = 0;
//	while (1)
//	{
//		// move plane by normal with step size
//		temp_best_plane_normal_dir.point += glm::normalize(temp_best_plane_normal_dir.normal) * step_size;
//		float symm_score = generate_symmetry_score(mesh, &temp_best_plane_normal_dir);
//		if (symm_score < prev_sym_score)
//		{
//			sym_score_normal_dir = symm_score;
//			temp_best_plane_normal_dir.point -= glm::normalize(temp_best_plane_normal_dir.normal) * step_size;
//			break;
//		}
//		prev_sym_score = symm_score;
//	}
//	prev_sym_score = symmetry_correspondence_score;
//	Plane temp_best_plane_opposite_normal_dir = best_plane;
//	float sym_score_opposite_normal_dir = 0;
//	while (1)
//	{
//		// move plane by normal with step size
//		temp_best_plane_opposite_normal_dir.point += glm::normalize(temp_best_plane_opposite_normal_dir.normal) * step_size;
//		float symm_score = generate_symmetry_score(mesh, &temp_best_plane_opposite_normal_dir);
//		if (symm_score < prev_sym_score)
//		{
//			sym_score_opposite_normal_dir = symm_score;
//			temp_best_plane_opposite_normal_dir.point -= glm::normalize(temp_best_plane_opposite_normal_dir.normal) * step_size;
//			break;
//		}
//		prev_sym_score = symm_score;
//	}
//	if (sym_score_normal_dir > sym_score_opposite_normal_dir)
//	{
//		return temp_best_plane_normal_dir;
//	}
//	else
//	{
//		return temp_best_plane_opposite_normal_dir;
//	}
//	return rot_planes[symmetry_correspondence_index];
//#pragma endregion 
	//just return the one, easy way 
	 return planes[smallest_dist_index];

#pragma region iteration
	//Plane best_plane_A = planes[smallest_dist_index];
	//
	//
	//// last part
	//// do this one more time for second plane

	//// median
	//std::sort(distances[smallest_dist_index].begin(), distances[smallest_dist_index].end());
	//float median_di = distances[smallest_dist_index][distances[smallest_dist_index].size() / 2];
	//float c = 1.5;
	//float sigma = c * median_di;

	//// gave every point a weight
	//std::vector<float> weights(N);
	//for (size_t i = 0; i < N; i++)
	//{
	//	if (distances[smallest_dist_index][i] >= sigma)
	//	{
	//		weights[i] = 0;
	//	}
	//	else
	//	{
	//		weights[i] = (2 * pow(sigma, 2)) / (pow(pow(distances[smallest_dist_index][i], 2) + pow(sigma, 2), 2));
	//	}
	//}

	//// redo symmetry plane with new weights
	//Co(0, 0) = 0;
	//Co(0, 1) = 0;
	//Co(0, 2) = 0;
	//Co(1, 0) = 0;
	//Co(1, 1) = 0;
	//Co(1, 2) = 0;
	//Co(2, 0) = 0;
	//Co(2, 1) = 0;
	//Co(2, 2) = 0;
	//for (size_t i = 0; i < N; i++)
	//{
	//	glm::vec3 pi_m;
	//	pi_m = mesh.vertices[i] - m;
	//	Eigen::VectorXd pi(3);
	//	pi(0) = pi_m.x;
	//	pi(1) = pi_m.y;
	//	pi(2) = pi_m.z;

	//	pi = pi * weights[i];
	//	Eigen::MatrixXd  Co_i = pi * pi.transpose();
	//	Co = Co + Co_i;
	//}
	//Co = Co / s;

	//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esb(Co);
	//eigen_vecs = esb.eigenvectors().real();
	//eigen_values = esb.eigenvalues().real();
	//
	//Plane planesB[3];
	//planesB[0].normal = glm::vec3(eigen_vecs.col(0).real()(0), eigen_vecs.col(0).real()(1), eigen_vecs.col(0).real()(2));
	//planesB[0].point = m;
	//planesB[1].normal = glm::vec3(eigen_vecs.col(1).real()(0), eigen_vecs.col(1).real()(1), eigen_vecs.col(1).real()(2));
	//planesB[1].point = m;
	//planesB[2].normal = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	//planesB[2].point = m;

	//distances[0].clear();
	//distances[1].clear();
	//distances[2].clear();

	//weights.clear();
	//minimum_distances[0] = 0; //for each plane 
	//minimum_distances[1] = 0; 
	//minimum_distances[2] = 0; 
	//for (size_t p = 0; p < 3; p++) //for all planes 
	//{
	//	float minimum_of_di = INFINITY;
	//	for (size_t i = 0; i < mesh.vertices.size(); i++)
	//	{
	//		// generate sir
	//		glm::vec3 s_ir = symmetry_point_from_plane(&planesB[p], &mesh.vertices[i]);
	//		float s_ir_status = get_point_status_from_plane(&planesB[p], &s_ir);
	//		for (size_t j = 0; j < mesh.vertices.size(); j++)
	//		{
	//			if (s_ir_status * get_point_status_from_plane(&planesB[p], &mesh.vertices[j]) >= 0)
	//			{
	//				float dist = glm::distance(s_ir, mesh.vertices[j]);
	//				if (dist < minimum_of_di)
	//				{
	//					minimum_of_di = dist;
	//				}
	//			}
	//		}
	//		distances[p].push_back(minimum_of_di); //for median calculation
	//		minimum_distances[p] += minimum_of_di;
	//	}
	//	//minimum_distances[p] = minimum_of_di;
	//}
	//smallest_dist_index = 0;
	//smallest_dist = INFINITY;
	//for (size_t i = 0; i < 3; i++)
	//{
	//	if (smallest_dist > minimum_distances[i])
	//	{
	//		smallest_dist = minimum_distances[i];
	//		smallest_dist_index = i;
	//	}
	//}
	//Plane best_plane_B =  planesB[smallest_dist_index];

	//glm::vec4 plane_coeff_A;
	//glm::vec4 plane_coeff_B;
	//get_coefficients_from_plane(best_plane_A , plane_coeff_A.x , plane_coeff_A.y , plane_coeff_A.z , plane_coeff_A.w);
	//get_coefficients_from_plane(best_plane_B , plane_coeff_B.x , plane_coeff_B.y , plane_coeff_B.z , plane_coeff_B.w);
	//
	//glm::vec4 dif_plane = plane_coeff_A - plane_coeff_B;
	//float plane_dif_epsilon = glm::length(dif_plane);
	//const float epsilon_const = 0.01;
	//
	//Plane plane_current = best_plane_B; // useless, just for initializations sake 
	//Plane plane_prev = best_plane_B;
	////base case 
	//if (plane_dif_epsilon < epsilon_const)
	//{
	//	return plane_prev;
	//}
	//while (plane_dif_epsilon < epsilon_const)
	//{
	//	std::sort(distances[smallest_dist_index].begin(), distances[smallest_dist_index].end());
	//	float median_di = distances[smallest_dist_index][distances[smallest_dist_index].size() / 2];
	//	float c = 1.5;
	//	float sigma = c * median_di;

	//	// gave every point a weight
	//	std::vector<float> weights(N);
	//	for (size_t i = 0; i < N; i++)
	//	{
	//		if (distances[smallest_dist_index][i] >= sigma)
	//		{
	//			weights[i] = 0;
	//		}
	//		else
	//		{
	//			weights[i] = (2 * pow(sigma, 2)) / (pow(pow(distances[smallest_dist_index][i], 2) + pow(sigma, 2), 2));
	//		}
	//	}

	//	// redo symmetry plane with new weights
	//	Co(0, 0) = 0;
	//	Co(0, 1) = 0;
	//	Co(0, 2) = 0;
	//	Co(1, 0) = 0;
	//	Co(1, 1) = 0;
	//	Co(1, 2) = 0;
	//	Co(2, 0) = 0;
	//	Co(2, 1) = 0;
	//	Co(2, 2) = 0;
	//	for (size_t i = 0; i < N; i++)
	//	{
	//		glm::vec3 pi_m;
	//		pi_m = mesh.vertices[i] - m;
	//		Eigen::VectorXd pi(3);
	//		pi(0) = pi_m.x;
	//		pi(1) = pi_m.y;
	//		pi(2) = pi_m.z;

	//		pi = pi * weights[i];
	//		Eigen::MatrixXd  Co_i = pi * pi.transpose();
	//		Co = Co + Co_i;
	//	}
	//	Co = Co / s;

	//	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esb(Co);
	//	eigen_vecs = esb.eigenvectors().real();
	//	eigen_values = esb.eigenvalues().real();

	//	Plane planesB[3];
	//	planesB[0].normal = glm::vec3(eigen_vecs.col(0).real()(0), eigen_vecs.col(0).real()(1), eigen_vecs.col(0).real()(2));
	//	planesB[0].point = m;
	//	planesB[1].normal = glm::vec3(eigen_vecs.col(1).real()(0), eigen_vecs.col(1).real()(1), eigen_vecs.col(1).real()(2));
	//	planesB[1].point = m;
	//	planesB[2].normal = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	//	planesB[2].point = m;

	//	distances[0].clear();
	//	distances[1].clear();
	//	distances[2].clear();

	//	minimum_distances[0] = 0; //for each plane 
	//	minimum_distances[1] = 0;
	//	minimum_distances[2] = 0;
	//	for (size_t p = 0; p < 3; p++) //for all planes 
	//	{
	//		float minimum_of_di = INFINITY;
	//		for (size_t i = 0; i < mesh.vertices.size(); i++)
	//		{
	//			// generate sir
	//			glm::vec3 s_ir = symmetry_point_from_plane(&planesB[p], &mesh.vertices[i]);
	//			float s_ir_status = get_point_status_from_plane(&planesB[p], &s_ir);
	//			for (size_t j = 0; j < mesh.vertices.size(); j++)
	//			{
	//				if (s_ir_status * get_point_status_from_plane(&planesB[p], &mesh.vertices[j]) >= 0)
	//				{
	//					float dist = glm::distance(s_ir, mesh.vertices[j]);
	//					if (dist < minimum_of_di)
	//					{
	//						minimum_of_di = dist;
	//					}
	//				}
	//			}
	//			distances[p].push_back(minimum_of_di); //for median calculation
	//			minimum_distances[p] += minimum_of_di;
	//		}
	//		//minimum_distances[p] = minimum_of_di;
	//	}
	//	smallest_dist_index = 0;
	//	smallest_dist = INFINITY;
	//	for (size_t i = 0; i < 3; i++)
	//	{
	//		if (smallest_dist > minimum_distances[i])
	//		{
	//			smallest_dist = minimum_distances[i];
	//			smallest_dist_index = i;
	//		}
	//	}
	//	Plane plane_current = planesB[smallest_dist_index];

	//	glm::vec4 plane_coeff_cur;
	//	glm::vec4 plane_coeff_prev;
	//	get_coefficients_from_plane(plane_current, plane_coeff_cur.x, plane_coeff_cur.y, plane_coeff_cur.z, plane_coeff_cur.w);
	//	get_coefficients_from_plane(plane_prev, plane_coeff_prev.x, plane_coeff_prev.y, plane_coeff_prev.z, plane_coeff_prev.w);

	//	glm::vec4 dif_plane = plane_coeff_cur - plane_coeff_prev;
	//	plane_dif_epsilon = glm::length(dif_plane);

	//	plane_prev = plane_current;
	//}

	//return plane_current;
#pragma endregion
}
/* We separate the mesh into two in order to get symmetry sets
* m1 represents the points where you get the + sign when you plug the vertices in the plane
* m2 represents the minus sign 
*/
void dom_sym_generate_two_separate_mesh_using_dominant_symmetry_plane(Plane plane, Mesh* mesh_to_be_separated, Mesh* m1, Mesh* m2 , std::vector<int>* indices_for_m1 , std::vector<int>* indices_for_m2)
{
	indices_for_m1->resize(mesh_to_be_separated->vertices.size());
	indices_for_m2->resize(mesh_to_be_separated->vertices.size());
	std::fill(indices_for_m1->begin(), indices_for_m1->end(), -1);
	std::fill(indices_for_m2->begin(), indices_for_m2->end(), -1);
	// separate vertices
	for (size_t i = 0; i < mesh_to_be_separated->vertices.size(); i++)
	{
		if (get_point_status_from_plane( &plane , &mesh_to_be_separated->vertices[i]) > 0 )
		{
			(*indices_for_m1)[i] = m1->vertices.size();

			m1->vertices.push_back(mesh_to_be_separated->vertices[i]);
		}
		else
		{
			(*indices_for_m2)[i] = m2->vertices.size();

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

			if (sign_i  * sign_j >= 0 )
			{
				if (sign_i >= 0)
				{
					adjacencies_m1.push_back(std::pair <int,float>((*indices_for_m1)[mesh_to_be_separated->adjacenies[i][j].first], glm::distance(mesh_to_be_separated->vertices[i], mesh_to_be_separated->vertices[j])));
				}
				else if (sign_i < 0)
				{
					adjacencies_m2.push_back(std::pair <int, float>((*indices_for_m2)[mesh_to_be_separated->adjacenies[i][j].first], glm::distance(mesh_to_be_separated->vertices[i], mesh_to_be_separated->vertices[j])));
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

		if (sign_i >= 0 && sign_i_1 >= 0 && sign_i_2 >= 0)
		{
			m1->triangles.push_back((*indices_for_m1)[mesh_to_be_separated->triangles[i]]);
			m1->triangles.push_back((*indices_for_m1)[mesh_to_be_separated->triangles[i+1]]);
			m1->triangles.push_back((*indices_for_m1)[mesh_to_be_separated->triangles[i+2]]);
		}
		else if (sign_i < 0 && sign_i_1 < 0 && sign_i_2 < 0)
		{
			m2->triangles.push_back((*indices_for_m2)[mesh_to_be_separated->triangles[i]]);
			m2->triangles.push_back((*indices_for_m2)[mesh_to_be_separated->triangles[i+1]]);
			m2->triangles.push_back((*indices_for_m2)[mesh_to_be_separated->triangles[i+2]]);
		}
	}

}

void match_two_meshes_with_fps(Mesh* selected_mesh ,  Mesh* m1, Mesh* m2, std::vector<int>* indices_for_m1, std::vector<int>* indices_for_m2 , int no_of_samples)
{
	std::vector<unsigned int> sampled_points_m1 =  furthest_point_sampling(m1, no_of_samples);
	std::vector<unsigned int> sampled_points_m2 =  furthest_point_sampling(m2, no_of_samples);
	std::vector<TrilateralDescriptor> trilateral_desc_m1 = get_trilateral_points_using_closest_pairs(m1, sampled_points_m1);
	std::vector<TrilateralDescriptor> trilateral_desc_m2 = get_trilateral_points_using_closest_pairs(m2, sampled_points_m2);
	std::vector<TrilateralDescriptor> trilateral_desc_original_mesh;

	float trilteralW1 = 1;
	float trilteralW2 = 1;
	float trilteralW3 = 1;
	float trilteralW4 = 1;
	//use mapped functions to return to original
	for (size_t i = 0; i < trilateral_desc_m1.size(); i++)
	{
		trilateral_desc_m1[i].p1 = (*indices_for_m1)[trilateral_desc_m1[i].p1];
		trilateral_desc_m1[i].p2 = (*indices_for_m1)[trilateral_desc_m1[i].p2];
		trilateral_desc_m1[i].p3 = (*indices_for_m1)[trilateral_desc_m1[i].p3];
		
		trilateral_desc_original_mesh.push_back(trilateral_desc_m1[i]);
	}

	for (size_t i = 0; i < trilateral_desc_m2.size(); i++)
	{
		trilateral_desc_m2[i].p1 = (*indices_for_m1)[trilateral_desc_m2[i].p1];
		trilateral_desc_m2[i].p2 = (*indices_for_m1)[trilateral_desc_m2[i].p2];
		trilateral_desc_m2[i].p3 = (*indices_for_m1)[trilateral_desc_m2[i].p3];
		trilateral_desc_original_mesh.push_back(trilateral_desc_m2[i]);
	}

	std::vector<std::pair<unsigned int,unsigned int>> pairs =  point_match_trilateral_weights(selected_mesh, trilateral_desc_original_mesh, trilteralW1, trilteralW2, trilteralW3);
	display_accuracy(selected_mesh ,  pairs);
}

std::vector<TrilateralDescriptor> match_two_meshes_with_fps(Mesh* selected_mesh, Plane* plane, int no_of_samples)
{
	std::vector<unsigned int> sampled_points = furthest_point_sampling(selected_mesh, no_of_samples);
	
	std::vector<unsigned int> sampled_points_plane_right;
	std::vector<unsigned int> sampled_points_plane_left;
	

	for (size_t i = 0; i < no_of_samples; i++)
	{
		float sign = get_point_status_from_plane(plane ,&selected_mesh->vertices[sampled_points[i]]);
		if (sign >=  0)
		{
			sampled_points_plane_right.push_back(sampled_points[i]);
		}
		else
		{
			sampled_points_plane_left.push_back(sampled_points[i]);
		}
	}
	
	std::vector<TrilateralDescriptor> trilateral_desc_vec_right = get_trilateral_points_using_closest_pairs(selected_mesh, sampled_points_plane_right);
	std::vector<TrilateralDescriptor> trilateral_desc_vec_left = get_trilateral_points_using_closest_pairs(selected_mesh, sampled_points_plane_left);
	
	std::vector<TrilateralDescriptor> trilateral_desc_vec_concat = trilateral_desc_vec_right;
	trilateral_desc_vec_concat.insert(trilateral_desc_vec_concat.end(), trilateral_desc_vec_left.begin(), trilateral_desc_vec_left.end());

	return trilateral_desc_vec_concat;
}

static Plane rotate_plane(Plane plane , float rotation_degree)
{
	Plane rotated_plane = plane;

	//rotate normal around point m ? 

	glm::mat4 rot_mat = glm::rotate(glm::mat4(1.0f), glm::radians(rotation_degree), plane.point);
	
	rotated_plane.normal = glm::vec3(glm::vec4(rotated_plane.normal,1) * rot_mat);

	return rotated_plane;
	
}