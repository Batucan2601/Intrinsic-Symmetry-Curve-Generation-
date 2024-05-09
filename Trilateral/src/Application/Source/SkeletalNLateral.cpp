#include "../Include/SkeletalNLateral.h"
#include "../Include/ShapeDiameter.h"
#include <eigen/Eigen/Dense>

#define SDF_PARAMETER 0.055f
static std::vector<float> pair_data;
static unsigned int vao;
static float maximum_sdf; // for sdf filter 

static std::vector<int> get_smallest_N(const std::vector<float>& distances, int N, int index , float constraint);
static bool compare_two_sdf_values(Mesh* m, Skeleton& skeleton, int index1, int index2, std::vector<float>& sdf_array);
static void get_maximum_sdf(std::vector<float>& sdf_array);

//only geodesic distances for skeletal end points 
SkeletalNLateral::SkeletalNLateral(Skeleton& skeleton, const std::vector<int>& point_indices, int N)
{
	// 1 - copy data 
	this->point_indices = point_indices;
	this->N = N;

	// 2- calculate geodesic distances
	std::vector<std::vector<float>> dijkstra_distances(N, std::vector<float>(N));
	std::vector<std::vector<int>> vertex_list;

	for (size_t i = 0; i < this->N; i++)
	{
		std::vector<int> vertex_list_i;
		std::vector<float> dijkstra_distances_i;
		skeleton_calculate_dijkstra(skeleton, skeleton.endPoints[this->point_indices[i]].index , vertex_list_i, dijkstra_distances_i);
		for (size_t j = i; j < this->N; j++)
		{
			if (i == j)
			{
				continue;
			}
			dijkstra_distances[i][j] = dijkstra_distances_i[skeleton.endPoints[this->point_indices[j]].index];
			dijkstra_distances[j][i] = dijkstra_distances_i[skeleton.endPoints[this->point_indices[j]].index];

			vertex_list.push_back(vertex_list_i);

		}
	}

	this->geodesic_distances = dijkstra_distances;
	this->predecessor_list = vertex_list;

	//get distance to the mid point
	std::vector<int> vertex_list_temp;
	std::vector<float> dijkstra_distances_temp;
	skeleton_calculate_dijkstra(skeleton, skeleton.mid_point_index, vertex_list_temp, dijkstra_distances_temp);

	//get minimum dist apart from itself
	int minimum_index = -1;
	float minimum_dist = INFINITY;
	for (size_t i = 0; i < dijkstra_distances_temp.size(); i++)
	{
		if (minimum_dist > dijkstra_distances_temp[i])
		{
			minimum_dist = dijkstra_distances_temp[i];
			minimum_index = i;
		}
	}

	this->geo_dist_to_skel_mid_point = minimum_dist;

}

float SkeletalNLateral_compareTwoSkeletalNLateral(SkeletalNLateral& nLateral1, SkeletalNLateral& nLateral2 , int N )
{
	//TODO YOU SHOULD ADD DISTANCE TO THE MIDPOINT
	//use geodesic distances
	// generate the first eigen vector
	Eigen::VectorXf nLateral1_vector(N + 1);
	for (size_t i = 0; i < N; i++)
	{
		nLateral1_vector(i) = nLateral1.geodesic_distances[0][i];
	}

	nLateral1_vector(N) = nLateral1.geo_dist_to_skel_mid_point;
	
	//generate permuation
	std::vector<unsigned int> permutation_vector;
	for (size_t i = 0; i < N; i++)
	{
		permutation_vector.push_back(i);
	}
	//all permutations for indices
	std::vector<std::vector<unsigned int>> all_permutations;
	std::vector<Eigen::VectorXf> nLateral2_vector_list;

	do {
		all_permutations.push_back(permutation_vector);
	} while (std::next_permutation(permutation_vector.begin(), permutation_vector.end()));

	for (size_t i = 0; i < all_permutations.size(); i++)
	{
		Eigen::VectorXf nLateral2_vector(N + 1 );
		for (size_t j = 0; j < all_permutations[i].size(); j++)
		{
			int index_1 = all_permutations[i][0];
			int index_2 = all_permutations[i][j];
			nLateral2_vector(j) = nLateral2.geodesic_distances[index_1][index_2];
		}
		nLateral2_vector(N) = nLateral2.geo_dist_to_skel_mid_point;
		nLateral2_vector_list.push_back(nLateral2_vector);
	}

	//now compare, find smallest and return the similarity
	float minimum_value = INFINITY;
	int minimum_index = -1;
	for (size_t i = 0; i < all_permutations.size(); i++)
	{
		float val = (nLateral1_vector - nLateral2_vector_list[i]).norm();
		if (minimum_value >  val )
		{
			minimum_value = val;
			minimum_index = i;
			
			//permutation vector
			/*for (size_t j = 0; j < all_permutations[i].s; j++)
			{

			}
			nLateral2.point_indices = */
		}
	}

	return minimum_value;
}

std::vector<std::pair<int, int>> SkeletalNLateral_compare_endpoints_with_SkeletalNlateral(Skeleton& skeleton, Mesh* m ,
int N, std::vector<float>& mesh_sdf_array)
{
	// declare the endpoint pairs 
	std::vector<std::pair<int, int>> end_point_pairs;
	//generate every Nlateral with closest points
	std::vector<SkeletalNLateral> skeletalNLateral_vec;
	
	get_maximum_sdf(mesh_sdf_array);
	for (size_t i = 0; i < skeleton.endPoints.size(); i++)
	{
		// get closest indices
		int current_index = skeleton.endPoints[i].index;
		std::vector<int> vertex_list;
		std::vector<float> dijkstra_distances;
		std::vector<int> smallest_endpoint_indices;
		skeleton_get_dijkstra_endpoints(skeleton,  current_index, vertex_list, dijkstra_distances);
		//get smallest N-1
		smallest_endpoint_indices = get_smallest_N(dijkstra_distances, N - 1, i , 1.0f);
		
		SkeletalNLateral skeletalNLateral(skeleton , smallest_endpoint_indices , N);
		
		skeletalNLateral_vec.push_back(skeletalNLateral);
	}

	for (size_t i = 0; i < skeletalNLateral_vec.size(); i++)
	{
		float best_val = INFINITY ;
		int best_index = -1;
		//compare each other
		for (size_t j = 0; j < skeletalNLateral_vec.size(); j++)
		{
			int end_point_indices_i = skeletalNLateral_vec[i].point_indices[0];
			int end_point_indices_j = skeletalNLateral_vec[j].point_indices[0];
			if (i == j)
			{
				continue;
			}
			//check if they do not have each other in point indices
			bool is_NLateralSame = false;
			for (size_t k = 0; k < N; k++)
			{
				if (skeletalNLateral_vec[i].point_indices[0] == skeletalNLateral_vec[j].point_indices[k])
				{
					is_NLateralSame = true;
					break;
				}
			}
			if(is_NLateralSame)
			{
				continue; 
			}
			// check the sdf
			int skeleton_index_i = skeleton.endPoints[end_point_indices_i].index;
			int skeleton_index_j = skeleton.endPoints[end_point_indices_j].index;
			int mesh_index1 = mesh_get_closest_index(m,skeleton.skeletonFormat[skeleton_index_i].point);
			int mesh_index2 = mesh_get_closest_index(m,skeleton.skeletonFormat[skeleton_index_j].point);
			bool is_sdf_similiar = compare_two_sdf_values(m, skeleton, end_point_indices_i, end_point_indices_j, mesh_sdf_array);
			
			if(!is_sdf_similiar)
			{
				continue; 
			}
			else
			{
				int debug; 
			}

			
			float norm =SkeletalNLateral_compareTwoSkeletalNLateral(skeletalNLateral_vec[i], skeletalNLateral_vec[j], N);
			if (best_val > norm /*&& skeletalNLateral_vec[i].point_indices[0] != skeletalNLateral_vec[j].point_indices[0]*/)
			{
				best_val = norm;
				best_index = j;
			}
		}
		end_point_pairs.push_back(std::pair<int, int>(skeletalNLateral_vec[i].point_indices[0] , skeletalNLateral_vec[best_index].point_indices[0]));
	}

	return end_point_pairs;
}

//constraint is used for selecting points with greater geodesic distances compared to constraint variable
static std::vector<int> get_smallest_N(const std::vector<float>& distances , int N , int index , float constraint)
{
	std::vector<int> smallest_indices;
	smallest_indices.push_back(index);
	for (size_t i = 0; i < N; i++)
	{
		int smallest_index = -1;
		float smallest_value = INFINITY;
		for (size_t j = 0; j < distances.size(); j++)
		{
			if (j == index) //itself 
			{
				continue; 
			}

			if (smallest_value > distances[j] && constraint > 0 &&  distances[j] > constraint )
			{
				//check if already included
				bool is_already_included = false;
				for (size_t k = 0; k < smallest_indices.size(); k++)
				{
					if (smallest_indices[k] == j)
					{
						is_already_included = true; 
						break;
					}
				}
				if (!is_already_included)
				{
					smallest_value = distances[j];
					smallest_index = j;
				}
			}

		}
		smallest_indices.push_back(smallest_index);
	}
	return smallest_indices; 
}



void SkeletalNLateral_generate_buffer(Skeleton& skeleton, std::vector<std::pair<int, int>>& pairs)
{

	unsigned int  vbo;
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);

	//forge it to a list

	for (size_t i = 0; i < pairs.size(); i++)
	{
		int endPoint1_index = pairs[i].first;
		int endPoint2_index = pairs[i].second;
		
		int point1_index = skeleton.endPoints[endPoint1_index].index;
		int point2_index = skeleton.endPoints[endPoint2_index].index;
		
		glm::vec3 p1 = skeleton.skeletonFormat[point1_index].point;
		glm::vec3 p2 = skeleton.skeletonFormat[point2_index].point;

		pair_data.push_back(p1.x);
		pair_data.push_back(p1.y);
		pair_data.push_back(p1.z);
		pair_data.push_back(0.0f);
		pair_data.push_back(255.0f);
		pair_data.push_back(0.0f);

		pair_data.push_back(p2.x);
		pair_data.push_back(p2.y);
		pair_data.push_back(p2.z);
		pair_data.push_back(0.0f);
		pair_data.push_back(255.0f);
		pair_data.push_back(0.0f);
	}

}
void SkeletalNLateral_buffer()
{
	glBindVertexArray(vao);
	glBufferData(GL_ARRAY_BUFFER, pair_data.size() * sizeof(float), &pair_data[0], GL_STATIC_DRAW);
	glBindVertexArray(0);
}
void SkeletalNLateral_draw(MeshFactory& mesh_fac , unsigned int shader_id)
{
	glBindVertexArray(vao);
	glm::mat4 model = mesh_fac.mesh_vec[0].model_mat;
	glm::mat4 MVP = mesh_fac.projection * mesh_fac.view * model;
	mesh_fac.mesh_vec[0].MVP = MVP;
	glUniformMatrix4fv(glGetUniformLocation(shader_id, "u_MVP"), 1, GL_FALSE, &MVP[0][0]);
	glDrawArrays(GL_LINES, 0, pair_data.size() / 6 );
	glBindVertexArray(0);

}

//index1 and index2 are skeleton indices, convert them to  mesh indices an compare their sdf
// the golden rule seems like if difference in SDF is more than 0.055MAX_SDF you return false 
static bool compare_two_sdf_values(Mesh* m, Skeleton& skeleton, int end_point_index1, int end_point_index2, std::vector<float>& sdf_array)
{

	if (fabs(sdf_array[end_point_index1] - sdf_array[end_point_index2])  > (SDF_PARAMETER * maximum_sdf) )
	{
		return false;
	}
	return true; 
}

static void get_maximum_sdf(std::vector<float>& sdf_array)
{
	int maximum_index = -1; 
	float maximum_dist = -INFINITY;
	for (size_t i = 0; i < sdf_array.size(); i++)
	{
		if (sdf_array[i] > maximum_dist)
		{
			maximum_dist = sdf_array[i];
			maximum_index = i; 
		}
	}
	maximum_sdf = maximum_dist;
}