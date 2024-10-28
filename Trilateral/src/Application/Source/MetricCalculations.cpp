#include "../Include/MetricCalculations.h"
#include "../Include/MeshFactory.h"
#include "../Include/Geodesic.h"

#define _USE_MATH_DEFINES
#include <math.h>
//v6'nin ground-truth simetrik noktasi v66 ise ve senin methodun v6->v77'ye gonderdiyse o zaman geodesic(v66, v77) costun olacak;
float Metric_get_geodesic_cost(TrilateralMesh* m, unsigned int point_index1, unsigned int calculated_index1_correspondence, bool isNormalized)
{
	unsigned int ground_truth = 0;

	ground_truth = m->ground_truth_symmetry_pairs[point_index1];

	std::vector<float> distance_matrix_for_ground_truth = Geodesic_dijkstra(*m, ground_truth);


	float dif_ground_vs_corresp = distance_matrix_for_ground_truth[calculated_index1_correspondence];
	
	//get the maximum of distance_matrix_for_ground_truth for normalization
	int max_index = -1;
	float max_value = -INFINITY;
	for (size_t i = 0; i < distance_matrix_for_ground_truth.size(); i++)
	{
		if (distance_matrix_for_ground_truth[i] > max_value)
		{
			max_value = distance_matrix_for_ground_truth[i];
			max_index = i; 
		}
	}
	if (isNormalized)
	{
		dif_ground_vs_corresp /= max_value; //return the value normalized
	}
	return dif_ground_vs_corresp; 
}
float Metric_get_geodesic_cost_with_list(TrilateralMesh* m, std::vector<unsigned int> point_indices, std::vector<unsigned int> calculated_index_correspondence_list)
{
	unsigned int N = point_indices.size();
	float error = 0;
	std::vector<unsigned int> ground_truth_indices; 
	for (size_t i = 0; i < N; i++)
	{
		float single_error = Metric_get_geodesic_cost( m , point_indices[i], calculated_index_correspondence_list[i] , true);
		error += single_error;
	}
	return error/N ; // return average error
}

float Metric_get_geodesic_cost(TrilateralMesh* m )
{
	unsigned int N = m->calculated_symmetry_pairs.size();
	float error = 0;
	
	// just in case 
	read_symmetry_format((char*)"../../Trilateral/TrilateralMesh/off/sym.txt", m);


	for (size_t i = 0; i < N; i++)
	{
		int index = m->calculated_symmetry_pairs[i].first;
		int calculated_pair = m->calculated_symmetry_pairs[i].second;

		float single_error = Metric_get_geodesic_cost(m, index, calculated_pair, true );
		error += single_error;
	}
	return error / N; // return average error
}

// to be fair it is really hard to sample correct indices, therefore these two metrics
// probably wont use in the end 
std::vector<std::pair<unsigned int , unsigned int> > Metric_get_correct_pairs(TrilateralMesh* m)
{
	std::vector<std::pair<unsigned int, unsigned int>> correct_pairs; 
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index = m->calculated_symmetry_pairs[i].first;
		int calculated_symmetry_index = m->calculated_symmetry_pairs[i].second;
		int ground_truth = m->ground_truth_symmetry_pairs[index];
		if (ground_truth == calculated_symmetry_index)
		{
			correct_pairs.push_back(std::pair<unsigned int, unsigned int >(index, calculated_symmetry_index));
		}
	}
	return correct_pairs;
}

std::vector<std::pair<unsigned int, unsigned int> > Metric_get_incorrect_pairs(TrilateralMesh* m)
{
	std::vector<std::pair<unsigned int, unsigned int>> incorrect_pairs;
	for (size_t i = 0; i < m->calculated_symmetry_pairs.size(); i++)
	{
		int index = m->calculated_symmetry_pairs[i].first;
		int calculated_symmetry_index = m->calculated_symmetry_pairs[i].second;
		int ground_truth = m->ground_truth_symmetry_pairs[index];
		if (ground_truth != calculated_symmetry_index)
		{
			incorrect_pairs.push_back(std::pair<unsigned int, unsigned int >(index, calculated_symmetry_index));
		}
	}
	return incorrect_pairs;
}


/*Correspondence Rate : The percentage of correspondences with geodesic dist
qance between the computed and ground - truth correspondence of the vertex less than
area(M)
PI*N.We take N = 20 as used in MT.
*/

float Metric_get_correspondance_rate(TrilateralMesh* m)
{
	float threshold = sqrtf(m->mesh_area / (M_PI * m->vertices.size()));
	float percentage = 0;
	int N = m->calculated_symmetry_pairs.size();
	int passed_point_count = 0;
	for (size_t i = 0; i < N; i++)
	{
		int index = m->calculated_symmetry_pairs[i].first;
		int calculated_index = m->calculated_symmetry_pairs[i].second;
		float geo_cost = Metric_get_geodesic_cost(m, index, calculated_index, false);
		 
		if (geo_cost < threshold)
		{
			passed_point_count++;
		}
	}

	percentage = passed_point_count / (float)N;
	return percentage;
}

static float s_gaussian_swipe_distance;
static int s_gaussian_point_no;
void Metric_set_gaussian(TrilateralMesh* m , int gaussian_point , float gaussian_dist )
{
	s_gaussian_point_no = gaussian_point;
	s_gaussian_swipe_distance = gaussian_dist;
}
float Metric_get_gaussian_dist(TrilateralMesh* m)
{
	return s_gaussian_swipe_distance;
}
int Metric_get_gaussian_point_no(TrilateralMesh* m)
{
	return s_gaussian_point_no;
}
static std::string getCurrentDate() {
	time_t now = time(0);
	tm* ltm = localtime(&now);
	char date[20];
	sprintf(date, "%04d-%02d-%02d-%02d-%02d-%02d", 1900 + ltm->tm_year, 1 + ltm->tm_mon, ltm->tm_mday , ltm->tm_hour , ltm->tm_min , ltm->tm_sec );
	return std::string(date);
}
void Metric_write_to_file(TrilateralMesh* m, const std::string& file_name)
{
	std::string concat_name = file_name;
	//just in case 
	read_symmetry_format((char*)"../../Trilateral/Mesh/off/sym.txt", m);
	
	std::ofstream out_file;
	// Specify the file name
	// attach date to file
	std::time_t currentTime = std::time(nullptr);
	// Convert to local time
	std::tm* localTime = std::localtime(&currentTime);
	// Convert to a readable format
	concat_name = concat_name + getCurrentDate() + ".txt";
	// Open the file
	out_file.open(concat_name);
	if (!out_file) {
		std::cerr << "Error opening file: " << std::endl;
		exit(1);
	}
	// 1 - write geodesic cost
	out_file << "Normalize geodesic cost === " << std::to_string(Metric_get_geodesic_cost(m)) << std::endl;
	out_file << "correspondance rate === " << std::to_string(Metric_get_correspondance_rate(m)) << std::endl;
	out_file << " gaussian point size == " << std::to_string(Metric_get_gaussian_point_no(m)) << std::endl;
	out_file << " gaussian sweep distance == " << std::to_string(Metric_get_gaussian_dist(m)) << std::endl;
	out_file.close();
}