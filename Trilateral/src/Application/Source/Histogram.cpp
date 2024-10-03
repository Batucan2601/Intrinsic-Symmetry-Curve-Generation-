#include "../Include/Histogram.h"
#include "../Include/Geodesic.h"
#include "../Include/ROI.h"
#define _USE_MATH_DEFINES
#include "math.h"


Histogram::Histogram(int size )
{
	for (size_t i = 0; i < size; i++)
	{
		this->histogram.push_back(0);
	}
}
Histogram::Histogram(std::vector<float> vec)
{
	this->histogram = vec; 
}
Histogram::Histogram()
{
}
//normalize to N  ( usually 1 ) 
void Histogram::normalize(float N)
{
	//normalize histogram.
	float histogram_sum = 0;
	for (size_t i = 0; i < this->histogram.size(); i++)
	{
		histogram_sum += this->histogram[i];
	}
	for (size_t i = 0; i < this->histogram.size(); i++)
	{
		this->histogram[i] /= histogram_sum;
		this->histogram[i] *= N;
	}

}
int Histogram::size()
{
	return  this->histogram.size();
}
float& Histogram::operator[](int index)
{
	return this->histogram[index];
}

float Histogram_L2Norm(const Histogram& h1, const Histogram& h2)
{
	if (h1.histogram.size() != h1.histogram.size())
	{
		return -1; 
	}
	Eigen::VectorXd  vector_h1 =  stdVectorToEigenVectorXd(h1.histogram);
	Eigen::VectorXd  vector_h2 =  stdVectorToEigenVectorXd(h2.histogram);

	return (vector_h1 - vector_h2).norm();
}

float Histogram_ChiSquareDistance(const Histogram& h1, const Histogram& h2)
{
	if (h1.histogram.size() != h2.histogram.size()) {
		std::cerr << "Error: Vectors must be of the same length." << std::endl;
		return -1;
	}
	float chiSquareDist = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		float numerator = std::pow(h1.histogram[i] - h2.histogram[i], 2);
		float denominator = h1.histogram[i] + h2.histogram[i];
		if (denominator != 0) {
			chiSquareDist += numerator / denominator;
		}
	}

	return chiSquareDist;
}

float Histogram_klDivergence(const Histogram& h1, const Histogram& h2 ) {
	// Check if both vectors are of the same length

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		sumP += h1.histogram[i];
		sumQ += h2.histogram[i];
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return -1;
	}

	double klDiv = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		if (h1.histogram[i] > 0.0 && h2.histogram[i] > 0.0) {
			klDiv += h1.histogram[i] * log(h1.histogram[i] / h2.histogram[i]);
		}
	}

	return klDiv;
}

float Histogram_earthMoversDistance(const Histogram& h1 , const Histogram& h2) {
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1;
	}

	// Check if the vectors are valid probability distributions
	double h1_sum = 0.0, h2_sum = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		h1_sum += h1.histogram[i];
		h2_sum += h2.histogram[i];
	}
	if (fabs(h1_sum - 1.0) > 1e-9 || fabs(h2_sum- 1.0) > 1e-9) {
		return -1;
	}

	// Compute the Earth Mover's Distance
	double cumulativeDifference = 0.0;
	double totalDistance = 0.0;

	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		cumulativeDifference += h1.histogram[i] - h2.histogram[i];
		totalDistance += fabs(cumulativeDifference);
	}

	return totalDistance;
}

float Histogram_bhattacharyyaDistance(const Histogram& h1, const Histogram& h2) 
{
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1; 
	}

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		sumP += h1.histogram[i];
		sumQ += h2.histogram[i];
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return  -1; 
	}

	// Calculate the Bhattacharyya coefficient
	double bc = 0.0;  // Bhattacharyya coefficient
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		if (h1.histogram[i] > 0.0 && h2.histogram[i] > 0.0) {
			bc += sqrt(h1.histogram[i] * h2.histogram[i]);
		}
	}

	// Return the Bhattacharyya distance
	return -log(bc);
}

float Histogram_jensenShannonDivergence(const Histogram& h1 , const Histogram& h2) 
{
	// Check if both vectors are of the same length
	if (h1.histogram.size() != h2.histogram.size()) {
		return -1; 
	}

	// Check if the vectors are valid probability distributions
	double sumP = 0.0, sumQ = 0.0;
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		sumP += h1.histogram[i];
		sumQ += h2.histogram[i];
	}
	if (fabs(sumP - 1.0) > 1e-9 || fabs(sumQ - 1.0) > 1e-9) {
		return -1; 
	}

	// Create the midpoint vector M = (P + Q) / 2
	std::vector<float> M(h1.histogram.size(), 0.0);
	for (size_t i = 0; i < h1.histogram.size(); ++i) {
		M[i] = 0.5 * (h1.histogram[i] + h2.histogram[i]);
	}

	// Compute the Jensen-Shannon Divergence
	Histogram M_h(M);
	float jsd = 0.5 * Histogram_klDivergence(h1,Histogram(M_h)) + 0.5 * Histogram_klDivergence(h2, Histogram(M_h));

	return jsd;
}

float Histogram_kolmogorovSmirnovTest(const Histogram& h1, const Histogram& h2) 
{
	if (h1.histogram.empty() || h2.histogram.empty()) {
		return -1;
	}

	// Combine both samples into a single sorted array
	std::vector<float> combined = h1.histogram;
	combined.insert(combined.end(),h2.histogram.begin(), h2.histogram.end());
	sort(combined.begin(), combined.end());

	// Calculate the CDF for each sample
	std::vector<double> cdf1(h1.histogram.size(), 0.0);
	std::vector<double> cdf2(h2.histogram.size(), 0.0);

	int n1 = h1.histogram.size();
	int n2 = h2.histogram.size();

	// Create step function counts for each sample
	int count1 = 0, count2 = 0;
	double dMax = 0.0;  // Maximum KS statistic

	// Calculate the cumulative step functions and the maximum deviation (KS statistic)
	for (const double& value : combined) {
		while (count1 < n1 && h1.histogram[count1] <= value) {
			count1++;
		}
		while (count2 < n2 && h2.histogram[count2] <= value) {
			count2++;
		}

		// Calculate the CDFs for each sample at the current point
		double cdf1Val = (double)count1 / n1;
		double cdf2Val = (double)count2 / n2;

		// Update the maximum deviation
		dMax = std::max(dMax, fabs(cdf1Val - cdf2Val));
	}

	return dMax;
}


Histogram histogram_roi_area_detailed(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited)
{
	//histogram to be returned 
	Histogram histogram(division_no);
	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, point_index1);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, point_index2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, point_index3);
	//find the maximum distance from is_visited and paths
	float max = -999;
	float min = 100000;
	std::vector<unsigned int> triangles_inside;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (is_visited[i] == INSIDE) // if vertex is visited 
		{
			for (size_t j = 0; j < path_1_2.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_1_2[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_1_2[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_1_2[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_1_2[j]]);
				}
			}
			for (size_t j = 0; j < path_1_3.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_1_3[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_1_3[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_1_3[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_1_3[j]]);
				}
			}
			for (size_t j = 0; j < path_2_3.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_2_3[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_2_3[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_2_3[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_2_3[j]]);
				}
			}
		}
	}
	//now recolor
	std::vector<glm::vec3> new_color_buffer;
	for (size_t i = 0; i < m->colors.size(); i++)
	{

		if (is_visited[i] == EDGE) //edge 
		{
			global_is_visited[i] = EDGE;
		}
		else if (is_visited[i] == OUTSIDE) //not visited 
		{
		}
		else if (is_visited[i] == INSIDE) // get the max distance
		{
			if (global_is_visited[i] != EDGE)
			{
				global_is_visited[i] = INSIDE;
			}
		}
	}
	float max_dist_inside = -1;
	float min_dist_inside = -1;
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		if (is_visited[index1] == INSIDE || is_visited[index2] == INSIDE || is_visited[index3] == INSIDE) //if any vertex is visited
		{
			triangles_inside.push_back(index1);
			triangles_inside.push_back(index2);
			triangles_inside.push_back(index3);
			if (distance_matrix_p1[index1] > max_dist_inside)
			{
				max_dist_inside = distance_matrix_p1[index1];
			}
			if (distance_matrix_p1[index2] > max_dist_inside)
			{
				max_dist_inside = distance_matrix_p1[index2];
			}
			if (distance_matrix_p1[index3] > max_dist_inside)
			{
				max_dist_inside = distance_matrix_p1[index3];
			}

			if (distance_matrix_p1[index1] < min_dist_inside)
			{
				min_dist_inside = distance_matrix_p1[index1];
			}
			if (distance_matrix_p1[index2] < min_dist_inside)
			{
				min_dist_inside = distance_matrix_p1[index2];
			}
			if (distance_matrix_p1[index3] < min_dist_inside)
			{
				min_dist_inside = distance_matrix_p1[index3];
			}
		}
	}
	float step = max_dist_inside / division_no;
	for (size_t i = 0; i < triangles_inside.size(); i += 3)
	{
		int i1 = triangles_inside[i];
		int i2 = triangles_inside[i + 1];
		int i3 = triangles_inside[i + 2];
		glm::vec3 p1 = m->vertices[i1];
		glm::vec3 p2 = m->vertices[i2];
		glm::vec3 p3 = m->vertices[i3];
		float i1_dist = distance_matrix_p1[i1];
		float i2_dist = distance_matrix_p1[i2];
		float i3_dist = distance_matrix_p1[i3];

		int step_no_i1 = i1_dist / step; // floor 
		int step_no_i2 = i2_dist / step; // floor 
		int step_no_i3 = i3_dist / step; // floor
		if (step_no_i1 == division_no)
		{
			step_no_i1--;
		}
		if (step_no_i2 == division_no)
		{
			step_no_i2--;
		}
		if (step_no_i3 == division_no)
		{
			step_no_i3--;
		}
		float triangle_area = compute_triangle_area(p1, p2, p3);
		if (step_no_i1 == step_no_i2 && step_no_i1 == step_no_i3)
		{
			histogram[step_no_i1] += triangle_area;
		}
		//one of them is in other step 
		else /*if (step_no_i1 == step_no_i2 && step_no_i1 != step_no_i3 &&
			 step_no_i1 == step_no_i3 && step_no_i1 != step_no_i2 &&
			 step_no_i2 == step_no_i3 && step_no_i2 != step_no_i3)*/
		{
			std::vector<std::pair<int, int>> steps;
			steps.push_back(std::pair<int, int>(step_no_i1, i1));
			steps.push_back(std::pair<int, int>(step_no_i2, i2));
			steps.push_back(std::pair<int, int>(step_no_i3, i3));
			CoreType_sort_by_value(steps);
			//they are now in ascending order
			// 1.1.1 small arc
			float dist_step_0 = distance_matrix_p1[steps[0].second];
			//find the distance to other hist
			float small_r = (step * (steps[0].first + 1)) - dist_step_0;
			glm::vec3 edge1 = m->vertices[steps[1].second] - m->vertices[steps[0].second];
			glm::vec3 edge2 = m->vertices[steps[2].second] - m->vertices[steps[0].second];
			float cosine = glm::dot(edge1, edge2) / (glm::length(edge1) * glm::length(edge2));
			float small_radian = acos(cosine); //radian

			float area_arc_small = M_PI * small_r * small_r * small_radian / (2 * M_PI);

			histogram[steps[0].first] += area_arc_small;


			// 1.1.2 closest arc
			//a miscalculation but think it like also an arc although it is a reverse arc
			float dist_step_2 = distance_matrix_p1[steps[2].second];
			float big_r = dist_step_2 - (steps[2].first * step);
			edge1 = m->vertices[steps[1].second] - m->vertices[steps[2].second];
			edge2 = m->vertices[steps[0].second] - m->vertices[steps[2].second];
			cosine = glm::dot(edge1, edge2) / (glm::length(edge1) * glm::length(edge2));
			small_radian = acos(cosine); //radian
			float area_arc_big = M_PI * big_r * big_r * small_radian / (2 * M_PI);
			histogram[steps[2].first] += area_arc_big;


			float area_left = triangle_area - (area_arc_small + area_arc_big);

			int steps_left = steps[2].first - (steps[0].first + 1);
			if (steps_left == 0)
			{
				histogram[steps[0].first] += area_left;
				histogram[steps[1].first] += area_left;
			}
			else
			{
				int total_fraction = steps_left * steps_left;
				for (size_t step_no = steps[0].first + 1; step_no < steps[2].first; step_no++)
				{
					int fraction_no = step_no - (steps[0].first);
					float fraction_percentage = (pow(fraction_no, 2) - pow(fraction_no - 1, 2)) / fraction_no;

					histogram[step_no] += fraction_percentage * area_left;
				}
			}



		}

	}

	histogram.normalize(1);
	
	return histogram;

}



Histogram  Histogram_triangle_area(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, int division_no,
	std::vector<int> is_visited, std::vector<int>& global_is_visited)
{
	Histogram histogram(division_no);
	std::vector<int> path_1_2 = Geodesic_between_two_points(*m, point_index1, point_index2);
	std::vector<int> path_1_3 = Geodesic_between_two_points(*m, point_index1, point_index3);
	std::vector<int> path_2_3 = Geodesic_between_two_points(*m, point_index2, point_index3);
	std::vector<float> distance_matrix_p1 = Geodesic_dijkstra(*m, point_index1);
	std::vector<float> distance_matrix_p2 = Geodesic_dijkstra(*m, point_index2);
	std::vector<float> distance_matrix_p3 = Geodesic_dijkstra(*m, point_index3);
	//find the maximum distance from is_visited and paths
	float max = -999;
	float min = 100000;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (is_visited[i] == INSIDE) // if vertex is visited 
		{
			for (size_t j = 0; j < path_1_2.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_1_2[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_1_2[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_1_2[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_1_2[j]]);
				}
			}
			for (size_t j = 0; j < path_1_3.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_1_3[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_1_3[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_1_3[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_1_3[j]]);
				}
			}
			for (size_t j = 0; j < path_2_3.size(); j++)
			{
				if (glm::length(m->vertices[i] - m->vertices[path_2_3[j]]) > max)
				{
					max = glm::length(m->vertices[i] - m->vertices[path_2_3[j]]);
				}
				if (glm::length(m->vertices[i] - m->vertices[path_2_3[j]]) < min)
				{
					min = glm::length(m->vertices[i] - m->vertices[path_2_3[j]]);
				}
			}
		}
	}
	float step = max / division_no;

	//now recolor
	std::vector<glm::vec3> new_color_buffer;
	for (size_t i = 0; i < m->colors.size(); i++)
	{

		if (is_visited[i] == EDGE) //edge 
		{
			//new_color_buffer.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
			//m->colors[i] = glm::vec3(1.0f, 0.0f, 0.0f);
			global_is_visited[i] = EDGE;
		}
		else if (is_visited[i] == OUTSIDE) //not visited 
		{
			//new_color_buffer.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			//m->colors[i] = glm::vec3(0.0f, 0.0f, 0.0f);

		}
		else if (is_visited[i] == INSIDE) // get the max distance
		{
			//global_is_visited[i] = INSIDE;

			if (global_is_visited[i] != EDGE)
			{
				global_is_visited[i] = INSIDE;
			}
		}
	}

	//m->colors = new_color_buffer;

	//getting the numbers
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		if (is_visited[m->triangles[i]] != 0 || is_visited[m->triangles[i + 1]] != 0 || is_visited[m->triangles[i + 2]] != 0) //if any vertex is visited
		{
			glm::vec3 p1 = m->vertices[m->triangles[i]];
			glm::vec3 p2 = m->vertices[m->triangles[i + 1]];
			glm::vec3 p3 = m->vertices[m->triangles[i + 2]];

			float min_dist_from_p1 = std::min(distance_matrix_p1[m->triangles[i]], distance_matrix_p2[m->triangles[i]]);
			min_dist_from_p1 = std::min(min_dist_from_p1, distance_matrix_p3[m->triangles[i]]);

			float min_dist_from_p2 = std::min(distance_matrix_p1[m->triangles[i + 1]], distance_matrix_p2[m->triangles[i + 1]]);
			min_dist_from_p2 = std::min(min_dist_from_p2, distance_matrix_p3[m->triangles[i + 1]]);

			float min_dist_from_p3 = std::min(distance_matrix_p1[m->triangles[i + 2]], distance_matrix_p2[m->triangles[i + 2]]);
			min_dist_from_p3 = std::min(min_dist_from_p3, distance_matrix_p3[m->triangles[i + 2]]);

			float area = compute_triangle_area(p1, p2, p3);

			int hist_no_p1 = min_dist_from_p1 / step;
			int hist_no_p2 = min_dist_from_p2 / step;
			int hist_no_p3 = min_dist_from_p3 / step;

			if (hist_no_p1 > division_no)
			{
				hist_no_p1 = division_no;
			}
			if (hist_no_p2 > division_no)
			{
				hist_no_p2 = division_no;
			}
			if (hist_no_p3 > division_no)
			{
				hist_no_p3 = division_no;
			}
			histogram[hist_no_p1] += area;
			histogram[hist_no_p2] += area;
			histogram[hist_no_p3] += area;
		}
	}

	//normalize histogram.
	float histogram_sum = 0;
	for (size_t i = 0; i < histogram.size(); i++)
	{
		histogram_sum += histogram[i];
	}
	for (size_t i = 0; i < histogram.size(); i++)
	{
		histogram[i] /= histogram_sum;
	}
	return histogram;
}