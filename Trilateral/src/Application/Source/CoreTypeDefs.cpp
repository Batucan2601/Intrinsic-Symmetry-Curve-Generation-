#include "../Include/CoreTypeDefs.h"
#include <src/Application/Include/TrilateralMap.h>

Mesh generate_mesh_from_plane( Plane* plane, glm::vec3 * m)
{
	// a(x - x0) + b(y - y0) + c(z - z0) = 0
	// only increase in y and z 
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;
	glm::vec3 p4;

	//increase x by two  increase y by two
	p1.x = m->x + 200.0f;
	p1.y = m->y + 200.0f;
	p1.z = -(plane->normal.x * (p1.x - m->x) + plane->normal.y * (p1.y - m->y)) / plane->normal.z + m->z;

	// increase x by two  decrease y by two 
	p2.x = m->x + 200.0f;
	p2.y = m->y - 200.0f;
	p2.z = -(plane->normal.x * (p2.x - m->x) + plane->normal.y * (p2.y - m->y)) / plane->normal.z + m->z;

	// decrease x by two  decrease y by two 
	p3.x = m->x - 200.0f ;
	p3.y = m->y - 200.0f;
	p3.z = -(plane->normal.x * (p3.x - m->x) + plane->normal.y * (p3.y - m->y)) / plane->normal.z + m->z;

	// decrease x by two increase y by two 
	p4.x = m->x - 200.0f;
	p4.y = m->y + 200.0f;
	p4.z = -(plane->normal.x * (p4.x - m->x) + plane->normal.y * (p4.y - m->y)) / plane->normal.z + m->z;

	Mesh plane_mesh(&p1,&p2,&p3,&p4); 
	return plane_mesh;
}

float get_point_status_from_plane(Plane* plane, glm::vec3* point)
{
	// negative if behind
	// positive if front
	// 0 if inside

	//plane equation 
	float result = plane->normal.x * (plane->point.x - point->x) + plane->normal.y * (plane->point.y - point->y) + plane->normal.z * (plane->point.z - point->z );

	return result;

}

glm::vec3 project_point_to_plane(Plane* plane,  glm::vec3* point)
{
	glm::vec3 v = *point - plane->point;
	float dist = glm::dot(v, plane->normal);
	return *point - dist * plane->normal; 
}
glm::vec3 symmetry_point_from_plane(Plane* plane, glm::vec3* point)
{
	glm::vec3 projected_point = project_point_to_plane(plane, point);
	//halfway there
	glm::vec3 sym_point = projected_point + (projected_point - (*point));

	return sym_point;
}
Plane generate_plane_from_two_vectors(const glm::vec3& vec1, const glm::vec3& vec2)
{
	Plane p; 
	p.normal = glm::normalize(glm::cross(vec1, vec2));
	p.point = glm::normalize(vec1); 

	return p; 
}
Plane generate_plane_from_three_points(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3)
{
	Plane p;

	glm::vec3 v1 = p2 - p1;
	glm::vec3 v2 = p3 - p1;
	glm::vec3 normal = glm::normalize(glm::cross(v1, v2));
	//float d = -glm::dot(normal, p1);
	p.normal = normal; 
	p.point = p1; 
	return p;
}

// Check if two line segments intersect
static bool is_line_segments_intersect(const glm::vec3& p1, const glm::vec3&  q1, const glm::vec3&  p2, const glm::vec3&  q2) {
	glm::vec3 v1(q1.x - p1.x, q1.y - p1.y, q1.z - p1.z);
	glm::vec3 v2(q2.x - p2.x, q2.y - p2.y, q2.z - p2.z);

	glm::vec3 v = glm::vec3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);

	double dot1 = v.x * v1.x + v.y * v1.y + v.z * v1.z;
	double dot2 = v.x * v2.x + v.y * v2.y + v.z * v2.z;

	double len1 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
	double len2 = v2.x * v2.x + v2.y * v2.y + v2.z * v2.z;

	double det = len1 * len2 - dot1 * dot2;

	if (fabs(det) < 1e-6) {
		// Line segments are parallel or collinear
		return false;
	}

	double s = (dot2 * len1 - dot1 * dot2) / det;
	double t = (dot1 * len2 - dot1 * dot2) / det;

	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
		// Line segments intersect within their bounds
		return true;
	}

	return false;
}

// Calculate the intersection point of two line segments
static glm::vec3 find_intersection_point(const glm::vec3& p1, const glm::vec3&  q1, const glm::vec3&  p2, const glm::vec3& q2) {
	glm::vec3  v1(q1.x - p1.x, q1.y - p1.y, q1.z - p1.z);
	glm::vec3  v2(q2.x - p2.x, q2.y - p2.y, q2.z - p2.z);

	glm::vec3  v = glm::vec3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);

	double dot1 = v.x * v1.x + v.y * v1.y + v.z * v1.z;
	double dot2 = v.x * v2.x + v.y * v2.y + v.z * v2.z;

	double len1 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
	double len2 = v2.x * v2.x + v2.y * v2.y + v2.z * v2.z;

	double det = len1 * len2 - dot1 * dot2;

	double s = (dot2 * len1 - dot1 * dot2) / det;

	double x = p1.x + s * v1.x;
	double y = p1.y + s * v1.y;
	double z = p1.z + s * v1.z;

	return glm::vec3(x, y, z);
}

//check if they intersect given that they are in same plane 
bool is_triangles_intersect(const glm::vec3& p1_1, const glm::vec3& p1_2, const glm::vec3& p1_3, const glm::vec3& p2_1, const glm::vec3& p2_2, const glm::vec3& p2_3)
{
	glm::vec3 vec1 = p1_2 - p1_1;
	glm::vec3 vec2 = p1_3 - p1_1;

	glm::vec3 cross1 = glm::cross(vec1, vec2);
	float area = glm::length(cross1) / 2;

	//chekc other 3 points
	
	//1 - select p2_1 
	float triangle_area_p2_1 = 0;
	vec1 = p1_1 - p2_1;
	vec2 = p1_2 - p2_1;
	triangle_area_p2_1  += glm::length(glm::cross(vec1 , vec2)) / 2;
	vec1 = p1_1 - p2_1;
	vec2 = p1_3 - p2_1;
	triangle_area_p2_1 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p1_2 - p2_1;
	vec2 = p1_3 - p2_1;
	triangle_area_p2_1 += glm::length(glm::cross(vec1, vec2)) / 2;
	
	if (abs(triangle_area_p2_1 - area) < 1e-5)
	{
		//inside
		return true; 
	}
	//2 - select p2_2 
	float triangle_area_p2_2 = 0;
	vec1 = p1_1 - p2_2;
	vec2 = p1_2 - p2_2;
	triangle_area_p2_2 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p1_1 - p2_2;
	vec2 = p1_3 - p2_2;
	triangle_area_p2_2 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p1_2 - p2_2;
	vec2 = p1_3 - p2_2;
	triangle_area_p2_2 += glm::length(glm::cross(vec1, vec2)) / 2;

	if (abs(triangle_area_p2_2 - area) < 1e-5)
	{
		//inside
		return true;
	}

	//1 - select p2_3 
	float triangle_area_p2_3 = 0;
	vec1 = p1_1 - p2_3;
	vec2 = p1_2 - p2_3;
	triangle_area_p2_3 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p1_1 - p2_3;
	vec2 = p1_3 - p2_3;
	triangle_area_p2_3 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p1_2 - p2_3;
	vec2 = p1_3 - p2_3;
	triangle_area_p2_3 += glm::length(glm::cross(vec1, vec2)) / 2;

	if (abs(triangle_area_p2_3 - area) < 1e-5)
	{
		//inside
		return true;
	}


	//check the opposite, but for one point 
	//1 - select p2_3 
	float triangle_area_p1_1 = 0;
	vec1 = p2_1 - p1_1;
	vec2 = p2_2 - p1_1;
	triangle_area_p1_1 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p2_1 - p1_1;
	vec2 = p2_3 - p1_1;
	triangle_area_p1_1 += glm::length(glm::cross(vec1, vec2)) / 2;
	vec1 = p2_2 - p1_1;
	vec2 = p2_3 - p1_1;
	triangle_area_p1_1 += glm::length(glm::cross(vec1, vec2)) / 2;

	vec2 = p2_3 - p2_1;
	vec1 = p2_2 - p2_1;

	cross1 = glm::cross(vec1, vec2);
	area = glm::length(cross1) / 2;

	if (abs(triangle_area_p1_1 - area) < 1e-5)
	{
		//inside
		return true;
	}
	
	
	return false; 
}

float area_of_triangle_intersection(const glm::vec3& p1_1, const glm::vec3& p1_2, const glm::vec3& p1_3, const glm::vec3& p2_1, const glm::vec3& p2_2, const glm::vec3& p2_3)
{
	std::vector<glm::vec3> convex_point_list;
	// try 3 edge of triangle 1 for 3 other edge of triangle 2
	bool is_intersect = false;
	int hit_no = 0; 
	// 1 - p1_1-p1_2  
	// 1.1 - p1_1-p1_2  vs p2_1 p2_2
	if (is_line_segments_intersect(p1_1, p1_2, p2_1, p2_2))
	{
		convex_point_list.push_back(find_intersection_point(p1_1, p1_2, p2_1, p2_2));
	}
	//1.2- p1_1 p1_2 vs p2_1 p2_3 
	if (is_line_segments_intersect(p1_1, p1_2, p2_1, p2_3))
	{
		convex_point_list.push_back(find_intersection_point(p1_1, p1_2, p2_1, p2_3));
	}
	//1.2- p1_1 p1_2 vs p2_2 p2_3 
	if (is_line_segments_intersect(p1_1, p1_2, p2_2, p2_3))
	{
		convex_point_list.push_back(find_intersection_point(p1_1, p1_2, p2_2, p2_3));
	}

	// 2 - p1_1-p1_3  
	// 1.1 - p1_1-p1_3  vs p2_1 p2_2
	if (is_line_segments_intersect(p1_1, p1_3, p2_1, p2_2))
	{
		convex_point_list.push_back(find_intersection_point(p1_1, p1_3, p2_1, p2_2));
	}
	//1.2- p1_1 p1_2 vs p2_1 p2_3 
	if (is_line_segments_intersect(p1_1, p1_3, p2_1, p2_3))
	{
		convex_point_list.push_back(find_intersection_point(p1_1, p1_3, p2_1, p2_3));
	}
	//1.2- p1_1 p1_2 vs p2_2 p2_3 
	if (is_line_segments_intersect(p1_1, p1_2, p2_2, p2_3))
	{
		convex_point_list.push_back(find_intersection_point(p1_1, p1_3, p2_2, p2_3));
	}

	// 1 - p1_2-p1_3  
	// 1.1 - p1_1-p1_2  vs p2_1 p2_2
	if (is_line_segments_intersect(p1_2, p1_3, p2_1, p2_2))
	{
		convex_point_list.push_back(find_intersection_point(p1_2, p1_3, p2_1, p2_2));
	}
	//1.2- p1_1 p1_2 vs p2_1 p2_3 
	if (is_line_segments_intersect(p1_2, p1_3, p2_1, p2_3))
	{
		convex_point_list.push_back(find_intersection_point(p1_2, p1_3, p2_1, p2_3));
	}
	//1.2- p1_1 p1_2 vs p2_2 p2_3 
	if (is_line_segments_intersect(p1_2, p1_3, p2_2, p2_3))
	{
		convex_point_list.push_back(find_intersection_point(p1_2, p1_3, p2_2, p2_3));
	}

	// cases

	// no intersection
	//either inside or outside 
	if (convex_point_list.size() == 0)
	{
		//create 
	}
	return 1;
}

void get_coefficients_from_plane(const Plane& plane, float& A, float& B, float& C, float& D)
{
	A = plane.normal.x;
	B = plane.normal.y;
	C = plane.normal.z;

	D = -1 * (plane.normal.x * plane.point.x + plane.normal.y * plane.point.y + plane.normal.z * plane.point.z);
}

std::vector<int> getNumberFromString(std::string s)
{
	std::stringstream str_strm;
	str_strm << s; //convert the string s into stringstream
	std::string temp_str;
	std::vector<int> num_vec; 
	int temp_int;
	while (!str_strm.eof()) {
		str_strm >> temp_str; //take words into temp_str one by one
		if (std::stringstream(temp_str) >> temp_int) { //try to convert string to int
			num_vec.push_back(temp_int);
		}
		temp_str = ""; //clear temp string
	}
	return num_vec;
}

Eigen::VectorXd stdVectorToEigenVectorXd(const std::vector<float>& std_vec) {
	Eigen::VectorXd eigen_vec(std_vec.size());
	for (size_t i = 0; i < std_vec.size(); ++i) {
		eigen_vec[i] = std_vec[i];
	}
	return eigen_vec;
}

float distancePointToLine(const glm::vec3& point, const glm::vec3& linePoint1, const glm::vec3& linePoint2) 
{
	// Calculate the line direction vector
	glm::vec3 lineDir = linePoint2 - linePoint1;

	// Calculate the vector from linePoint1 to the given point
	glm::vec3 pointVec = point - linePoint1;

	// Project pointVec onto lineDir to find the projection length
	float projectionLength = glm::dot(pointVec, lineDir) /  std::pow(glm::length(lineDir) ,2 );

	// Calculate the closest point on the line to the given point
	glm::vec3 closestPoint = linePoint1 + projectionLength * lineDir;

	// Calculate the distance from the closest point to the given point
	float distance = glm::length(point - closestPoint);

	return distance;
}


float compute_triangle_area(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 vec1 = p2 - p1;
	glm::vec3 vec2 = p3 - p1;

	glm::vec3 cross1 = glm::cross(vec1, vec2);
	float length = glm::length(cross1) / 2;
	return length;
}


float permutation_return_smallest_dif(Eigen::VectorXf vec1, Eigen::VectorXf vec2, int N)
{
	//generate permuation
	std::vector<unsigned int> permutation_vector;
	for (size_t i = 0; i < N; i++)
	{
		permutation_vector.push_back(i);
	}
	//all permutations for indices
	std::vector<std::vector<unsigned int>> all_permutations;
	std::vector<Eigen::VectorXf> vec1_list;
	std::vector<Eigen::VectorXf> vec2_list;

	do {
		all_permutations.push_back(permutation_vector);
	} while (std::next_permutation(permutation_vector.begin(), permutation_vector.end()));

	for (size_t i = 0; i < all_permutations.size(); i++)
	{
		Eigen::VectorXf vec1_p(N);
		Eigen::VectorXf vec2_p(N);
		for (size_t j = 0; j < N; j++)
		{
			vec1_p(j) = vec1(all_permutations[i][j]);
			vec2_p(j) = vec2(all_permutations[i][j]);
		}
		vec1_list.push_back(vec1_p);
		vec2_list.push_back(vec2_p);
	}
	int minimum_index = -1;
	float minimum_value = INFINITY;
	for (size_t i = 0; i < all_permutations.size(); i++)
	{
		for (size_t j = 0; j < all_permutations.size(); j++)
		{
			float diff = (vec1_list[i] - vec2_list[j]).norm();
			if (minimum_value > diff)
			{
				minimum_value = diff;
			}
		}
	}

	return minimum_value;
}
Plane generate_plane_from_formula(const float& A, const float& B, const float& C, const float& D)
{
	Plane p; 
	p.normal = glm::vec3(A, B, C);
	float x = 0;
	float y = 0; //C*z + D = 0 
	float z;
	if (C != 0)
	{
		z = -D / C; // z = -D / C 
	}

	p.point = glm::vec3(x, y, z);

	return p; 
}