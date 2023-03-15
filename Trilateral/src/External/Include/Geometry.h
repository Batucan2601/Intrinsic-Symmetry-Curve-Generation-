#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <iostream>
#pragma once
struct parametric
{
	glm::vec3 t; //t[0] x t[1] y t[2] z 
	glm::vec3 d;
};
static parametric create_line(glm::vec3 p1, glm::vec3 p2)
{
	parametric p;
	p.t = p1 - p2;
	p.d = p2;
	return p;
}
parametric parametric_line_2_planes(glm::vec4 plane1 , glm::vec4 plane2 ) //ax by cz = d
{
	glm::vec3 plane1_valx = glm::vec3(-1 * plane1[1] / plane1[0], -1 * plane1[2] / plane1[0], plane1[3] / plane1[0]); // x = (-by-cz+d ) / a   
	glm::vec3 plane2_valx = glm::vec3(-1 * plane2[1] / plane1[0], -1 * plane2[2] / plane1[0], plane2[3] / plane1[0]);
	

	// do x1 = x2 
	glm::vec2 valz = glm::vec2(-1 * plane1_valx[0] + plane2_valx[0], -1* plane1_valx[2] + plane2_valx[2]) / (plane1_valx[1] - plane2_valx[1]); // z = cy + bd 

	//get the parametric eqaution of x now
	glm::vec2 valx = glm::vec2(plane1_valx[0] + plane1_valx[1] * (valz[0]) , plane1_valx[2] + plane1_valx[1] * (valz[1])) ; // y , d according to plane1
	//now we now every variables in terms of y which is t now
	parametric param;
	param.t = glm::vec3( valx[0], 1 , valz[0]);
	param.d = glm::vec3(valx[1], 0, valz[1]);
	std::cout << " param " << std::endl;
	std::cout << param.t[0] << "  " << param.t[1] << "  " << param.t[2] << std::endl;
	std::cout << param.d[0] << "  " << param.d[1] << "  " << param.d[2] << std::endl;
	return param; 
}

float distance_between_point_and_line(glm::vec3 p1 , parametric &line )
{
	parametric perpendecular_line; 
	/*perpendecular_line.t[0] = line.t[0];
	perpendecular_line.t[1] = line.t[1];
	perpendecular_line.t[2] = line.t[2];*/

	perpendecular_line.d[0] = line.d[0] - p1[0];
	perpendecular_line.d[1] = line.d[1] - p1[1];
	perpendecular_line.d[2] = line.d[2] - p1[2];

	float cross_length = glm::length(glm::cross(perpendecular_line.d, p1));
	return  cross_length / glm::length(line.t);
	////do the dot product
	//parametric result; 
	//result.d[0] = perpendecular_line.d[0] * perpendecular_line.t[0];
	//result.d[1] = perpendecular_line.d[1] * perpendecular_line.t[1];
	//result.d[2] = perpendecular_line.d[2] * perpendecular_line.t[2];

	//result.t[0] = perpendecular_line.t[0] * perpendecular_line.t[0];
	//result.t[1] = perpendecular_line.t[1] * perpendecular_line.t[1];
	//result.t[2] = perpendecular_line.t[2] * perpendecular_line.t[2];

	//float new_t = -1 * (result.d[0] + result.d[1] + result.d[2]) / (result.t[0] + result.t[1] + result.t[2]);
	////problem not solved t = res 
	//glm::vec3 new_point(result.d[0] + new_t * result.t[0], result.d[1] + new_t * result.t[1], result.d[2] + new_t * result.t[2]);
	//	//distance
	//float dist = glm::distance(p1, new_point);
	//return dist;
}
static void sort_glm(std::vector<glm::vec3>& points, int index) // index 0 x 1 y 2 z
{
	//insertion sort so few numbers : ( 
	for (size_t i = 0; i < points.size(); i++)
	{
		float lowest = INFINITY;
		int lowest_index = 0;
		for (size_t j = i; j < points.size(); j++)
		{
			if (lowest > points[j][index])
			{
				lowest = points[j][index];
				lowest_index = j;
			}
		}
		//swap
		glm::vec3 temp = points[i];
		points[i] = points[lowest_index];
		points[lowest_index] = temp;
	}

}
static void polygonise(std::vector<glm::vec3> &points )
{
	//get the  leftmost and rigthmost  point 
	glm::vec3 rigthmost = points[0];
	glm::vec3 leftmost = points[points.size() - 1 ];
	std::vector<glm::vec3> up_points;
	std::vector<glm::vec3> down_points;

	parametric line = create_line(leftmost , rigthmost);
	
	// now push every point
	//push first
	up_points.push_back(points[0]);
	//push all others
	for (size_t i = 1; i < points.size() - 1; i++)
	{
		//debunk parametric 
		float t = (points[i].x - line.d[0] )/ line.t[0] ;
		float y = t * line.t[1] + line.d[1];
		if (y <= points[i].y)
		{
			up_points.push_back(points[i]);
		}
		else
		{
			down_points.push_back(points[i]);

		}
	}
	// push last
	up_points.push_back(points[points.size()-1 ]);
	//clear
	points.clear();
	//concat
	for (size_t i = 0; i < up_points.size(); i++)
	{
		points.push_back(up_points[i]);
	}
	for (size_t i = 0; i < down_points.size(); i++)
	{
		points.push_back(down_points[i]);
	}
	
}
static float triangle_area(glm::vec3 p1 , glm::vec3 p2 , glm::vec3 p3 )
{
	glm::vec3 edge1 =  p1 - p2;
	glm::vec3 edge2 =  p3 - p2; 

	return glm::length(glm::cross(edge1 , edge2) / 2.0f );
}
static float easy_triangulation(std::vector<glm::vec3>& points)
{
	// really noob triangulation
	// get the middle point and triangulate
	//not recommended for further applications
	glm::vec3 medium_point = glm::vec3(0.0f , 0.0f , 0.0f);
	for (size_t i = 0; i < points.size(); i++)
	{
		medium_point[0] += points[i][0];
		medium_point[1] += points[i][1];
		medium_point[2] += points[i][2];

	}
	medium_point /= points.size();
	
	float total_area = 0.0f;
	// start constructing triangles 
	for (size_t i = 0; i < points.size()-1; i++)
	{
		total_area += triangle_area(medium_point, points[i], points[i + 1]);
	}
	//the last one
	total_area += triangle_area(medium_point, points[points.size() - 1], points[0]);
	return total_area;
}

static float calculate_area(std::vector<glm::vec3>& points)
{
	// sort them with x 
	// write sort code yourself : ( 
	sort_glm(points, 0);
	// clockwise polygon starting from left 
	polygonise(points);
	return easy_triangulation(points);
}

static bool calculate_line_segment_intersection(glm::vec3 A , glm::vec3 B , glm::vec3 C , glm::vec3 D , glm::vec3 & intersection_point )
{
	float s = (((D.y - B.y) * (A.x - B.x)) - ((D.x - B.x) * (A.y - B.y))) / (((A.x - B.x) * (D.y - C.y)) - ((A.y - B.y) * (D.x - C.x)));
	float t = ( (D.x - B.x) - s * (D.x - C.x) )  / (A.x - B.x);
	float third_eq = t * (A.z - B.z) + s * (D.z - C.z);
	if (glm::abs(third_eq - (D.z - B.z) ) < 0.05f) // error payi
	{
		if (t > 0.0f && t < 1.0f && s < 1.0f && s > 0.0f)
		{
			intersection_point = A * t + (1 - t) * B;
			return true; 
		}
		return false; 
	}
	return false; 
}
//static glm::vec3 calculate_line_segment_intersection(parametric p, parametric q)
//{
//	glm::vec3 q_p = q.d - p.d;
//	glm::vec3 r_s = glm::cross(p.t, q.t);
//	if (r_s[0] == 0 && r_s[1] == 0 && r_s[2] == 0)
//	{
//		return glm::vec3(-1, -1, -1);
//	}
//	glm::vec3 res_t = glm::cross(q_p, q.t) / r_s;
//	if (res_t[0] == 0 && res_t[1] == 0 && res_t[2] == 0)
//	{
//		return glm::vec3(-1, -1, -1);
//	}
//	if ((res_t[0] == INFINITY || res_t[0] == -1 * INFINITY) || (res_t[1] == INFINITY || res_t[1] == -1 * INFINITY) || (res_t[2] == INFINITY || res_t[2] == -1 * INFINITY))
//	{
//		return glm::vec3(-1, -1, -1);
//	}
//	return glm::vec3(p.d[0] + p.t[0] * res_t[0] , p.d[1] + p.t[1] * res_t[1], p.d[2] + p.t[2] * res_t[2]);
//}


//origin is the central point
//should get the 4 points around the origin
std::vector<glm::vec3> get_elliptic_points(glm::vec4 plane, glm::vec3 origin, glm::vec3 point2, glm::vec2 slope, float radius)
{
	//create the vector
	std::vector<glm::vec3> vector;
	//an identity matrix to do calculation on
	glm::mat4 identity = glm::mat4(1);
	//get the 4 points
	glm::vec4 top_point = glm::vec4(origin.x, origin.y + radius, origin.z, 1.0f);
	glm::vec4 down_point = glm::vec4(origin.x, origin.y - radius, origin.z, 1.0f);
	glm::vec4 left_point = glm::vec4(origin.x + radius, origin.y, origin.z, 1.0f);
	glm::vec4 right_point = glm::vec4(origin.x - radius, origin.y, origin.z, 1.0f);
	//rotate all
	//identity = glm::rotate(identity, glm::radians(slope.x), glm::vec3(1.0f, 0.0f, 0.0f));
	//identity = glm::rotate(identity, glm::radians(180 - slope.y), glm::vec3(0.0f, 1.0f, 0.0f));
	top_point = identity * top_point;
	down_point = identity * down_point;
	left_point = identity * left_point;
	right_point = identity * right_point;

	glm::vec3 intersection_point_up = line_plane_intersection(plane, glm::vec3(top_point.x, top_point.y, top_point.z), glm::vec3(top_point.x + point2.x, top_point.y + point2.y, top_point.z + point2.z));
	glm::vec3 intersection_point_down = line_plane_intersection(plane, glm::vec3(down_point.x, down_point.y, down_point.z), glm::vec3(down_point.x + point2.x, down_point.y + point2.y, down_point.z + point2.z));
	glm::vec3 intersection_point_left = line_plane_intersection(plane, glm::vec3(left_point.x, left_point.y, left_point.z), glm::vec3(left_point.x + point2.x, left_point.y + point2.y, left_point.z + point2.z));
	glm::vec3 intersection_point_right = line_plane_intersection(plane, glm::vec3(right_point.x, right_point.y, right_point.z), glm::vec3(right_point.x + point2.x, right_point.y + point2.y, right_point.z + point2.z));

	vector.push_back(intersection_point_up);
	vector.push_back(intersection_point_down);
	vector.push_back(intersection_point_left);
	vector.push_back(intersection_point_right);
	// kolaya kac, ucgen kenarlari gibi dusun 
	return vector;
}

static std::vector<int> calculate_sensors(glm::vec3 p1 , glm::vec3 p2, glm::vec3 p3 , int sensor_no , std::vector<glm::vec3> & hit_points )
{
	std::vector<glm::vec3> sensor_vec; 
	float dist = p3[0] - p2[0];
	for (size_t i = 0; i < sensor_no; i++)
	{
		glm::vec3 new_p_left = glm::vec3(  p2.x + ((dist / sensor_no) * i )  , 0.0f , p2.z );
		glm::vec3 new_p_right = glm::vec3(p2.x + ((dist / sensor_no) * (i+1) ) , 0.0f, p2.z);
		sensor_vec.push_back(p1);
		sensor_vec.push_back(new_p_left);
		sensor_vec.push_back(new_p_right);
	}

	// get the points  
	std::vector<int> sensors; 
	for (size_t i = 0; i < hit_points.size(); i++)
	{
		for (size_t j = 0; j < sensor_vec.size() / 3 ; j++)
		{
			if (is_point_in_triangle(hit_points[i], sensor_vec[j * 3 + 0], sensor_vec[j * 3 + 1], sensor_vec[j * 3 + 2]))
			{
				sensors.push_back(j);
			}
		}
	}
	return sensors;
}
// intersections between triangle and the ellipse 
std::vector<glm::vec3>  find_intersections(std::vector<glm::vec3> ellipse,  glm::vec3 p1 , glm::vec3 p2 , glm::vec3 p3  , int sensor_no , int &biggest_sensor , int& smallest_sensor , float & total_area )

{
	//create 3 lines
	parametric edge_1 = create_line(p1 , p2 );
	parametric edge_2 = create_line(p1, p3);
	parametric edge_3 = create_line(p2, p3);

	//create 4 lines
	parametric ellipse_up_right = create_line(ellipse[0], ellipse[3]);
	parametric ellipse_down_right = create_line(ellipse[1], ellipse[3]);
	parametric ellipse_up_left = create_line(ellipse[0], ellipse[2]);
	parametric ellipse_down_left = create_line(ellipse[1], ellipse[2]);

	//get intersection points
	std::vector<glm::vec3> point_vec;
	glm::vec3 point = glm::vec3(0.0f , 0.0f , 0.0f) ; 
	if (calculate_line_segment_intersection(p1, p2, ellipse[0], ellipse[3], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p2, ellipse[1], ellipse[3], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p2, ellipse[0], ellipse[2], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p2, ellipse[1], ellipse[2], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p3, ellipse[0], ellipse[3], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p3, ellipse[1], ellipse[3], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p3, ellipse[0], ellipse[2], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p1, p3, ellipse[1], ellipse[2], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p3, p2, ellipse[0], ellipse[3], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p3, p2, ellipse[1], ellipse[3], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p3, p2, ellipse[0], ellipse[2], point))
	{
		point_vec.push_back(point);
	}
	if (calculate_line_segment_intersection(p3, p2, ellipse[1], ellipse[2], point))
	{
		point_vec.push_back(point);
	}
	// point in triangle
	if (is_point_in_triangle(ellipse[0], p1, p2, p3))
	{
		point_vec.push_back(ellipse[0]);
	}
	if (is_point_in_triangle(ellipse[1], p1, p2, p3))
	{
		point_vec.push_back(ellipse[1]);
	}
	if (is_point_in_triangle(ellipse[2], p1, p2, p3))
	{
		point_vec.push_back(ellipse[2]);
	}
	if (is_point_in_triangle(ellipse[3], p1, p2, p3))
	{
		point_vec.push_back(ellipse[3]);
	}

	std::vector<int> sensors = calculate_sensors(p1, p2, p3, sensor_no , point_vec);
	//get the smallest and biggest 
	
	for (size_t i = 0; i < sensors.size(); i++)
	{
		if (biggest_sensor < sensors[i])
		{
			biggest_sensor = sensors[i];
		}
		if (smallest_sensor > sensors[i])
		{
			smallest_sensor = sensors[i];
		}
	}
	;
	//calculate area 
	if (smallest_sensor != 10000)
	{
		total_area = calculate_area(point_vec);
	}


	return point_vec;
}

