#pragma once
#include <random>

glm::vec3 generate_random_point( )
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 10.0);
	float k = dist(rd);
	//find x 
	std::random_device generator_x;
	std::mt19937 mt_x(generator_x());
	std::uniform_real_distribution<double> distribution_x(min_x_range, max_x_range);
	float x = distribution_x(generator_x);
	//find y
	std::random_device generator_y;
	std::mt19937 mt_y(generator_y());
	std::uniform_real_distribution<double> distribution_y(min_y_range, max_y_range);
	float y = distribution_y(generator_y);
	//find z
	std::random_device generator_z;
	std::mt19937 mt_z(generator_z());
	std::uniform_real_distribution<double> distribution_z(min_z_range, max_z_range);
	float z = distribution_z(generator_z);
	
	glm::vec3 vec = glm::vec3(x, y, z);
	return vec; 
}

void generate_direction(glm::vec3 & random_point )
{
	float y_max = global_variables["arac yuksekligi"];
	float y_min = 0.0f;
	float x_max = global_variables["arac uzunlugu"];
	float x_min = 0.0f;

	glm::vec3 point_y_max, point_y_min, point_x_max, point_x_min; 
	point_y_max = glm::vec3(random_point.x, y_max, 0.0f);
	point_y_min = glm::vec3(random_point.x, y_min, 0.0f);
	point_x_max = glm::vec3(x_max, random_point.y, 0.0f);
	point_x_min = glm::vec3(x_min, random_point.y, 0.0f);

	float tan_y_max = (y_max - random_point.y) / (0.0f - random_point.z );
	float tan_y_min = (random_point.y - y_min ) / (random_point.z - 0.0f );

	float tan_x_max = (x_max - random_point.x) / (0.0f - random_point.z );
	float tan_x_min = (random_point.x - x_min) / (random_point.z - 0.0f );

	// get the arc
	float y_max_deg = glm::atan(tan_y_max);
	float y_min_deg = glm::atan(tan_y_min);
	float x_max_deg = glm::atan(tan_x_max);
	float x_min_deg = glm::atan(tan_x_min);

	//generate two numbers between
	//find x deg 
	std::default_random_engine generator_x;
	std::uniform_real_distribution<float> distribution_x(x_min_deg, x_max_deg);
	float x_deg = distribution_x(generator_x);

	//find y
	std::default_random_engine generator_y;
	std::uniform_real_distribution<float> distribution_y(y_min_deg, y_max_deg);
	float y_deg = distribution_y(generator_y);

	//setup rocket
	global_variables["roket gelis acisi x"] = x_deg; 
	global_variables["roket gelis acisi y"] = y_deg;
	global_variables["roket pozisyonu x"] = random_point.x;
	global_variables["roket pozisyonu y"] = random_point.y;
	global_variables["roket pozisyonu z"] = random_point.z;

}