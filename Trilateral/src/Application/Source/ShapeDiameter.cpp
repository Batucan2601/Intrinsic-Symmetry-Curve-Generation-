#include "../Include/ShapeDiameter.h"
#include "../Include/Ray.h"
#include "../Include/CoreTypeDefs.h"
#include "glm/glm.hpp"
#define _USE_MATH_DEFINES //for pi 
#include <math.h>
#include <random>
#include <glm/ext/quaternion_trigonometric.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

static void calculateRightUpVectors(const glm::vec3& normal, const glm::vec3& pointOnPlane, glm::vec3& right, glm::vec3& up) {
	// Assume the plane's "forward" vector is the normal
	glm::vec3 forward = glm::normalize(normal);

	// Choose an arbitrary "right" vector
	glm::vec3 arbitraryRight = glm::vec3(1.0f, 0.0f, 0.0f);

	// Calculate the right vector by taking the cross product of the normal and the arbitrary right vector
	right = glm::normalize(glm::cross(forward, arbitraryRight));

	// Calculate the up vector by taking the cross product of the right vector and the normal
	up = glm::normalize(glm::cross(right, normal));
}

void ShapeDiameter_calculate(TrilateralMesh* mesh,  std::vector<unsigned int> indices , std::vector<float>& shape_diameter)
{

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> random_r(1, CONE_RADIUS * CONE_RADIUS_MULTIPLIER); // distribution in range [1, 6]

	// save each of the distances
	// for each endpoint
	for (size_t i = 0; i < indices.size(); i++)
	{
		//
		std::vector<float> ray_distances; 
		//generate ray
		TrilateralRay ray;
		ray.origin = mesh->vertices[indices[i]];
		ray.direction = 0.01f * (mesh->normals[indices[i]] * -1.0f) ; //inverse of normal
		glm::vec3 cone_center = ray.origin + ray.direction * (float)CONE_HEIGHT;
		for (size_t j = 0; j < NUMBER_OF_RAYS; j++)
		{
			glm::vec3 right;
			glm::vec3 left;
			float r = CONE_RADIUS * (float)std::sqrt(random_r(rng)) / (float)std::sqrt(CONE_RADIUS * CONE_RADIUS_MULTIPLIER);
			float theta = (float)std::sqrt(random_r(rng)) / (float)std::sqrt(CONE_RADIUS * CONE_RADIUS_MULTIPLIER) * 2 * M_PI;
			//generate right and up vector
			calculateRightUpVectors(  ray.direction , cone_center ,right , left  );

			TrilateralRay random_ray;
			random_ray.origin = ray.origin;
			glm::vec3 point_on_circle = cone_center + right * r * cos(theta) + left * r * sin(theta);
			random_ray.direction = glm::normalize(point_on_circle - random_ray.origin);
			glm::vec3 hit_point;

			int closest_triangle_index = -1;
			float closest_dist = INFINITY;
			for (size_t k = 0; k < mesh->triangles.size(); k+=3)
			{
				
				glm::vec3 p1 = mesh->vertices[mesh->triangles[k]];
				glm::vec3 p2 = mesh->vertices[mesh->triangles[k+1]];
				glm::vec3 p3 = mesh->vertices[mesh->triangles[k+2]];
				if (ray_triangle_intersection(random_ray, p1, p2, p3, hit_point))
				{
					float dist= glm::distance( ray.origin , hit_point);
					if (dist < closest_dist)
					{
						closest_triangle_index = k;
						closest_dist  = dist;

					}
				}

			}
			ray_distances.push_back(closest_dist);

			
		}
		//find the standart deviation 

		// 1 -calculate mean
		float mean = 0;
		for (size_t j = 0; j < ray_distances.size(); j++)
		{
			mean += ray_distances[j];
		}
		mean /= ray_distances.size();
		float variance = 0.0f;
		for (size_t j = 0; j < ray_distances.size(); j++)
		{
			variance += std::pow(ray_distances[j] - mean, 2);
		}
		variance /= ray_distances.size();

		float std_deviation = std::sqrt(variance);

		std::vector<float> standart_deviation_distances;
		for (size_t j = 0; j < ray_distances.size(); j++)
		{
			if (mean - std_deviation <= ray_distances[j] && mean + std_deviation >= ray_distances[j])
			{
				standart_deviation_distances.push_back(ray_distances[j]);
			}
		}
		
		float result = 0;
		// now weight them and generate a final ray
		for (size_t j = 0; j < standart_deviation_distances.size(); j++)
		{
			result += standart_deviation_distances[j];
		}
		result = result / standart_deviation_distances.size();

		shape_diameter.push_back(result);
	}

	//normalize sdf
	float max = 0; 
	for (size_t i = 0; i < shape_diameter.size(); i++)
	{
		if (max < shape_diameter[i])
		{
			max = shape_diameter[i];
		}
	}
	for (size_t i = 0; i < shape_diameter.size(); i++)
	{
		shape_diameter[i] /= max;
	}
}

float ShapeDiameter_calculate_simple(TrilateralMesh* mesh, unsigned int index)
{
	float min_distance = INFINITY;
	for (size_t i = 0; i < mesh->triangles.size(); i+=3)
	{
		int index1 = mesh->triangles[i];
		int index2 = mesh->triangles[i+1];
		int index3 = mesh->triangles[i+2];
		TrilateralRay ray; 
		TrilateralRay ray_reverse; 
		ray.origin = mesh->vertices[index];
		ray.direction = -mesh->normals[index];
		ray_reverse.origin  = mesh->vertices[index];
		ray_reverse.direction  = mesh->normals[index];
		glm::vec3 hit_point1;
		glm::vec3 hit_point2;
		bool is_hit = ray_triangle_intersection(ray, mesh->vertices[index1],
		mesh->vertices[index2], mesh->vertices[index3], hit_point1);
		bool is_hit_reverse = ray_triangle_intersection(ray_reverse, mesh->vertices[index1],
			mesh->vertices[index2], mesh->vertices[index3], hit_point2);
		if (is_hit)
		{
			float distance = glm::distance(ray.origin, hit_point1);
			if (distance < min_distance && distance != 0 )
			{
				min_distance = distance;
			}
		}
		else if (is_hit_reverse)
		{
			float distance = glm::distance(ray.origin, hit_point2);
			if (distance < min_distance && distance != 0)
			{
				min_distance = distance;
			}
		}
		
	}
	return min_distance;
}
float ShapeDiameter_calculate_simple_max_dif(TrilateralMesh* mesh, std::vector<unsigned int>& indices)
{
	std::vector<float> results;
	for (size_t i = 0; i < indices.size(); i++)
	{
		float sdf = ShapeDiameter_calculate_simple(mesh, indices[i]);
		results.push_back(sdf);
	}
	float maximum = -INFINITY; 
	float minimum = INFINITY; 
	for (size_t i = 0; i < results.size(); i++)
	{
		if (results[i] > maximum)
		{
			maximum = results[i];
		}
		if (results[i] < minimum)
		{
			minimum = results[i];
		}
	}
	return maximum-minimum; 
}
float computeSDF_index(TrilateralMesh* m,int index , int numRays = 10) 
{
	float sdf_value = 0 ;

	glm::vec3& point = m->vertices[index];
	glm::vec3& normal = m->normals[index];
	float sdf = 0.0f;

	// Shoot multiple rays in random directions around the normal
	std::vector<TrilateralRay> rays;
	std::vector<float> ray_distances;
	for (int j = 0; j < numRays; ++j) 
	{
		// Slight random perturbation to the normal
		glm::vec3 rayDir = normal; /* + glm::normalize(glm::vec3((std::rand() % 100 - 50) / 100.0f,
			(std::rand() % 100 - 50) / 100.0f,
			(std::rand() % 100 - 50) / 100.0f));*/
		int rot_x = (std::rand() % 30) - 15;
		int rot_y = (std::rand() % 30) - 15;
		glm::vec3 right, up;
		TrilateralRay ray;
		ray.origin = m->vertices[index];
		ray.direction = -rayDir;
		calculateRightUpVectors(ray.direction, ray.origin, right, up);
		glm::quat quaternion_right = glm::angleAxis((float)rot_x, up);
		glm::quat quaternion_up = glm::angleAxis((float)rot_y, right);
		//glm::quat quaternion_right{} = glm::angleAxis(rot_y, right);
		ray.direction = quaternion_right * ray.direction;
		ray.direction = quaternion_up * ray.direction;
		ray.direction = glm::normalize(ray.direction);
		//rotate direction
		// rotate 
		// Find intersection distance
		float smallest_dist = INFINITY;
		for (size_t k = 0; k < m->triangles.size(); k += 3)
		{
			glm::vec3 p1 = m->vertices[m->triangles[k]];
			glm::vec3 p2 = m->vertices[m->triangles[k + 1]];
			glm::vec3 p3 = m->vertices[m->triangles[k + 2]];
			glm::vec3 hit_point;
			bool is_hit = ray_triangle_intersection(ray, p1, p2, p3, hit_point);
			if (is_hit)
			{
				float distance = glm::distance(m->vertices[index], hit_point);
				if (distance < smallest_dist && distance > 1e-12)
				{
					smallest_dist = distance;
				}
			}

		}
		if (smallest_dist != INFINITY && !std::isnan(smallest_dist))
		{
			ray_distances.push_back(smallest_dist);
			rays.push_back(ray);
		}
		//sdf += smallest_dist;

	}
	// 1 -calculate mean
	float mean = 0;
	for (size_t j = 0; j < ray_distances.size(); j++)
	{
		mean += ray_distances[j];
	}
	mean /= ray_distances.size();
	float variance = 0.0f;
	for (size_t j = 0; j < ray_distances.size(); j++)
	{
		variance += std::pow(ray_distances[j] - mean, 2);
	}
	variance /= ray_distances.size();

	float std_deviation = std::sqrt(variance);

	std::vector<float> standart_deviation_distances;
	for (size_t j = 0; j < ray_distances.size(); j++)
	{
		if (mean - std_deviation <= ray_distances[j] && mean + std_deviation >= ray_distances[j])
		{
			standart_deviation_distances.push_back(ray_distances[j]);
		}
	}
	// Average the distances
	for (size_t j = 0; j < standart_deviation_distances.size(); j++)
	{
		sdf += standart_deviation_distances[j];
	}
	if (standart_deviation_distances.size() > 0)
	{
		sdf_value  = sdf / standart_deviation_distances.size();
	}

	return sdf_value;
}

std::vector<float> computeSDF(TrilateralMesh* m,  int numRays = 10) {
	std::vector<float> sdfValues(m->vertices.size());

	for (size_t i = 0; i < m->vertices.size(); ++i) {
		glm::vec3& point = m->vertices[i];
		glm::vec3& normal = m->normals[i];
		float sdf = 0.0f;

		// Shoot multiple rays in random directions around the normal
		std::vector<TrilateralRay> rays;
		std::vector<float> ray_distances;
		for (int j = 0; j < numRays; ++j) {
			// Slight random perturbation to the normal
			glm::vec3 rayDir = normal; /* + glm::normalize(glm::vec3((std::rand() % 100 - 50) / 100.0f,
				(std::rand() % 100 - 50) / 100.0f,
				(std::rand() % 100 - 50) / 100.0f));*/
			int rot_x = (std::rand() % 30) - 15; 
			int rot_y = (std::rand() % 30) - 15;
			glm::vec3 right, up; 
			TrilateralRay ray;
			ray.origin = m->vertices[i];
			ray.direction = -rayDir;
			calculateRightUpVectors(ray.direction, ray.origin, right, up);
			glm::quat quaternion_right = glm::angleAxis((float)rot_x, up);
			glm::quat quaternion_up = glm::angleAxis((float)rot_y, right);
			//glm::quat quaternion_right{} = glm::angleAxis(rot_y, right);
			ray.direction = quaternion_right * ray.direction;
			ray.direction = quaternion_up * ray.direction;
			ray.direction = glm::normalize(ray.direction);
			//rotate direction
			// rotate 
			// Find intersection distance
			float smallest_dist = INFINITY; 
			for (size_t k = 0; k < m->triangles.size(); k+=3)
			{
				glm::vec3 p1 = m->vertices[m->triangles[k]];
				glm::vec3 p2 = m->vertices[m->triangles[k+1]];
				glm::vec3 p3 = m->vertices[m->triangles[k+2]];
				glm::vec3 hit_point;
				bool is_hit= ray_triangle_intersection(ray, p1,p2,p3, hit_point);
				if (is_hit)
				{
					float distance = glm::distance(m->vertices[i], hit_point);
					if (distance < smallest_dist && distance > 1e-12)
					{
						smallest_dist = distance;
					}
				}

			}
			if (smallest_dist != INFINITY &&  !std::isnan(smallest_dist))
			{
				ray_distances.push_back(smallest_dist);
				rays.push_back(ray);
			}
			//sdf += smallest_dist;
			 
		}
		// 1 -calculate mean
		float mean = 0;
		for (size_t j = 0; j < ray_distances.size(); j++)
		{
			mean += ray_distances[j];
		}
		mean /= ray_distances.size();
		float variance = 0.0f;
		for (size_t j = 0; j < ray_distances.size(); j++)
		{
			variance += std::pow(ray_distances[j] - mean, 2);
		}
		variance /= ray_distances.size();

		float std_deviation = std::sqrt(variance);

		std::vector<float> standart_deviation_distances;
		for (size_t j = 0; j < ray_distances.size(); j++)
		{
			if (mean - std_deviation <= ray_distances[j] && mean + std_deviation >= ray_distances[j])
			{
				standart_deviation_distances.push_back(ray_distances[j]);
			}
		}
		// Average the distances
		for (size_t j = 0; j  < standart_deviation_distances.size(); j++)
		{
			sdf += standart_deviation_distances[j];
		}
		if (standart_deviation_distances.size() > 0)
		{
			sdfValues[i] = sdf / standart_deviation_distances.size();
		}
	}
	m->sdf = sdfValues;
	return sdfValues;
}

void ShapeDiameter_color(TrilateralMesh* m)
{
	std::vector<float> sdf_values;
	float max = -INFINITY;
	float min = INFINITY;
	std::vector<unsigned int> indices;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		indices.push_back(i);
	}
	//ShapeDiameter_calculate(m, indices, sdf_values);
	sdf_values = computeSDF(m, 100);
	/*for (size_t i = 0; i < sdf_values.size(); i++)
	{
		if (sdf_values[i] == INFINITY)
		{
			std::vector<float> smooth_values; 
			//smooth it
			for (size_t j = 0; j < m->adjacenies[i].size(); j++)
			{
				if (sdf_values[j] != INFINITY)
				{
					smooth_values.push_back(sdf_values[j]);
				}
			}
			float sum = 0;
			for (size_t j = 0; j < smooth_values.size(); j++)
			{
				sum += smooth_values[j];
			}
			sum = sum / smooth_values.size();
			sdf_values[i] = sum; 
		}
	}*/

	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		float val = sdf_values[i];
		if (val > max && val != INFINITY )
		{
			max = val;
		}
		if (val < min)
		{
			min = val;
		}
		sdf_values.push_back(val);
	}
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		glm::vec3 color = CoreType_getColor(sdf_values[i] , min , max );
		m->raylib_mesh.colors[i * 4] = color[0] * 255;
		m->raylib_mesh.colors[i * 4 + 1] = color[1] * 255;
		m->raylib_mesh.colors[i * 4 + 2] = color[2] * 255;
		m->raylib_mesh.colors[i * 4 + 3] = 255;
	}
	m->update_raylib_mesh();
}
