#include "../Include/ShapeDiameter.h"
#include <src/Application/Include/Ray.h>
#include "glm/glm.hpp"
#define _USE_MATH_DEFINES //for pi 
#include <math.h>
#include <random>


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

void ShapeDiameter_calculate(Mesh* mesh,  std::vector<unsigned int> indices , std::vector<float>& shape_diameter)
{

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> random_r (1, CONE_RADIUS * CONE_RADIUS_MULTIPLIER); // distribution in range [1, 6]

	// save each of the distances
	// for each endpoint
	for (size_t i = 0; i < indices.size(); i++)
	{
		//
		std::vector<float> ray_distances; 
		//generate ray
		Ray ray; 
		ray.origin = mesh->vertices[indices[i]];
		ray.direction = mesh->normals[indices[i]] * -1.0f; //inverse of normal
		glm::vec3 cone_center = ray.origin + ray.direction * (float)CONE_HEIGHT;
		for (size_t j = 0; j < NUMBER_OF_RAYS; j++)
		{
			glm::vec3 right;
			glm::vec3 left;
			float r = CONE_RADIUS * (float)std::sqrt(random_r(rng)) / (float)std::sqrt(CONE_RADIUS * CONE_RADIUS_MULTIPLIER);
			float theta = (float)std::sqrt(random_r(rng)) / (float)std::sqrt(CONE_RADIUS * CONE_RADIUS_MULTIPLIER) * 2 * M_PI;
			//generate right and up vector
			calculateRightUpVectors(  ray.direction , cone_center ,right , left  );

			Ray random_ray;
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
					if (dist < dist)
					{
						closest_triangle_index = k;
						closest_dist  = dist;

					}
				}

			}
			//compute the closest again
			glm::vec3 p1 = mesh->vertices[mesh->triangles[closest_triangle_index]];
			glm::vec3 p2 = mesh->vertices[mesh->triangles[closest_triangle_index + 1]];
			glm::vec3 p3 = mesh->vertices[mesh->triangles[closest_triangle_index + 2]];
			ray_triangle_intersection(random_ray, p1, p2, p3, hit_point);
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

