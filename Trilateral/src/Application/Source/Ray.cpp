#include "../Include/Ray.h"
#include "../Include/CoreTypeDefs.h"
#include "glm/glm.hpp"

bool ray_triangle_intersection(Ray& ray, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	float A, B, C, D;
    float t;
    const float FLOAT_COMPARISON = 1e-6;
    Plane plane= generate_plane_from_two_vectors(p2 - p1, p3 - p1);
	get_coefficients_from_plane(plane, A, B, C, D);
    float denom = glm::dot(plane.normal, ray.direction);
    if (denom > 1e-6) {
        glm::vec3 p0l0 = plane.point - ray.origin;
        t = glm::dot(p0l0, plane.normal) / denom;
    }
    if (t < 0)
    {
        return false;
    }
    //now check if inside a triangle 

    //barycentric coord

    glm::vec3 Q = ray.origin + ray.direction * t;

    //generate 3 triangles
    float triangle1_area = glm::length(glm::cross(p1 - Q, p2 - Q));  //ABQ
    float triangle2_area = glm::length(glm::cross(p1 - Q, p3 - Q));  //ABQ
    float triangle3_area = glm::length(glm::cross(p2 - Q, p3 - Q));  //ABQ
    
    float main_triangle = glm::length(glm::cross(p1 - Q, p2 - Q));

    if (abs(main_triangle - (triangle1_area + triangle2_area + triangle3_area) ) <= FLOAT_COMPARISON)
    {
        return true;
    }
    return false; 
}