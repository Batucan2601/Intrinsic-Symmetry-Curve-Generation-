#include "../Include/Ray.h"
#include "../Include/CoreTypeDefs.h"
#include "glm/glm.hpp"

bool ray_triangle_intersection(TrilateralRay& ray, glm::vec3 vertex0, glm::vec3 vertex1, glm::vec3 vertex2, glm::vec3& hitpoint)
{
    const float EPSILON = 1e-6f; // Adjusted for better numerical precision
    glm::vec3 edge1 = vertex1 - vertex0;
    glm::vec3 edge2 = vertex2 - vertex0;

    glm::vec3 h = glm::cross(ray.direction, edge2);
    float a = glm::dot(edge1, h);

    // Check if the ray is parallel to the triangle
    if (std::abs(a) < EPSILON) {
        return false;
    }

    float f = 1.0f / a;
    glm::vec3 s = ray.origin - vertex0;
    float u = f * glm::dot(s, h);

    // Check if the intersection is outside the triangle
    if (u < 0.0f || u > 1.0f) {
        return false;
    }

    glm::vec3 q = glm::cross(s, edge1);
    float v = f * glm::dot(ray.direction, q);

    // Check if the intersection is outside the triangle
    if (v < 0.0f || u + v > 1.0f) {
        return false;
    }

    // Calculate the distance to the intersection point
    float t = f * glm::dot(edge2, q);

    if (t > EPSILON) { // Intersection exists in the positive ray direction
        hitpoint = ray.origin + t * ray.direction;
        return true;
    }

    return false; // Intersection behind the ray origin or too close
}