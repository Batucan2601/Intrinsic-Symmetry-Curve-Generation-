#include "../Include/Ray.h"
#include "../Include/CoreTypeDefs.h"
#include "glm/glm.hpp"

bool ray_triangle_intersection(TrilateralRay& ray, glm::vec3 vertex0, glm::vec3 vertex1, glm::vec3 vertex2, glm::vec3& hitpoint)
{
    glm::vec3 edge1, edge2, h, s, q;
    float a, f, u, v;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = glm::cross(ray.direction, edge2);
    a = glm::dot(edge1, h);

    if (a > -0.00001 && a < 0.00001)
        return false;

    f = 1.0 / a;
    s = ray.origin - vertex0;
    u = f * glm::dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    q = glm::cross(s, edge1);
    v = f * glm::dot(ray.direction, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    float t = f * glm::dot(edge2, q);
    if (t > 0) { // ray intersection
        hitpoint.x = ray.origin.x + ray.direction.x * t;
        hitpoint.y = ray.origin.y + ray.direction.y * t;
        hitpoint.z = ray.origin.z + ray.direction.z * t;
        return true;
    }

    return false;
}