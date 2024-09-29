#include "../Include/Ray.h"
#include "../Include/CoreTypeDefs.h"
#include "glm/glm.hpp"

bool ray_triangle_intersection(TrilateralRay& ray, glm::vec3 vertex0, glm::vec3 vertex1, glm::vec3 vertex2, glm::vec3& hitpoint)
{
    glm::vec3 edge1, edge2, h, s, q;
    float a, f, u, v;

    edge1.x = vertex1.x - vertex0.x;
    edge1.y = vertex1.y - vertex0.y;
    edge1.z = vertex1.z - vertex0.z;
    edge2.x = vertex2.x - vertex0.x;
    edge2.y = vertex2.y - vertex0.y;
    edge2.z = vertex2.z - vertex0.z;

    h = glm::cross(ray.direction, edge2);
    a = glm::dot(edge1, h);

    if (a > -0.00001 && a < 0.00001)
        return false;

    f = 1.0 / a;
    s.x = ray.origin.x - vertex0.x;
    s.y = ray.origin.y - vertex0.y;
    s.z = ray.origin.z - vertex0.z;
    u = f * glm::dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    q = glm::cross(s, edge1);
    v = f * glm::dot(ray.direction, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    float t = f * glm::dot(edge2, q);
    if (t > 0.00001) { // ray intersection
        hitpoint.x = ray.origin.x + ray.direction.x * t;
        hitpoint.y = ray.origin.y + ray.direction.y * t;
        hitpoint.z = ray.origin.z + ray.direction.z * t;
        return true;
    }

    return false;
}