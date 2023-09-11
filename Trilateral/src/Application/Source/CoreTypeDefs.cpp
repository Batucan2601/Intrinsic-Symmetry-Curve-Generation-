#include "../Include/CoreTypeDefs.h"

Mesh generate_mesh_from_plane( Plane* plane, glm::vec3 * m)
{
	// a(x - x0) + b(y - y0) + c(z - z0) = 0
	// only increase in y and z 
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;
	glm::vec3 p4;

	//increase x by two  increase y by two
	p1.x = m->x + 20.0f;
	p1.y = m->y + 20.0f;
	p1.z = (plane->normal.x * (p1.x - m->x) + plane->normal.y * (p1.y - m->y)) / plane->normal.z + m->z;

	// increase x by two  decrease y by two 
	p2.x = m->x + 20.0f;
	p2.y = m->y - 20.0f;
	p2.z = (plane->normal.x * (p2.x - m->x) + plane->normal.y * (p2.y - m->y)) / plane->normal.z + m->z;

	// decrease x by two  decrease y by two 
	p3.x = m->x - 20.0f ;
	p3.y = m->y - 20.0f;
	p3.z = (plane->normal.x * (p3.x - m->x) + plane->normal.y * (p3.y - m->y)) / plane->normal.z + m->z;

	// decrease x by two increase y by two 
	p4.x = m->x - 20.0f;
	p4.y = m->y + 20.0f;
	p4.z = (plane->normal.x * (p4.x - m->x) + plane->normal.y * (p4.y - m->y)) / plane->normal.z + m->z;

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