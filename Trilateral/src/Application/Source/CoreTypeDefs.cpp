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
	p1.x = m->x + 2.0f;
	p1.y = m->y + 2.0f;
	p1.z = (plane->normal.x * (p1.x - m->x) + plane->normal.y * (p1.y - m->y)) / plane->normal.z + m->z;

	// increase x by two  decrease y by two 
	p2.x = m->x + 2.0f;
	p2.y = m->y - 2.0f;
	p2.z = (plane->normal.x * (p2.x - m->x) + plane->normal.y * (p2.y - m->y)) / plane->normal.z + m->z;

	// decrease x by two  decrease y by two 
	p3.x = m->x - 2.0f ;
	p3.y = m->y - 2.0f;
	p3.z = (plane->normal.x * (p3.x - m->x) + plane->normal.y * (p3.y - m->y)) / plane->normal.z + m->z;

	// decrease x by two increase y by two 
	p4.x = m->x - 2.0f;
	p4.y = m->y + 2.0f;
	p4.z = (plane->normal.x * (p4.x - m->x) + plane->normal.y * (p4.y - m->y)) / plane->normal.z + m->z;

	Mesh plane_mesh(&p1,&p2,&p3,&p4); 
	return plane_mesh;
}