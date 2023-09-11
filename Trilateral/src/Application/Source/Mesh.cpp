#include "../Include/Mesh.h"
#include "happly.h"
#include "glm/gtc/matrix_transform.hpp"


Mesh::Mesh()
{

}
Mesh::Mesh(char* filename)
{
	model_mat = glm::mat4(1.0f);
	std::string string_file_name = (std::string)filename;
	if (string_file_name.find(".off") != std::string::npos)
	{
		off_format = true;
		read_off_format(filename);

	}
	else if (string_file_name.find(".ply") != std::string::npos)
	{
		ply_format = true;
		read_ply_format(filename);

	}
}
void Mesh::read_ply_format(char* filename)
{
	happly::PLYData plyIn(filename);
	std::vector<std::array<double, 3U > > vertex_pos = plyIn.getVertexPositions();
	std::vector<std::vector<size_t>> face_indices = plyIn.getFaceIndices();

	for (size_t i = 0; i < vertex_pos.size(); i++)
	{
		this->vertices.push_back(glm::vec3(vertex_pos[i][0], vertex_pos[i][1], vertex_pos[i][2]));
		this->colors.push_back(glm::vec3(0.0f, 0.0f, 0.0f));

	}

	// fill the needed vectors
	std::vector<std::pair<int, float>> temp_vector;
	for (size_t i = 0; i < vertex_pos.size(); i++)
	{
		adjacenies.push_back(temp_vector);
	}
	for (size_t i = 0; i < face_indices.size(); i++)
	{
		Edge e1;
		Edge e2;
		Edge e3;

		// edge 1
		e1.p1 = this->vertices[(float)face_indices[i][0]];
		e1.p2 = this->vertices[(float)face_indices[i][1]];
		e1.distance = glm::distance(e1.p1, e1.p2);
		// edge 1
		e2.p1 = this->vertices[(float)face_indices[i][0]];
		e2.p2 = this->vertices[(float)face_indices[i][2]];
		e2.distance = glm::distance(e2.p1, e2.p2);

		// edge 3
		// edge 1
		e3.p1 = this->vertices[(float)face_indices[i][1]];
		e3.p2 = this->vertices[(float)face_indices[i][2]];
		e3.distance = glm::distance(e3.p1, e3.p2);

		this->edges.push_back(e1);
		this->edges.push_back(e2);
		this->edges.push_back(e3);

		//edge1 
		this->vertex_indices.push_back((int)face_indices[i][0]);
		this->vertex_indices.push_back((int)face_indices[i][1]);
		//edge2 
		this->vertex_indices.push_back((int)face_indices[i][0]);
		this->vertex_indices.push_back((int)face_indices[i][2]);
		//edge3 
		this->vertex_indices.push_back(face_indices[i][1]);
		this->vertex_indices.push_back(face_indices[i][2]);


		//for first 
		adjacenies[face_indices[i][0]].push_back(std::make_pair(face_indices[i][1], e1.distance));
		adjacenies[face_indices[i][0]].push_back(std::make_pair(face_indices[i][2], e2.distance));
		//for second 
		adjacenies[face_indices[i][1]].push_back(std::make_pair(face_indices[i][0], e1.distance));
		adjacenies[face_indices[i][1]].push_back(std::make_pair(face_indices[i][2], e3.distance));
		//for third
		adjacenies[face_indices[i][2]].push_back(std::make_pair(face_indices[i][1], e3.distance));
		adjacenies[face_indices[i][2]].push_back(std::make_pair(face_indices[i][0], e2.distance));

		//triangles
		triangles.push_back(face_indices[i][0]);
		triangles.push_back(face_indices[i][1]);
		triangles.push_back(face_indices[i][2]);
	}

	//neighbourhood search 
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		std::vector<unsigned int> neighbour_vec;
		this->neighbours.push_back(neighbour_vec);
	}
	for (size_t i = 0; i < this->triangles.size(); i+= 3 )
	{
		bool p1_exists = false;
		bool p2_exists = false;
		bool p3_exists = false;
		unsigned int p1 = this->triangles[i];
		unsigned int p2 = this->triangles[i+1];
		unsigned int p3 = this->triangles[i+2];
		
		//for p1 if p2 or p3 does not exists add them
		for (size_t j = 0; j < this->neighbours[p1].size(); j++)
		{
			if (this->neighbours[p1][j] == p2)
			{
				p2_exists = true;
			}
			if (this->neighbours[p1][j] == p3)
			{
				p3_exists = true;
			}
		}
		if (!p2_exists)
		{
			this->neighbours[p1].push_back(p2);
		}
		if (!p3_exists)
		{
			this->neighbours[p1].push_back(p3);
		}

		p2_exists = false;
		p3_exists = false;

		//for p1 if p2 or p3 does not exists add them
		for (size_t j = 0; j < this->neighbours[p2].size(); j++)
		{
			if (this->neighbours[p2][j] == p1)
			{
				p1_exists = true;
			}
			if (this->neighbours[p2][j] == p3)
			{
				p3_exists = true;
			}
		}
		if (!p1_exists)
		{
			this->neighbours[p2].push_back(p1);
		}
		if (!p3_exists)
		{
			this->neighbours[p2].push_back(p3);
		}

		p1_exists = false;
		p3_exists = false;

		//for p1 if p2 or p3 does not exists add them
		for (size_t j = 0; j < this->neighbours[p3].size(); j++)
		{
			if (this->neighbours[p3][j] == p1)
			{
				p1_exists = true;
			}
			if (this->neighbours[p3][j] == p2)
			{
				p2_exists = true;
			}
		}
		if (!p1_exists)
		{
			this->neighbours[p3].push_back(p1);
		}
		if (!p2_exists)
		{
			this->neighbours[p3].push_back(p2);
		}

		p1_exists = false;
		p2_exists = false;
	}
}
void Mesh::read_off_format(char* filename)
{
	// standart readeing 
	errno = 0;

	FILE* fPtr = fopen(filename, "r");
	if (errno != 0)
	{
		perror("Error occurred while opening file.\n");
		exit(1);
	}
	char str[334];
	fscanf(fPtr, "%s", str);
	
	int nVerts, nTris, n, i = 0;
	float x, y, z;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nTris, &n);
	// fill the needed vectors
	std::vector<std::pair<int, float>> temp_vector;
	for (size_t i = 0; i < nVerts; i++)
	{

		adjacenies.push_back(temp_vector);
	}

	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		glm::vec3 point = glm::vec3(x, y, z);
		vertices.push_back(point);
		this->colors.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
	}
	while (fscanf(fPtr, "%d", &i) != EOF)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		Edge e1;
		Edge e2;
		Edge e3;
		// edge 1
		e1.p1 = this->vertices[x];
		e1.p2 = this->vertices[y];
		e1.distance = glm::distance(e1.p1, e1.p2);
		// edge 1
		e2.p1 = this->vertices[x];
		e2.p2 = this->vertices[z];
		e2.distance = glm::distance(e2.p1, e2.p2);

		// edge 3
		// edge 1
		e3.p1 = this->vertices[y];
		e3.p2 = this->vertices[z];
		e3.distance = glm::distance(e3.p1, e3.p2);

		edges.push_back(e1);
		edges.push_back(e2);
		edges.push_back(e3);


		//indices of vertices 
		//edge1 
		vertex_indices.push_back(x);
		vertex_indices.push_back(y);
		//edge2 
		vertex_indices.push_back(x);
		vertex_indices.push_back(z);
		//edge3 
		vertex_indices.push_back(y);
		vertex_indices.push_back(z);

		//adjacencies

		//for first 
		adjacenies[x].push_back(std::make_pair(y, e1.distance));
		adjacenies[x].push_back(std::make_pair(z, e2.distance));
		//for second 
		adjacenies[y].push_back(std::make_pair(x, e1.distance));
		adjacenies[y].push_back(std::make_pair(z, e3.distance));
		//for third
		adjacenies[z].push_back(std::make_pair(x, e2.distance));
		adjacenies[z].push_back(std::make_pair(y, e3.distance));

		//triangles
		triangles.push_back(x);
		triangles.push_back(y);
		triangles.push_back(z);


	}
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		std::vector<unsigned int> neighbour_vec;
		this->neighbours.push_back(neighbour_vec);
	}
	for (size_t i = 0; i < this->triangles.size(); i += 3)
	{
		bool p1_exists = false;
		bool p2_exists = false;
		bool p3_exists = false;
		unsigned int p1 = this->triangles[i];
		unsigned int p2 = this->triangles[i + 1];
		unsigned int p3 = this->triangles[i + 2];

		//for p1 if p2 or p3 does not exists add them
		for (size_t j = 0; j < this->neighbours[p1].size(); j++)
		{
			if (this->neighbours[p1][j] == p2)
			{
				p2_exists = true;
			}
			if (this->neighbours[p1][j] == p3)
			{
				p3_exists = true;
			}
		}
		if (!p2_exists)
		{
			this->neighbours[p1].push_back(p2);
		}
		if (!p3_exists)
		{
			this->neighbours[p1].push_back(p3);
		}

		p2_exists = false;
		p3_exists = false;

		//for p1 if p2 or p3 does not exists add them
		for (size_t j = 0; j < this->neighbours[p2].size(); j++)
		{
			if (this->neighbours[p2][j] == p1)
			{
				p1_exists = true;
			}
			if (this->neighbours[p2][j] == p3)
			{
				p3_exists = true;
			}
		}
		if (!p1_exists)
		{
			this->neighbours[p2].push_back(p1);
		}
		if (!p3_exists)
		{
			this->neighbours[p2].push_back(p3);
		}

		p1_exists = false;
		p3_exists = false;

		//for p1 if p2 or p3 does not exists add them
		for (size_t j = 0; j < this->neighbours[p3].size(); j++)
		{
			if (this->neighbours[p3][j] == p1)
			{
				p1_exists = true;
			}
			if (this->neighbours[p3][j] == p2)
			{
				p2_exists = true;
			}
		}
		if (!p1_exists)
		{
			this->neighbours[p3].push_back(p1);
		}
		if (!p2_exists)
		{
			this->neighbours[p3].push_back(p2);
		}

		p1_exists = false;
		p2_exists = false;
	}
	fclose(fPtr);
}
//plane constructor 
Mesh::Mesh(glm::vec3* p1, glm::vec3* p2, glm::vec3* p3, glm::vec3* p4)
{
	this->model_mat = glm::mat4(1.0f);

	this->vertices.push_back(*p1 );
	this->vertices.push_back(*p2);
	this->vertices.push_back(*p3);
	this->vertices.push_back(*p4);

	this->triangles.push_back(0);
	this->triangles.push_back(1);
	this->triangles.push_back(2);

	this->triangles.push_back(2);
	this->triangles.push_back(3);
	this->triangles.push_back(0);

	this->colors.push_back(glm::vec3(1.0f, 0.5f, 0.5f));
	this->colors.push_back(glm::vec3(1.0f, 0.5f, 0.5f));
	this->colors.push_back(glm::vec3(1.0f, 0.5f, 0.5f));
	this->colors.push_back(glm::vec3(1.0f, 0.5f, 0.5f));

}
glm::mat4 Mesh::move_mesh(glm::vec3 direction)
{
	glm::mat4 temp = glm::translate(glm::mat4(1.0f), direction);
	return temp;
}
glm::mat4 Mesh::scale_mesh(glm::vec3 scale)
{
	glm::mat4 temp = glm::scale(glm::mat4(1.0f), scale);
	return temp;
}


void read_symmetry_format(char* filename, Mesh* m)
{
	std::ifstream symFile(filename);
	double number = 0;
	int index = 0;
	while (symFile >> number)
	{
		std::pair<unsigned int, unsigned int> sym_pair;
		sym_pair.first = index;
		sym_pair.second = (unsigned int)number-1;
		index++;
		m->symmetry_pairs.push_back(sym_pair);
	}
	/*while (symFile >> number)
	{
		std::pair<unsigned int, unsigned int> sym_pair;
		sym_pair.first = (unsigned int)number; 
		symFile >> number;
		sym_pair.second= (unsigned int)number;
		m->symmetry_pairs.push_back(sym_pair);
	}*/

}