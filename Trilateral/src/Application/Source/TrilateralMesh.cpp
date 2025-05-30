#include "../Include/TrilateralMesh.h"
#include "happly.h"
#include "glm/gtc/matrix_transform.hpp"
#include <glm/gtc/type_ptr.hpp>
#include "../Include/CoreTypeDefs.h"
#include "../Include/ShapeDiameter.h"
#include "../Include/Geodesic.h"

static void calculate_areas(TrilateralMesh* m);
static float compute_triangle_area(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
static void calculate_mesh_area(TrilateralMesh* m);

TrilateralMesh::TrilateralMesh()
{

}
TrilateralMesh::TrilateralMesh(char* filename)
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
	this->file_name = filename;
	//get without slashes
	int slash_index = -1;
	for (size_t i = 0; i < this->file_name.size(); i++)
	{
		if (this->file_name[i] == '/' || this->file_name[i] == '\\')
		{
			slash_index = i;
		}
	}
	this->file_name = this->file_name.substr( slash_index+1, this->file_name.size());
	generate_normals();
	calculate_areas(this);
	calculate_mesh_area(this);
	this->generate_raylib_mesh();
	this->calculate_PCA();
}
void TrilateralMesh::read_ply_format(char* filename)
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
static int countFloatsInLine(const std::string& line) {
	std::stringstream ss(line);
	float value;
	int count = 0;

	// Try to extract floats from the line
	while (ss >> value) {
		count++;
		// Clear the bad input state if there's any non-float
		ss.clear();
		// Ignore any non-numeric characters after the float
		ss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
	}

	return count;
}
void TrilateralMesh::read_off_format(char* filename)
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
	float r, g, b;
	r = 0; g = 0; b = 0;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nTris, &n);
	// fill the needed vectors
	std::vector<std::pair<int, float>> temp_vector;
	for (size_t i = 0; i < nVerts; i++)
	{
		adjacenies.push_back(temp_vector);
	}
	bool isColorExist = false;
	int originalPos = ftell(fPtr);
	char buffer[100];
	if (fgets(buffer, sizeof(buffer), fPtr)) {
		printf("Temporarily read line: %s", buffer);
	}
	fseek(fPtr, originalPos-2, SEEK_SET);
	int float_count = countFloatsInLine(buffer);
	if (float_count == 6)
	{
		isColorExist = true;
	}
	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		if (isColorExist)
		{
			fscanf(fPtr, "%f %f %f", &r, &g, &b);
		}
		glm::vec3 point = glm::vec3(x, y, z);
		vertices.push_back(point);
		this->colors.push_back(glm::vec3(r,g,b));
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
		bool is_x = false;
		bool is_y = false;
		bool is_z = false;
		for (size_t i = 0; i < adjacenies[x].size(); i++)
		{
			if (adjacenies[x][i].first == (int)y)
			{
				is_y = true; 
			}
			if (adjacenies[x][i].first == (int)z)
			{
				is_z = true;
			}
		}
		if (!is_y)
		{
			adjacenies[x].push_back(std::make_pair(y, e1.distance));
		}
		if (!is_z)
		{
			adjacenies[x].push_back(std::make_pair(z, e2.distance));
		}
		is_x = false;
		is_y = false;
		is_z = false;
		for (size_t i = 0; i < adjacenies[y].size(); i++)
		{
			if (adjacenies[y][i].first == (int)x)
			{
				is_x = true;
			}
			if (adjacenies[y][i].first == (int)z)
			{
				is_z = true;
			}
		}
		//for second 
		if (!is_x)
		{
			adjacenies[y].push_back(std::make_pair(x, e1.distance));
		}
		if (!is_z)
		{
			adjacenies[y].push_back(std::make_pair(z, e3.distance));
		}
		is_x = false;
		is_y = false;
		is_z = false;
		for (size_t i = 0; i < adjacenies[z].size(); i++)
		{
			if (adjacenies[z][i].first == (int)x)
			{
				is_x = true;
			}
			if (adjacenies[z][i].first == (int)y)
			{
				is_z = true;
			}
		}
		//for third
		if (!is_x)
		adjacenies[z].push_back(std::make_pair(x, e2.distance));

		if (!is_y)
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
TrilateralMesh::TrilateralMesh(glm::vec3* p1, glm::vec3* p2, glm::vec3* p3, glm::vec3* p4)
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


	this->raylib_mesh.vertices = (float*)&this->vertices[0];
	this->raylib_mesh.animVertices = NULL;
	this->raylib_mesh.normals = NULL;
	this->raylib_mesh.tangents = NULL;
	this->raylib_mesh.texcoords2 = NULL;
	this->raylib_mesh.boneIds = NULL;
	this->raylib_mesh.boneWeights = NULL;

	this->raylib_mesh.vertexCount = this->vertices.size();
	this->raylib_mesh.triangleCount = this->triangles.size() / 3;
	//copy indices
	this->raylib_mesh.indices = (unsigned short*)malloc(this->triangles.size() * sizeof(unsigned short));
	for (size_t i = 0; i < this->triangles.size(); i++)
	{
		this->raylib_mesh.indices[i] = this->triangles[i];
	}
	std::vector<unsigned char > zeroes(this->vertices.size() * 4, 0);
	this->raylib_mesh.colors = (unsigned char*)malloc(this->vertices.size() * 4 * 1);
	this->raylib_mesh.texcoords = (float*)malloc(this->vertices.size() * 2 * 4);
	this->raylib_mesh.vaoId = 0;

}
glm::mat4 TrilateralMesh::move_mesh(glm::vec3 direction)
{
	glm::mat4 temp = glm::translate(glm::mat4(1.0f), direction);
	return temp;
}
glm::mat4 TrilateralMesh::scale_mesh(glm::vec3 scale)
{
	glm::mat4 temp = glm::scale(glm::mat4(1.0f), scale);
	return temp;
}

void TrilateralMesh::generate_raylib_mesh()
{

	this->raylib_mesh.vertices = (float*)malloc(this->vertices.size() * 3 * sizeof(float));
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		this->raylib_mesh.vertices[i * 3] = this->vertices[i].x;
		this->raylib_mesh.vertices[i * 3 + 1] = this->vertices[i].y;
		this->raylib_mesh.vertices[i * 3 + 2] = this->vertices[i].z;
	}
	this->raylib_mesh.animVertices = NULL;
	this->raylib_mesh.boneMatrices = NULL;
	this->raylib_mesh.animNormals = NULL;
	this->raylib_mesh.normals = (float*)malloc(this->vertices.size() * 3 * sizeof(float));
	for (size_t i = 0; i < this->normals.size(); i++)
	{
		glm::vec3 normal = this->normals[i];
		this->raylib_mesh.normals[i * 3] = normal.x;
		this->raylib_mesh.normals[i * 3 + 1] = normal.y;
		this->raylib_mesh.normals[i * 3 + 2] = normal.z;
	}
	this->raylib_mesh.tangents = NULL;
	this->raylib_mesh.texcoords2 = NULL;
	this->raylib_mesh.boneIds = NULL;
	this->raylib_mesh.boneWeights = NULL;

	this->raylib_mesh.vertexCount = this->vertices.size();
	this->raylib_mesh.triangleCount = this->triangles.size() / 3;
	//copy indices
	this->raylib_mesh.indices = (unsigned short*)malloc(this->triangles.size() * sizeof(unsigned short));
	for (size_t i = 0; i < this->triangles.size(); i++)
	{
		this->raylib_mesh.indices[i] = this->triangles[i];
	}
	std::vector<unsigned char > zeroes(this->vertices.size() * 4, 0);
	this->raylib_mesh.colors = (unsigned char*)malloc(this->vertices.size() * 4 * 1);
	this->raylib_mesh.texcoords = (float*)malloc(this->vertices.size() * 2 * 4);
	this->raylib_mesh.vaoId = 0;

	UploadMesh(&this->raylib_mesh, false);
	this->color_all(LIGHTGRAY);
}
void TrilateralMesh::update_raylib_mesh()
{
	//color update
	UpdateMeshBuffer(this->raylib_mesh, 3, this->raylib_mesh.colors, this->vertices.size() * 4 * sizeof(unsigned char), 0); // Buffer index 3 is for color
	
}
void read_symmetry_format(char* filename, TrilateralMesh* m)
{
	std::ifstream symFile(filename);
	double number = 0;
	int index = 0;
	if (!symFile.is_open())
	{
		std::cout << " failed to open the file " << std::endl; 
		return;
	}
	while (symFile >> number)
	{
		std::pair<unsigned int, unsigned int> sym_pair;
		sym_pair.first = index;
		sym_pair.second = (unsigned int)number-1;
		index++;
		m->symmetry_pairs.push_back(sym_pair);

	}
	m->ground_truth_symmetry_pairs = std::vector<unsigned int>(m->symmetry_pairs.size());
	for (size_t i = 0; i < m->symmetry_pairs.size(); i++)
	{
		m->ground_truth_symmetry_pairs[m->symmetry_pairs[i].first] = m->symmetry_pairs[i].second;
	}
	for (size_t i = 0; i < m->symmetry_pairs.size(); i++)
	{
		m->ground_truth_symmetry_pairs[m->symmetry_pairs[i].second] = m->symmetry_pairs[i].first;
	}	
	return;
}

//gives normal for each point 
void TrilateralMesh::generate_normals()
{
	std::vector<glm::vec3> normals(this->vertices.size() , glm::vec3(0.0f,0.0f,0.0f));
	std::vector<float> each_point_count(this->vertices.size(), 0);

	for (size_t i = 0; i < this->triangles.size(); i+=3)
	{
		glm::vec3 edge1 = this->vertices[this->triangles[i+1]] - this->vertices[this->triangles[i]];
		glm::vec3 edge2 = this->vertices[this->triangles[i+2]] - this->vertices[this->triangles[i]];
		glm::vec3 normal = glm::cross(edge1 , edge2);
		each_point_count[this->triangles[i]] += 1;
		each_point_count[this->triangles[i+1]] += 1;
		each_point_count[this->triangles[i+2]] += 1;
		
		normals[this->triangles[i]] = normals[this->triangles[i]] + normal;
		normals[this->triangles[i+1]] = normals[this->triangles[i+1]] + normal;
		normals[this->triangles[i+2]] = normals[this->triangles[i+2]] + normal;
	}
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		normals[i] = normals[i] / each_point_count[i];
		normals[i] = glm::normalize(normals[i]);
	}
	this->normals = normals;

	//now the part where we display the normals
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		//start of the line is vertex itself
		normals_display.push_back(this->vertices[i].x);
		normals_display.push_back(this->vertices[i].y);
		normals_display.push_back(this->vertices[i].z);
		
		normals_display.push_back(0);
		normals_display.push_back(255);
		normals_display.push_back(0);

		// end part is the normal
		glm::vec3 normal = this->vertices[i] + this->normals[i];
		normals_display.push_back(normal.x);
		normals_display.push_back(normal.y);
		normals_display.push_back(normal.z);

		normals_display.push_back(0);
		normals_display.push_back(255);
		normals_display.push_back(0);

	}

}


unsigned int mesh_get_closest_index(TrilateralMesh* m, const glm::vec3& point)
{
	unsigned int closest_distance_index = -1;
	float closest_distance = INFINITY;
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		float dist = glm::distance(m->vertices[i], point );
		if (dist < closest_distance)
		{
			closest_distance = dist; 
			closest_distance_index = i;
		}
	}
	return closest_distance_index;
}


static float compute_triangle_area(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 vec1 = p2 - p1;
	glm::vec3 vec2 = p3 - p1;

	glm::vec3 cross1 = glm::cross(vec1, vec2);
	float length = glm::length(cross1) / 2;
	return length;
}
static void calculate_areas(TrilateralMesh* m)
{
	m->areas = std::vector<float>(m->vertices.size(), 0);
	for (size_t i = 0; i < m->triangles.size(); i+=3)
	{
		glm::vec3 p1 = m->vertices[m->triangles[i]];
		glm::vec3 p2 = m->vertices[m->triangles[i+1]];
		glm::vec3 p3 = m->vertices[m->triangles[i+2]];
		float area = compute_triangle_area(p1, p2, p3);
		m->areas[m->triangles[i]] += area;
		m->areas[m->triangles[i+1]] += area;
		m->areas[m->triangles[i+2]] += area;
	}
}
void TrilateralMesh::calculate_sdf()
{
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		float sdf = ShapeDiameter_calculate_simple(this, i);
		this->sdf.push_back(sdf);
	}
}
static void calculate_mesh_area(TrilateralMesh* m )
{
	m->mesh_area = 0;
	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		glm::vec3 p1 = m->vertices[m->triangles[i]];
		glm::vec3 p2 = m->vertices[m->triangles[i + 1]];
		glm::vec3 p3 = m->vertices[m->triangles[i + 2]];
		float area = compute_triangle_area(p1, p2, p3);
		m->mesh_area += area; 
	}
}

glm::vec3 mesh_generate_weighted_mid_point(TrilateralMesh* m) //supposed to be best way
{
	glm::vec3 midpoint = glm::vec3(0, 0, 0);
	float biggest_triangle = -INFINITY; 

	for (size_t i = 0; i < m->triangles.size(); i+=3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		float tri_area = compute_triangle_area(m->vertices[index1], m->vertices[index2], m->vertices[index3]);

		if (tri_area > biggest_triangle)
		{
			biggest_triangle = tri_area;
		}
	}

	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		float tri_area = compute_triangle_area(m->vertices[index1], m->vertices[index2], m->vertices[index3]);

		midpoint = midpoint + m->vertices[index1] * tri_area / biggest_triangle;
		midpoint = midpoint + m->vertices[index2] * tri_area / biggest_triangle;
		midpoint = midpoint + m->vertices[index3] * tri_area / biggest_triangle;
	}

	midpoint /= m->triangles.size();
	return midpoint;
}

// computes all of the areas of triangle for each point
// normalizes for each point
std::vector<float> mesh_point_surfel_normalized(TrilateralMesh* m)
{
	std::vector<float> triangles_normalized(m->vertices.size() , 0);
	glm::vec3 midpoint = glm::vec3(0, 0, 0);
	float biggest_triangle = -INFINITY;

	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		float tri_area = compute_triangle_area(m->vertices[index1], m->vertices[index2], m->vertices[index3]);

		if (tri_area > biggest_triangle)
		{
			biggest_triangle = tri_area;
		}
	}

	for (size_t i = 0; i < m->triangles.size(); i += 3)
	{
		int index1 = m->triangles[i];
		int index2 = m->triangles[i + 1];
		int index3 = m->triangles[i + 2];
		float tri_area = compute_triangle_area(m->vertices[index1], m->vertices[index2], m->vertices[index3]);

		triangles_normalized[index1] += tri_area / biggest_triangle;
		triangles_normalized[index2] += tri_area / biggest_triangle;
		triangles_normalized[index3] += tri_area / biggest_triangle;
	}

	return triangles_normalized;
}


void TrilateralMesh::color_all(Color color)
{
	std::vector<unsigned int> points; 
	for (size_t i = 0; i < this->vertices.size(); i++)
	{
		points.push_back(i);
	}
	this->color_points(points, color);
}
void TrilateralMesh::color_points(std::vector<unsigned int>& points , Color color)
{
	for (size_t i = 0; i < points.size(); i++)
	{
		this->raylib_mesh.colors[points[i] * 4] = color.r;
		this->raylib_mesh.colors[points[i] * 4 + 1] = color.g;
		this->raylib_mesh.colors[points[i] * 4 + 2] = color.b;
		this->raylib_mesh.colors[points[i] * 4 + 3] = color.a;
	}
	this->update_raylib_mesh();
}

void TrilateralMesh::calculate_PCA()
{
	// generate PCA weights are same and 1 for now 
	float s = this->vertices.size();
	glm::vec3 m(0.0f, 0.0f, 0.0f);
	int N = this->vertices.size();
	m = plane_point;

	Eigen::MatrixXd Co(3, 3);

	Co(0, 0) = 0;
	Co(0, 1) = 0;
	Co(0, 2) = 0;
	Co(1, 0) = 0;
	Co(1, 1) = 0;
	Co(1, 2) = 0;
	Co(2, 0) = 0;
	Co(2, 1) = 0;
	Co(2, 2) = 0;
	for (size_t i = 0; i < N; i++)
	{
		glm::vec3 pi_m;
		pi_m = this->vertices[i] - m;
		Eigen::VectorXd pi(3);
		pi(0) = pi_m.x;
		pi(1) = pi_m.y;
		pi(2) = pi_m.z;

		Eigen::MatrixXd  Co_i = pi * pi.transpose();
		Co = Co + Co_i;
	}
	Co = Co / s;

	//// get the eigenvectors 
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Co);
	Eigen::MatrixXd eigen_vecs = es.eigenvectors().real();
	Eigen::VectorXd eigen_values = es.eigenvalues().real();

	double biggest_value = -INFINITY;
	int biggest_index = -1;
	//get the best eigen value
	for (size_t i = 0; i < eigen_values.rows(); i++)
	{
		if (biggest_value < (float)eigen_values(i))
		{
			biggest_value = (float)eigen_values(i);
			biggest_index = i;
		}
	}
	

	// generate the 3 planes
	//std::vector<float> eigenvalues = { eigen_values.col(0).row(0) , eigen_values.col(0).row(1) , eigen_values.col(0).row(2) };
	this->PCA = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	this->PCA = glm::normalize(PCA);
	/*planes[0].point = m;
	planes[1].normal = glm::vec3(eigen_vecs.col(1).real()(0), eigen_vecs.col(1).real()(1), eigen_vecs.col(1).real()(2));
	planes[1].point = m;
	planes[2].normal = glm::vec3(eigen_vecs.col(2).real()(0), eigen_vecs.col(2).real()(1), eigen_vecs.col(2).real()(2));
	planes[2].point = m;*/
}

void TrilateralMesh::calculate_SDF(int num_rays , float angle )
{
	this->sdf = computeSDF(this, num_rays, angle);
}

void TrilateralMesh::color_midpoints(Color color )
{
	std::vector<unsigned int>midpoints; 
	for (size_t i = 0; i < this->calculated_symmetry_pairs.size(); i++)
	{
		int index1 = this->calculated_symmetry_pairs[i].first;
		int index2 = this->calculated_symmetry_pairs[i].second;
		int mid  = Geodesic_get_midpoint_from_path(this, index1, index2);
		midpoints.push_back(mid);
	}

	this->color_points(midpoints, color);
}