#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include "glm/glm.hpp"
using std::ofstream;
using std::ifstream;
#define INFINITE 10000000

struct Edge
{
	glm::vec3 p1;
	glm::vec3 p2;
	float distance;
};
class Mesh
{
public:
	std::vector<glm::vec3> vertices; //fill with constructor  
	std::vector<glm::vec3> colors; //fill with constructor  

	std::vector<Edge> edges; //fill with constructor 
	std::vector<int> vertex_indices; //fill with constructor 
	std::vector<std::vector<std::pair<int, float>>> adjacenies; // all adjacensies for all vertices
	std::vector<unsigned int > triangles;
	std::vector<std::vector<unsigned int>> neighbours; 
	std::vector<std::pair<unsigned int, unsigned int>> symmetry_pairs; // for ground truth 
	std::vector<unsigned int> symmetry_pairs_map; // for ground truth hashmap for pairs
	std::vector<std::pair<unsigned int, unsigned int>> calculated_symmetry_pairs; // for what we calculated
	std::vector<glm::vec3> normals; // fil with constructor, actual vertices
	std::vector<float> normals_display; // vertex information of normals
	glm::mat4 model_mat;
	std::string file_name; 
	//VAO object
	unsigned int vao_normals; 
	//project 
	bool off_format = false;
	bool ply_format = false;
	//MVP
	glm::mat4 MVP;
	Mesh();
	Mesh(char* filename);
	Mesh(glm::vec3 *p1 , glm::vec3* p2 , glm::vec3* p3 , glm::vec3* p4);
	glm::mat4 move_mesh(glm::vec3 direction);
	glm::mat4 scale_mesh(glm::vec3 scale);

	private:
	void read_ply_format(char* filename);
	void read_off_format(char* filename);
	void generate_normals();
};

void read_symmetry_format(char* filename, Mesh* m);
