#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include "glm/glm.hpp"
#include "raylib.h"
using std::ofstream;
using std::ifstream;
#define INFINITY 10000000

struct Edge
{
	glm::vec3 p1;
	glm::vec3 p2;
	float distance;
};
class TrilateralMesh
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
	std::vector<unsigned int> ground_truth_symmetry_pairs; // for ground truth hashmap for pairs
	std::vector<std::pair<unsigned int, unsigned int>> calculated_symmetry_pairs; // for what we calculated
	std::vector<glm::vec3> normals; // fil with constructor, actual vertices
	std::vector<float> normals_display; // vertex information of normals
	std::vector<float> areas; // vertex information of normals
	std::vector<float> normalized_heat_kernel_signature; // vertex information of normals
	glm::mat4 model_mat;
	std::string file_name;
	Mesh raylib_mesh;
	//VAO object
	unsigned int vao_normals; 
	//project 
	bool off_format = false;
	bool ply_format = false;
	//MVP
	glm::mat4 MVP;
	//area
	float mesh_area; 
	TrilateralMesh();
	TrilateralMesh(char* filename);
	TrilateralMesh(glm::vec3 *p1 , glm::vec3* p2 , glm::vec3* p3 , glm::vec3* p4);
	glm::mat4 move_mesh(glm::vec3 direction);
	glm::mat4 scale_mesh(glm::vec3 scale);

	private:
	void read_ply_format(char* filename);
	void read_off_format(char* filename);
	void generate_normals();
};
typedef struct
{
	std::vector<float> skeleton_points;
	std::vector<unsigned int> skeleton_indices;
}SkeletonMesh;
void read_symmetry_format(char* filename, TrilateralMesh* m);
unsigned int mesh_get_closest_index(TrilateralMesh* m, const glm::vec3& point);
glm::vec3 mesh_generate_weighted_mid_point(TrilateralMesh* m); //supposed to be best way 
std::vector<float> mesh_point_surfel_normalized(TrilateralMesh*m );