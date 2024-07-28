#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <vector>
#include <glm/ext/vector_float3.hpp>
#include <glm/geometric.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "happly.h"
#include "../Include/Mesh.h"


#define INFINITE 10000000


#pragma region functions


typedef struct
{
	std::vector<float> point_pairs;
	glm::mat4 model_mat;
}MeshPointPairs;
#pragma endregion functions
class MeshFactory
{
public:
	MeshFactory();
	~MeshFactory();
	void add_mesh(Mesh &m);
	void buffer_meshes();
	void draw_meshes();
	void draw_mesh(int mesh_no);
	void draw_lines();
	void create_lines(int mesh_no1 , int mesh_no2 , int  no_of_lines);
	void remove_mesh(int mesh_no);
	void remove_all();
	void add_all();
	void add_lines_between_meshes(int mesh_index1, int mesh_index2, int p1, int p2);
	void get_camera_and_projection(glm::mat4 view_, glm::mat4  projection_);
	void remove_lines();
	std::vector<Mesh> mesh_vec;
	std::vector<MeshPointPairs> mesh_point_pairs; //point paris for each mesh
	std::vector<glm::vec3> points;
	std::vector<glm::vec3> colors;
	std::vector<int> point_indices;
	std::vector<std::pair<int, int>> correspondence_lines;
	//skeleton part 

	
	
	// for skeleton detection from pose  video-to-pose
	SkeletonMesh mesh_skeleton_vec;
	unsigned int skeleton_VAO;  // skeleton

											 // 
	glm::mat4 view;
	glm::mat4 projection;
	// line 
private:
	
};

