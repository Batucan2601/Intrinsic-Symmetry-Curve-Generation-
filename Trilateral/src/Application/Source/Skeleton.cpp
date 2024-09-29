#include "../Include/Skeleton.h"
#include "../Include/Mesh.h"
#include "../Include/ShapeDiameter.h"
#include "../Include/MeshFactory.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <GL/glew.h>
#include <unordered_set>
#include <eigen/Eigen/Dense>

#pragma region VIDEO-TO-POSE 3D repo functions legacy for now
using std::ifstream;
std::vector<float> generate_bounding_box(std::string file_name)
{
	std::vector<float> bounding_box(6,0); 
	int bounding_box_index = 0;
	int line_count = 0;
	ifstream indata; // indata is like cin
	
	indata.open("../../Trilateral/Mesh/off/KIDS_skeleton/" + file_name);
	if (!indata)
	{
		return bounding_box;
	}
	std::string line;
	while (std::getline(indata, line))
	{
		if (line_count > 31 /*16*/)
		{
			//read the values
			std::stringstream ss(line.substr(6));
			
			if (!(ss >> bounding_box[bounding_box_index++]))
			{
				int err = 1;
			}
		}
		line_count += 1;


	}
	indata.close();

	return bounding_box;
}

std::map<std::string, glm::vec3 > generate_skeleton_keypoints(std::string file_name)
{

	std::map<std::string, glm::vec3 > keypoints;
	int bounding_box_index = 0;
	int line_count = 0;
	ifstream indata; // indata is like cin

	indata.open("../../Trilateral/Mesh/off/KIDS_skeleton/" + file_name);
	if (!indata)
	{
		return keypoints;
	}
	std::string line;
	while (std::getline(indata, line))
	{
		if (line_count <= 16)
		{
			//read the values
			auto two_dot_pos = line.find(':');

			std::string key = line.substr(0, two_dot_pos);

			auto square_brackets_start_pos = line.find('[');
			auto square_brackets_end_pos = line.find(']');

			std::string val = line.substr(square_brackets_start_pos+1, square_brackets_end_pos-1);

			std::stringstream ss(val);
			glm::vec3 point; 
			ss >> point.x;
			ss >> point.y;
			ss >> point.z;

			keypoints[key] = point; 

		}
		else
		{
			break; 
		}
		line_count += 1;


	}
	indata.close();

	return keypoints;
}

void match_skeleton_keypoints( MeshFactory& meshFactory ,Mesh* m , std::vector<float>& skeleton_bounding_box, std::map<std::string, glm::vec3>& keypoints)
{
	// 1- find the box in world space
	std::vector<float> bb_mesh;
	/*
	1 - min_x
	2 - min_y
	3 - min_z
	4 - max_x
	5 - max_y
	6 - max_z
	*/
	bb_mesh.push_back(INFINITY);
	bb_mesh.push_back(INFINITY);
	bb_mesh.push_back(INFINITY);
	bb_mesh.push_back(-INFINITY);
	bb_mesh.push_back(-INFINITY);
	bb_mesh.push_back(-INFINITY);
	
	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (m->vertices[i].x < bb_mesh[0]) //xmin 
		{
			bb_mesh[0] = m->vertices[i].x;
		}
		if (m->vertices[i].y < bb_mesh[1]) //ymin 
		{
			bb_mesh[1] = m->vertices[i].y;
		}
		if (m->vertices[i].z < bb_mesh[2]) //zmin
		{
			bb_mesh[2] = m->vertices[i].z;
		}
		if (m->vertices[i].x > bb_mesh[3]) //xmax
		{
			bb_mesh[3] = m->vertices[i].x;
		}
		if (m->vertices[i].y > bb_mesh[4]) //ymax
		{
			bb_mesh[4] = m->vertices[i].y;
		}
		if (m->vertices[i].z > bb_mesh[5]) //zmax
		{
			bb_mesh[5] = m->vertices[i].z;

		}
	}
	std::cout << " original mesh bb " << std::endl;
	for (size_t i = 0; i < 6 ; i++)
	{
		std::cout << bb_mesh[i] << std::endl;
	}
	std::cout << " sekeleton bb " << std::endl;
	for (size_t i = 0; i < 6; i++)
	{
		std::cout << skeleton_bounding_box[i] << std::endl;
	}

	//assume we start from lower left corner (x_min y_min z_min)
	glm::vec3 start_point(skeleton_bounding_box[0], skeleton_bounding_box[1], skeleton_bounding_box[2]);

	// converts 1 unit of skeleton's object space to world space in the application I guess
	glm::vec3 skeleton_to_world_space_axis((bb_mesh[3] - bb_mesh[0]) / (skeleton_bounding_box[0] - skeleton_bounding_box[3]),
		(bb_mesh[4] - bb_mesh[1]) / (skeleton_bounding_box[1] - skeleton_bounding_box[4]),
		(bb_mesh[5] - bb_mesh[2]) / (skeleton_bounding_box[5] - skeleton_bounding_box[2]));

	

	//find the smallest dif to point and select it as l ankle

	std::map<std::string, glm::vec3>::iterator it = keypoints.begin();


	// change of axis 

	// 1 - swap -x and z 
	/*it = keypoints.begin();
	while (it != keypoints.end())
	{
		glm::vec3 temp = it->second;

		float temp_f = temp.y;
		temp.y = temp.z;
		temp.z = temp_f;
		it->second = temp;
		it++;
	}*/

	it = keypoints.begin();
	while (it != keypoints.end())
	{
		glm::vec3 key_point_pos = keypoints[it->first];

		glm::vec3 l_keypoint_subtract_minimum = key_point_pos - start_point;

		glm::vec3 pos = glm::vec3(l_keypoint_subtract_minimum.x * skeleton_to_world_space_axis.x, l_keypoint_subtract_minimum.y * skeleton_to_world_space_axis.y,
			l_keypoint_subtract_minimum.z * skeleton_to_world_space_axis.z);

		pos = pos + glm::vec3(bb_mesh[0], bb_mesh[1], bb_mesh[2]);
		float min_dif = INFINITY;
		int index = -1;
		for (size_t i = 0; i < m->vertices.size(); i++)
		{
			if (glm::distance(pos, m->vertices[i]) < min_dif)
			{
				min_dif = glm::distance(pos, m->vertices[i]);
				index = i;
			}
		}
		m->colors[index] = glm::vec3(255, 0, 0);
		std::cout << it->first << "  " <<  m->vertices[index].x << "  " << m->vertices[index].y << "  " << m->vertices[index].z << std::endl;

		it++;
	}

	// generate the skeleton

	// 1 -between eyes 
	/*meshFactory.mesh_skeleton_vec.push_back(keypoints["LEye"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LEye"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LEye"].z);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["REye"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["REye"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["REye"].z);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	// 2 -between eyes and Nose 

	meshFactory.mesh_skeleton_vec.push_back((keypoints["LEye"].x + keypoints["REye"].x)  / 2.0f );
	meshFactory.mesh_skeleton_vec.push_back((keypoints["LEye"].y + keypoints["REye"].y)  / 2.0f );
	meshFactory.mesh_skeleton_vec.push_back((keypoints["LEye"].z + keypoints["REye"].z)  / 2.0f );

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["Nose"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	// 3 -between Shoulders 
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	// 4 - Full arm 
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LElbow"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LElbow"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LElbow"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LElbow"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LElbow"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LElbow"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LWrist"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LWrist"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LWrist"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RElbow"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RElbow"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RElbow"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RElbow"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RElbow"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RElbow"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RWrist"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RWrist"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RWrist"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	// 5 - Between shoulders and between hips 

	meshFactory.mesh_skeleton_vec.push_back((keypoints["RShoulder"].x + keypoints["LShoulder"].x) / 2.0f);
	meshFactory.mesh_skeleton_vec.push_back((keypoints["RShoulder"].y + keypoints["LShoulder"].y) / 2.0f);
	meshFactory.mesh_skeleton_vec.push_back((keypoints["RShoulder"].z + keypoints["LShoulder"].z) / 2.0f);

	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back((keypoints["RHip"].x + keypoints["LHip"].x) / 2.0f);
	meshFactory.mesh_skeleton_vec.push_back((keypoints["RHip"].y + keypoints["LHip"].y) / 2.0f);
	meshFactory.mesh_skeleton_vec.push_back((keypoints["RHip"].z + keypoints["LHip"].z) / 2.0f);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	// 6 - HIPS

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RHip"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RHip"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RHip"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LHip"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LHip"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LHip"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	// full legs

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RHip"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RHip"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RHip"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RKnee"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RKnee"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RKnee"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RKnee"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RKnee"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RKnee"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RAnkle"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RAnkle"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RAnkle"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);


	meshFactory.mesh_skeleton_vec.push_back(keypoints["LHip"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LHip"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LHip"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LKnee"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LKnee"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LKnee"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LKnee"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LKnee"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LKnee"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["LAnkle"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LAnkle"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LAnkle"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(255.0f);*/


	skeleton_generate_buffer(meshFactory);
	skeleton_buffer(meshFactory);
	/*unsigned int skeleton_vao, skeleton_vbo;
	glGenVertexArrays(1, &skeleton_vao);
	glGenBuffers(1, &skeleton_vbo);

	glBindVertexArray(skeleton_vao);
	glBindBuffer(GL_ARRAY_BUFFER, skeleton_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glBufferData(GL_ARRAY_BUFFER, meshFactory.mesh_skeleton_vec.size() * sizeof(float), &meshFactory.mesh_skeleton_vec[0], GL_STATIC_DRAW);

	meshFactory.skeleton_VAO = skeleton_vao;

	glBindVertexArray(0);*/

}

std::vector<float> generate_skeleton_lines(std::string file_name)
{
	std::vector<float> skeleton_lines;
	int line_count = 0;
	ifstream indata; // indata is like cin
	indata.open("../../Trilateral/Mesh/off/KIDS_skeleton/" + file_name);
	if (!indata)
	{
		return skeleton_lines;
	}
	std::string line;
	while (std::getline(indata, line))
	{
		if (line_count <= 32)
		{
			//read the values


			std::stringstream ss(line);
			glm::vec3 point;
			ss >> point.x;
			ss >> point.y;
			ss >> point.z;

			skeleton_lines.push_back(  point.x);
			skeleton_lines.push_back(  point.y);
			skeleton_lines.push_back(  point.z);

			skeleton_lines.push_back(255.0f);
			skeleton_lines.push_back(255.0f);
			skeleton_lines.push_back(255.0f);

		}
		else
		{
			break;
		}
		line_count += 1;


	}
	indata.close();

	return skeleton_lines;
}

void match_skeleton_lines(MeshFactory& meshFactory, Mesh* m, std::vector<float>& skeleton_bounding_box, std::vector<float> skeleton_lines)
{
	meshFactory.mesh_skeleton_vec.skeleton_points = skeleton_lines;

	std::vector<float> bb_mesh;
	/*
	1 - min_x
	2 - min_y
	3 - min_z
	4 - max_x
	5 - max_y
	6 - max_z
	*/
	bb_mesh.push_back(INFINITY);
	bb_mesh.push_back(INFINITY);
	bb_mesh.push_back(INFINITY);
	bb_mesh.push_back(-INFINITY);
	bb_mesh.push_back(-INFINITY);
	bb_mesh.push_back(-INFINITY);

	for (size_t i = 0; i < m->vertices.size(); i++)
	{
		if (m->vertices[i].x < bb_mesh[0]) //xmin 
		{
			bb_mesh[0] = m->vertices[i].x;
		}
		if (m->vertices[i].y < bb_mesh[1]) //ymin 
		{
			bb_mesh[1] = m->vertices[i].y;
		}
		if (m->vertices[i].z < bb_mesh[2]) //zmin
		{
			bb_mesh[2] = m->vertices[i].z;
		}
		if (m->vertices[i].x > bb_mesh[3]) //xmax
		{
			bb_mesh[3] = m->vertices[i].x;
		}
		if (m->vertices[i].y > bb_mesh[4]) //ymax
		{
			bb_mesh[4] = m->vertices[i].y;
		}
		if (m->vertices[i].z > bb_mesh[5]) //zmax
		{
			bb_mesh[5] = m->vertices[i].z;

		}
	}
	std::cout << " original mesh bb " << std::endl;
	for (size_t i = 0; i < 6; i++)
	{
		std::cout << bb_mesh[i] << std::endl;
	}
	std::cout << " sekeleton bb " << std::endl;
	for (size_t i = 0; i < 6; i++)
	{
		std::cout << skeleton_bounding_box[i] << std::endl;
	}

	//assume we start from lower left corner (x_min y_min z_min)
	glm::vec3 start_point(skeleton_bounding_box[0], skeleton_bounding_box[1], skeleton_bounding_box[2]);

	//swap axises
	// x and z 
	/*for (size_t i = 0; i < skeleton_lines.size(); i+=6)
	{
		float temp = skeleton_lines[i];
		skeleton_lines[i]  = skeleton_lines[i + 2 ];
		skeleton_lines[i + 2] = temp; 
	}*/


	// converts 1 unit of skeleton's object space to world space in the application I guess
	glm::vec3 skeleton_to_world_space_axis((bb_mesh[3] - bb_mesh[0]) / (skeleton_bounding_box[3] - skeleton_bounding_box[0]),
		(bb_mesh[4] - bb_mesh[1]) / (skeleton_bounding_box[4] - skeleton_bounding_box[1]),
		(bb_mesh[5] - bb_mesh[2]) / (skeleton_bounding_box[5] - skeleton_bounding_box[2]));

	
	

	for( int i =0; i < skeleton_lines.size(); i += 6   ) // 3 pos 3 color 
	{
		glm::vec3 key_point_pos = glm::vec3( skeleton_lines[i], skeleton_lines[i+1], skeleton_lines[i+2]);

		glm::vec3 l_keypoint_subtract_minimum = key_point_pos - start_point;

		glm::vec3 pos = glm::vec3(l_keypoint_subtract_minimum.x * skeleton_to_world_space_axis.x, l_keypoint_subtract_minimum.y * skeleton_to_world_space_axis.y,
			l_keypoint_subtract_minimum.z * skeleton_to_world_space_axis.z);

		pos = pos + glm::vec3(bb_mesh[0], bb_mesh[1], bb_mesh[2]);
		//float min_dif = INFINITY;
		//int index = -1;
		//for (size_t j = 0; j < m->vertices.size(); j++)
		//{
		//	glm::vec3 offsetted_vertex = m->vertices[j] - glm::vec3(bb_mesh[0], bb_mesh[1], bb_mesh[2]);

		//	if (glm::distance(pos, /*m->vertices[j]*/ offsetted_vertex) < min_dif)
		//	{
		//		min_dif = glm::distance(pos,/* m->vertices[j]*/ offsetted_vertex);
		//		index = j;
		//	}
		//}
		//m->colors[index] = glm::vec3(255, 0, 0);
	}

	std::vector<float> offsetted_skeleton_lines; 
	for (size_t i = 0; i < skeleton_lines.size(); i+=6)
	{
		glm::vec3 key_point_pos = glm::vec3(skeleton_lines[i], skeleton_lines[i + 1], skeleton_lines[i + 2]);

		start_point = glm::vec3((skeleton_bounding_box[0] + skeleton_bounding_box[3]) / 2, (skeleton_bounding_box[1] + skeleton_bounding_box[4]) / 2, (skeleton_bounding_box[2] + skeleton_bounding_box[5]) / 2);

		glm::vec3 skeleton_line_offsetted = key_point_pos - start_point;

		glm::vec3 pos = glm::vec3(skeleton_line_offsetted.x * skeleton_to_world_space_axis.x, skeleton_line_offsetted.y * skeleton_to_world_space_axis.y,
			skeleton_line_offsetted.z * skeleton_to_world_space_axis.z);

		//new approach
		//glm::vec3 pos = skeleton_line_offsetted * (glm::length(skeleton_to_world_space_axis));
		offsetted_skeleton_lines.push_back(pos.x);
		offsetted_skeleton_lines.push_back(pos.y);
		offsetted_skeleton_lines.push_back(pos.z);

		offsetted_skeleton_lines.push_back(255);
		offsetted_skeleton_lines.push_back(0);
		offsetted_skeleton_lines.push_back(0);

		float min_dif = INFINITY;
		int index = -1;
		for (size_t j = 0; j < m->vertices.size(); j++)
		{
			glm::vec3 offsetted_vertex = m->vertices[j] - glm::vec3((bb_mesh[0] + bb_mesh[3] ) / 2 , (bb_mesh[1] + bb_mesh[4]) / 2, (bb_mesh[2] + bb_mesh[5]) / 2);

			if (glm::distance(pos, /*m->vertices[j]*/ offsetted_vertex) < min_dif)
			{
				min_dif = glm::distance(pos,/* m->vertices[j]*/ offsetted_vertex);
				index = j;
			}
		}
		m->colors[index] = glm::vec3(255, 0, 0);

	}

	unsigned int skeleton_vao, skeleton_vbo;
	glGenVertexArrays(1, &skeleton_vao);
	glGenBuffers(1, &skeleton_vbo);

	glBindVertexArray(skeleton_vao);
	glBindBuffer(GL_ARRAY_BUFFER, skeleton_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glBufferData(GL_ARRAY_BUFFER, meshFactory.mesh_skeleton_vec.skeleton_points.size() * sizeof(float), /*&meshFactory.mesh_skeleton_vec[0]*/&offsetted_skeleton_lines[0], GL_STATIC_DRAW);

	meshFactory.skeleton_VAO = skeleton_vao;

	glBindVertexArray(0);


}
#pragma endregion VIDEO-TO-POSE 3D repo functions legacy for now




Skeleton skeleton_read_swc_file(MeshFactory& meshFactory,std::string file_name)
{
	Skeleton skeleton; 
	std::vector<SkeletonFormat> skeletonPoints;
	ifstream indata; // indata is like cin
	indata.open("../../Trilateral/Mesh/off/KIDS_skeleton/" + file_name);
	if (!indata)
	{
		return skeleton;
	}
	std::string line;
	unsigned int point_index = 0;
	while (std::getline(indata, line))
	{
		if (line.find('#') == std::string::npos && line.size() > 1 )
		{
			SkeletonFormat point; 
			SkeletonEndPoint skeleton_end_point; 
			float point_no;
			POINT_LABEL label_point;
			int soma_point;
			int parent; 
			float radius; 
			//read the points 
			std::stringstream ss(line);

			ss >> point_no;
			ss >> soma_point;

			ss >> point.point.x;
			ss >> point.point.y;
			ss >> point.point.z;

			ss >> radius; 

			ss >> point.parent;


			point.label = (POINT_LABEL)soma_point;
			skeletonPoints.push_back(point);

			if (point.label == END)
			{
				skeleton_end_point.index = point_index;
				skeleton.endPoints.push_back(skeleton_end_point);
			}
			
			point_index++;
		}

	}

	//adjacencies
	// just to be sure make them all false
	std::vector<std::vector<unsigned int>> adjacencies(skeletonPoints.size());
 	for (size_t i = 0; i < skeletonPoints.size(); i++)
	{
		for (size_t j = 0; j < skeletonPoints.size(); j++)
		{
			if (skeletonPoints[j].parent == i )
			{
				adjacencies[i].push_back(j);
				adjacencies[j].push_back(i);

			}
		}
	}

	std::vector<float> skeleton_points; 
	std::vector<unsigned int> skeleton_indices; 
	// another pass needed to generate lines with according parent
	for (size_t i = 0; i < skeletonPoints.size(); i++)
	{
		glm::vec3 point = skeletonPoints[i].point;
		int parent = skeletonPoints[i].parent;

		skeleton_points.push_back(point.x);
		skeleton_points.push_back(point.y);
		skeleton_points.push_back(point.z);
			
		if (skeletonPoints[i].label == END)
		{
			skeleton_points.push_back(255.0f);
			skeleton_points.push_back(0.0f);
			skeleton_points.push_back(0.0f);
		}
		else if (skeletonPoints[i].parent == -1)
		{
			skeleton_points.push_back(255.0f);
			skeleton_points.push_back(255.0f);
			skeleton_points.push_back(255.0f);
		}
		else
		{
			skeleton_points.push_back(255.0f);
			skeleton_points.push_back(255.0f);
			skeleton_points.push_back(255.0f);
		}

		if (parent >= 0)
		{
			skeleton_indices.push_back(i);
			skeleton_indices.push_back(parent);
		}

	}
	

	
	//lastly get midpoint and get the closest vertex
	glm::vec3 mid_point(0.0f, 0.0f, 0.0f);
	Mesh* m = &meshFactory.mesh_vec[0];
	mid_point =  mesh_generate_weighted_mid_point(m);

	//check the closest vertex in skeleton
	float minimum_dist = INFINITY;
	int minimum_index = -1;
	for (size_t i = 0; i < skeletonPoints.size(); i++)
	{
		float dist = glm::distance(skeletonPoints[i].point, mid_point);
		if ( minimum_dist > dist)
		{
			minimum_dist = dist; 
			minimum_index = i;
		}
	}


	//color mid point white
	skeleton_points[minimum_index * 6 + 3] = 0.0f;
	skeleton_points[minimum_index * 6 + 4] = 255.0f;
	skeleton_points[minimum_index * 6 + 5] = 0.0f;


	skeleton.skeletonFormat = skeletonPoints;
	skeleton.adjacencies = adjacencies;
	skeleton.mid_point_index = minimum_index;
	skeleton.skeleton_mid_point = mid_point;


	meshFactory.mesh_skeleton_vec.skeleton_points = skeleton_points;
	meshFactory.mesh_skeleton_vec.skeleton_indices = skeleton_indices;

	skeleton_generate_buffer(meshFactory);
	skeleton_buffer(meshFactory);

	return skeleton;
}


void skeleton_calculate_dijkstra(Skeleton skeleton, int index1,
	std::vector<int>& vertex_list, std::vector<float>& dijkstra_distances)
{
	int N = skeleton.skeletonFormat.size();
	//run a dijkstra from index1 
	std::vector<float> distances(N, INFINITY);
	distances[index1] = 0.0f;
	
	std::vector<int> discovered_vertices;
	std::vector<bool> is_vertex_discovered( N , false);
	std::vector<int> predecessors(N, -1);


	// initila nei
	discovered_vertices.push_back(index1);
	is_vertex_discovered[index1] = true;
	int no_of_discovered = 1;
	while (no_of_discovered < N )
	{
		float minimum_distance = INFINITY;
		int minimum_distance_index = -1;
		int discovered_vertex = -1; 
		//check adjacent non_discovered points
		for (size_t i = 0; i < discovered_vertices.size(); i++)
		{
			int discovered_vertex_adjacency_size = skeleton.adjacencies[discovered_vertices[i]].size();
			//check adjacencies
			for (size_t j = 0; j < discovered_vertex_adjacency_size; j++)
			{
				int discovered_vertex_adjaceny = skeleton.adjacencies[discovered_vertices[i]][j];
				if( !is_vertex_discovered[discovered_vertex_adjaceny])
				{
					//check distance and compare
					float dist = glm::distance(skeleton.skeletonFormat[discovered_vertices[i]].point, skeleton.skeletonFormat[discovered_vertex_adjaceny].point);
					if (dist < minimum_distance)
					{
						minimum_distance = dist; 
						minimum_distance_index = discovered_vertex_adjaceny;
						discovered_vertex = discovered_vertices[i];
					}
				}
			}
		}

		distances[minimum_distance_index] = distances[discovered_vertex] + minimum_distance;;
		predecessors[minimum_distance_index] = discovered_vertex;

		is_vertex_discovered[minimum_distance_index] = true;	
		discovered_vertices.push_back(minimum_distance_index);

		//change the weights
		for (size_t i = 0; i < skeleton.adjacencies[minimum_distance_index].size(); i++)
		{
			//if (is_vertex_discovered[skeleton.adjacencies[minimum_distance_index][i]])
			{
				glm::vec3 p_index1 = skeleton.skeletonFormat[minimum_distance_index].point;
				glm::vec3 p = skeleton.skeletonFormat[skeleton.adjacencies[minimum_distance_index][i]].point;
				float edge_dist = glm::distance(p_index1, p);
				if ((edge_dist + distances[minimum_distance_index]) < distances[skeleton.adjacencies[minimum_distance_index][i]])
				{
					distances[skeleton.adjacencies[minimum_distance_index][i]] = (edge_dist + distances[minimum_distance_index]);
					predecessors[skeleton.adjacencies[minimum_distance_index][i]] = minimum_distance_index;
				}
				else
				{
					int a = 1;
				}
			}
			
		}
		no_of_discovered += 1;
	}
	vertex_list = predecessors;
	dijkstra_distances = distances;
	return; 
}

//std::vector<std::vector<float>> skeleton_distances_table(std::vector<SkeletonFormat> skeletonFormat)
//{
//	// run dijkstra
//	
//	//dijkstra table
//	std::vector<std::vector<float>> distances_table;
//	//create an adjacency array
//	std::vector<std::vector<bool>> adjacencies;
//	//create an distances array
//	std::vector<std::vector<float>> distances;
//
//	//init dijkstra table
//	for (size_t i = 0; i < skeletonFormat.size(); i++)
//	{
//		std::vector<float> distance_table_i(skeletonFormat.size() , INFINITY);
//		distance_table_i[i] = 0;
//		distances_table.push_back(distance_table_i);
//	}
//	// just to be sure make them all false
//	for (size_t i = 0; i < skeletonFormat.size(); i++)
//	{
//		std::vector<bool> adjacency_i;
//		for (size_t j = 0; j < skeletonFormat.size(); j++)
//		{
//			adjacency_i.push_back(false);
//		}
//		adjacencies.push_back(adjacency_i);
//	}
//
//	//fill adjacencies
//	for (size_t i = 0; i < skeletonFormat.size(); i++)
//	{
//		int parent_index = skeletonFormat[i].parent;
//		adjacencies[i][parent_index] = true; 
//		adjacencies[parent_index][i] = true; 
//	}
//	//fill distances
//	for (size_t i = 0; i < skeletonFormat.size(); i++)
//	{
//		std::vector<float> distances_i;
//		glm::vec3 pi = skeletonFormat[i].point;
//		for (size_t j = 0; j < skeletonFormat.size(); j++)
//		{
//			if (adjacencies[i][j])
//			{
//				glm::vec3 pj = skeletonFormat[skeletonFormat[i].parent].point;
//				float distance = glm::distance(pi,pj);
//				distances_i.push_back(distance);
//			}
//			else
//			{
//				distances_i.push_back(INFINITY);
//
//			}
//		}
//	}
//	
//	//step 2 do dijkstra for ALL vertices
//	for (size_t i = 0; i < skeletonFormat.size(); i++)
//	{
//		std::vector<bool> is_explored(skeletonFormat.size(), false);
//		is_explored[i] = true;
//		int number_of_explored_vertices = 1;
//		while (number_of_explored_vertices < skeletonFormat.size())
//		{
//			//select neighbours
//			for (size_t j = 0; j < skeletonFormat.size(); j++)
//			{
//				if (adjacencies[i][j]) // if adjacent
//				{
//					if (!is_explored[j])
//					{
//						is_explored[j] = true;
//						distances_table[i][j] = glm::distance(skeletonFormat[i].point, skeletonFormat[j].point);
//
//					}
//				}
//			}
//			number_of_explored_vertices++;
//		}
//	}
//	
//}

void skeleton_generate_backbone(MeshFactory& meshFac, Skeleton skeleton, unsigned int mesh_index,BackBone& best_backbone ,
std::vector<unsigned int>& best_right_points , std::vector<unsigned int>& best_left_points )
{
	Mesh* mesh = &meshFac.mesh_vec[mesh_index];
	int N = skeleton.skeletonFormat.size();
	std::vector<unsigned int> end_point_indices;
	for (size_t i = 0; i < N; i++)
	{
		if (skeleton.skeletonFormat[i].label == END)
		{
			end_point_indices.push_back(i);
		}
	}
	int N_end_points = end_point_indices.size();
	std::vector<std::vector<int>> predecessor_list_for_end_points(N_end_points);
	std::vector<std::vector<float>> distance_matrix(N_end_points);


	//before generating please calculate  the sdf
	std::vector<unsigned int>mesh_index_vertices; 
	std::vector<float> shape_diameter_values;
	skeleton_calculate_closest_mesh_points(skeleton, &meshFac.mesh_vec[mesh_index], mesh_index_vertices);
	ShapeDiameter_calculate(&meshFac.mesh_vec[mesh_index], mesh_index_vertices, shape_diameter_values);
	// brute force generating backbone
	std::vector<BackBone> candidate_backbones; 
	for (size_t i = 0; i < N_end_points; i++)
	{
		int index1 = end_point_indices[i];
		std::vector<int> point_list;
		std::vector<float> distances_from_index1;

		skeleton_calculate_dijkstra(skeleton, index1, point_list, distances_from_index1);

		for (size_t j = i; j < N_end_points; j++)
		{
			//assume i and j are the end points of backbone

			//generate the path
			int index2 = end_point_indices[j];
			float backbone_length;
			std::vector<int> vertex_list_index1_index2;

			if (i != j)
			{ 
				skeleton_get_distance_and_vertex_list(skeleton, index1, index2, point_list,
					vertex_list_index1_index2, backbone_length);

				BackBone backbone; 
				backbone.start_index = index1;
				backbone.end_index   = index2;
				backbone.vertex_list = vertex_list_index1_index2;

				candidate_backbones.push_back(backbone);
				
			}
			predecessor_list_for_end_points[i] = point_list;

			distance_matrix[i] = distances_from_index1;

		}
	}
	std::vector<float> backbone_affinity_diffs; 
	std::vector<std::vector<std::pair<unsigned int, unsigned int>>> backbone_pairs_vec;
	std::vector<std::vector<unsigned int>> backbone_right_points;
	std::vector<std::vector<unsigned int>> backbone_left_points;
	//with each different backbone separate the endpoints into two
	for (size_t i = 0; i < candidate_backbones.size(); i++)
	{
		std::vector<unsigned int> right_points;
		std::vector<unsigned int> left_points;
		std::vector<NodeAffinityParams> right_points_node_params; 
		std::vector<NodeAffinityParams> left_points_node_params;


		//get the middle point of backbone
		for (size_t j = 0; j < N_end_points; j++)
		{
			//if not consisting in backbone 
			if (end_point_indices[j] != candidate_backbones[i].start_index && end_point_indices[j] != candidate_backbones[i].end_index)
			{
				int hitIndex = -1;
				float dist = -1.0f;
				std::vector<int> indices;
				skeleton_point_to_backbone(skeleton, candidate_backbones[i], end_point_indices[j], hitIndex, dist , indices , distance_matrix[j],
					predecessor_list_for_end_points[j]);

				//lets cross these vectors 
				// vector1 -> vector of hitIndex and it's parent on backbone chain
				// vector2 -> the n-2'th index on indices vectpr subtracted by n-1'th

				glm::vec3 hitIndexPoint(skeleton.skeletonFormat[hitIndex].point);
				glm::vec3 hitIndexBackBonePredecessor;

				// get the predecessor
				for (size_t k = 0; k < candidate_backbones[i].vertex_list.size(); k++)
				{
					if (candidate_backbones[i].vertex_list[k] == hitIndex)
					{
						hitIndexBackBonePredecessor = glm::vec3(skeleton.skeletonFormat[candidate_backbones[i].vertex_list[k - 1]].point);
						break;
					}
				}
				glm::vec3 vec1 = hitIndexBackBonePredecessor - hitIndexPoint;



				NodeAffinityParams param;
				param.distance_to_backbone = dist;
				param.point_in_backbone = hitIndex;
#pragma region WRONG
				//https://math.stackexchange.com/questions/214187/point-on-the-left-or-right-side-of-a-plane-in-3d-space#:~:text=How%20do%20you%20define%20right%2Fleft%3F&text=To%20distinguish%20the%20two%20sides,%E2%80%93%20J.%20J.


				if (candidate_backbones.size() <= 2)
				{
					continue;
				}

				int right_count = 0;
				int left_count = 0;
				for (size_t k = 1; k < candidate_backbones[i].vertex_list.size() - 1 ; k++)
				{
					glm::vec3 A = skeleton.skeletonFormat[candidate_backbones[i].vertex_list[k-1]].point;
					glm::vec3 B = skeleton.skeletonFormat[candidate_backbones[i].vertex_list[k]].point;
					glm::vec3 C = skeleton.skeletonFormat[candidate_backbones[i].vertex_list[k + 1]].point;
					glm::vec3 X = mesh->vertices[end_point_indices[j]];

					glm::vec3 A_ = A - B;
					glm::vec3 C_ = C - B;
					glm::vec3 X_ = X - B;
					Eigen::Matrix3f sign_matrix;
					sign_matrix.row(0) = Eigen::Vector3f(A_.x, A_.y, A_.z);
					sign_matrix.row(1) = Eigen::Vector3f(C_.x, C_.y, C_.z);
					sign_matrix.row(2) = Eigen::Vector3f(X_.x, X_.y, X_.z);

					float det = sign_matrix.determinant();

					if (det >= 0)
					{
						right_count++;
					}
					else
					{
						left_count++;
					}
				}

			

				if (right_count >= left_count)
				{
					right_points.push_back(mesh_index_vertices[j]);
					right_points_node_params.push_back(param);
				}
				else
				{
					left_points.push_back(mesh_index_vertices[j]);
					left_points_node_params.push_back(param);
				}
 				//if (glm::dot(vec1, glm::normalize(skeleton.skeletonFormat[0].point - skeleton.skeletonFormat[1].point)) > 0)
				//{
				//	right_points.push_back(mesh_index_vertices[j]);
				//	right_points_node_params.push_back(param);
				//}
				//else
				//{
				//	left_points.push_back(mesh_index_vertices[j]);
				//	left_points_node_params.push_back(param);
				//}
#pragma endregion WRONG 
				

			}
		}

		
		if (right_points.size() == 0 || left_points.size() == 0)
		{
			continue;
		}
		//now calculate the node affinity Skeleton-IntrinsicSymmetrizationofShapes
		std::vector<glm::vec3> hit_points;
		std::vector<std::pair<unsigned int , unsigned int >> backbone_pairs;
		// 1-  db(i,i_ )
		// 2- dl

		if (left_points_node_params.size() > right_points_node_params.size())
		{
			std::vector<NodeAffinityParams> temp_param_vec;
			std::vector<unsigned int > temp_pair_vec;

			temp_param_vec = right_points_node_params;
			right_points_node_params = left_points_node_params;
			left_points_node_params = temp_param_vec;

			temp_pair_vec = right_points;
			right_points = left_points;
			left_points = temp_pair_vec;
		}
		int distance_matrix_start_index = -1;
		int distance_matrix_end_index = -1;
		for (size_t t = 0; t < end_point_indices.size(); t++)
		{
			if (candidate_backbones[i].start_index == end_point_indices[t])
			{
				distance_matrix_start_index = t;
			}
			else if (candidate_backbones[i].end_index == end_point_indices[t])
			{
				distance_matrix_end_index = t;
			}
		}
		float total_dif_on_backbone = 0; 
		for (size_t j = 0; j < right_points_node_params.size(); j++)
		{
			float dl_j = right_points_node_params[j].distance_to_backbone;
			glm::vec3 backbone_point_j(skeleton.skeletonFormat[right_points_node_params[j].point_in_backbone].point);
			
			int mimimum_index = -1;
			float minimum_node_affinity_diff = INFINITY;
			
			for (size_t k = 0; k < left_points_node_params.size(); k++)
			{
				float dl_k = left_points_node_params[k].distance_to_backbone;
				glm::vec3 backbone_point_k(skeleton.skeletonFormat[left_points_node_params[k].point_in_backbone].point);
				
				float db = glm::distance(backbone_point_j,backbone_point_k) / 
				distance_matrix[distance_matrix_start_index][left_points_node_params[k].point_in_backbone]; //difference between length in backbone hitpoints
				float dl = std::fabs(dl_j - dl_k) / (dl_j + dl_k) ; //difference between branch lengths

				int right_index = -1;
				int left_index = -1;
				//inefficient
				for (size_t t = 0; t < mesh_index_vertices.size(); t++)
				{
					if (mesh_index_vertices[t] == right_points[j])
					{
						right_index = t;
						break;
					}
				}//inefficient
				for (size_t t = 0; t < mesh_index_vertices.size(); t++)
				{
					if (mesh_index_vertices[t] == left_points[k])
					{
						left_index = t;
						break;
					}
				}
				float d_sdf = std::fabs(shape_diameter_values[right_index] - shape_diameter_values[left_index]);


				//lets try vectors instead
				glm::vec3 minimum_node_affinity_right(dl_j / (dl_j + dl_k), shape_diameter_values[right_index] , db);
				glm::vec3 minimum_node_affinity_left(dl_k / (dl_j + dl_k), shape_diameter_values[left_index] , db);
				
				float diff = glm::distance(minimum_node_affinity_right, minimum_node_affinity_left);
				if (diff < minimum_node_affinity_diff)
				{
					minimum_node_affinity_diff = diff;
					mimimum_index = k;
				}
				/*float diff = dl + db + d_sdf;
				if (diff < minimum_node_affinity_diff)
				{
					minimum_node_affinity_diff = diff; 
					mimimum_index = k;
				}*/
			}
			std::pair<unsigned int, unsigned int> point_index_pair;
			point_index_pair.first = right_points[j];
			point_index_pair.second = left_points[mimimum_index];
			backbone_pairs.push_back(point_index_pair);
			total_dif_on_backbone += minimum_node_affinity_diff;
		}
		total_dif_on_backbone /= right_points_node_params.size();
		backbone_affinity_diffs.push_back(total_dif_on_backbone);
		backbone_pairs_vec.push_back(backbone_pairs);
		backbone_right_points.push_back(right_points);
		backbone_left_points.push_back(left_points);

	}

	//check the best backbone diff
	int minimum_affinity_diff_index = -1;
	float minimum_affinity_diff = INFINITY;
	for (size_t i = 0; i < backbone_affinity_diffs.size(); i++)
	{
		if (minimum_affinity_diff > backbone_affinity_diffs[i])
		{
			minimum_affinity_diff = backbone_affinity_diffs[i];
			minimum_affinity_diff_index = i;
		}
	}
	
	// now we decided that the best backbone is minimum_affinity_diff_index
	best_backbone = candidate_backbones[minimum_affinity_diff_index];
	best_right_points = backbone_right_points[minimum_affinity_diff_index];
	best_left_points = backbone_left_points[minimum_affinity_diff_index];
	//best_backbone_pairs = backbone_pairs_vec[minimum_affinity_diff_index];
	//paint the backbone to blue
	//meshFac.mesh_skeleton_vec.skeleton_points.clear();
	glBindVertexArray(meshFac.skeleton_VAO);
	for (size_t i = 0; i < skeleton.skeletonFormat.size(); i++)
	{
		for (size_t j = 0; j < best_backbone.vertex_list.size(); j++)
		{
			if (i == best_backbone.vertex_list[j])
			{
				meshFac.mesh_skeleton_vec.skeleton_points[i * 6 + 3] = 0;
				meshFac.mesh_skeleton_vec.skeleton_points[i * 6 + 4] = 0;
				meshFac.mesh_skeleton_vec.skeleton_points[i * 6 + 5] = 255;
			}
		}
	}
	glBufferData(GL_ARRAY_BUFFER, meshFac.mesh_skeleton_vec.skeleton_points.size() * sizeof(float), &meshFac.mesh_skeleton_vec.skeleton_points[0], GL_STATIC_DRAW);


}
void skeleton_point_to_backbone(Skeleton skeleton, BackBone backbone, int index1, int& hitIndex, float& dist, std::vector<int>& indices, 
std::vector<float>& distance_matrix,std::vector<int>& predecessor_list)
{
	int N = skeleton.skeletonFormat.size();
	// 1- get the minimum distance among the backbone
	int minimum_dist_index= -1;
	float minimum_distance  = INFINITY;

	
	for (size_t i = 0; i < backbone.vertex_list.size(); i++)
	{
		if (minimum_distance > distance_matrix[backbone.vertex_list[i]])
		{
			minimum_distance = distance_matrix[backbone.vertex_list[i]];
			minimum_dist_index = backbone.vertex_list[i];
		}
	}

	skeleton_get_distance_and_vertex_list(skeleton, index1, minimum_dist_index, predecessor_list, indices, minimum_distance);
	hitIndex = minimum_dist_index;
	dist = minimum_distance;
	return;
}


void skeleton_get_distance_and_vertex_list(Skeleton&skeleton,
	int index1, int index2, std::vector<int>& predecessor_list, std::vector<int>& predecessor_index2_index1,float& geodesic_dist)
{
	int current_index = index2;
	
	geodesic_dist = 0;
	predecessor_index2_index1.clear();

	while (predecessor_list[current_index] != -1 ) //index1 
	{
		predecessor_index2_index1.push_back(current_index);

		int predecessor_index = predecessor_list[current_index];
		glm::vec3 p_current_index = skeleton.skeletonFormat[current_index].point;
		glm::vec3 p_predecessor = skeleton.skeletonFormat[predecessor_index].point;

		geodesic_dist += glm::distance(p_current_index, p_predecessor);

		current_index = predecessor_index;
	}
}



void skeleton_calculate_closest_mesh_points(Skeleton& skeleton, Mesh* m, std::vector<unsigned int >& mesh_vertex_indices)
{
	std::vector<unsigned int> end_points;
	skeleton_get_end_points(skeleton , end_points);
	for (size_t i = 0; i < end_points.size(); i++)
	{
		glm::vec3 point = skeleton.skeletonFormat[end_points[i]].point;
		int minimum_dist_index = -1;
		float minimum_distance = INFINITY;
		for (size_t j = 0; j < m->vertices.size(); j++)
		{
			float dist = glm::distance(point , m->vertices[j]);
			if (dist < minimum_distance)
			{
				minimum_distance = dist;
				minimum_dist_index = j;
			}
		}
		mesh_vertex_indices.push_back(minimum_dist_index);

	}
}
void skeleton_get_end_points(Skeleton& skeleton, std::vector<unsigned int >& end_vertex_indices)
{
	for (size_t i = 0; i < skeleton.skeletonFormat.size(); i++)
	{
		if (skeleton.skeletonFormat[i].label == END)
		{
			end_vertex_indices.push_back(i);
		}
	}
}
/*
	TODO: The N laterals should be grouped according to best backbione as a symettry plane
	1 - regroup N_laterals ( endpoints ) in according to the symmetry plane 
*/
void skeleton_get_N_Lateral_points(MeshFactory& m_factory, Skeleton& skeleton, unsigned int selected_mesh , BackBone& best_backbone ,
std::vector<std::pair<unsigned int, unsigned int>> best_backbone_point_pairs, std::vector<unsigned int>& right_mesh_indices,
std::vector<unsigned int>& left_mesh_indices)
{
	std::vector<unsigned int> left_skeleton_indices; 
	std::vector<unsigned int> right_skeleton_indices;
	std::vector<unsigned int> mesh_vertex_indices;
	

	for (size_t i = 0; i < best_backbone_point_pairs.size(); i++)
	{
		right_skeleton_indices.push_back(best_backbone_point_pairs[i].first);
		left_skeleton_indices.push_back(best_backbone_point_pairs[i].second);
	}
	//delete duplicates
	std::unordered_set<unsigned int> right_indices_unique(right_skeleton_indices.begin(), right_skeleton_indices.end());
	std::unordered_set<unsigned int> left_indices_unique(left_skeleton_indices.begin(), left_skeleton_indices.end());

	right_skeleton_indices.assign(right_indices_unique.begin(), right_indices_unique.end());
	left_skeleton_indices.assign(left_indices_unique.begin(), left_indices_unique.end());

	//now convert back to mesh points from skeleton
	skeleton_calculate_closest_mesh_points(skeleton, &m_factory.mesh_vec[selected_mesh], mesh_vertex_indices);

	// this generates every skeleton point do a  n for loop to get the needed indices for mesh
	for (size_t i = 0; i < right_skeleton_indices.size(); i++)
	{
		right_mesh_indices.push_back(right_skeleton_indices[i]);
	}
	// this generates every skeleton point do  a  n for loop to get the needed indices for mesh
	for (size_t i = 0; i < left_skeleton_indices.size(); i++)
	{
		left_mesh_indices.push_back(left_skeleton_indices[i]);
	}
}




void skeleton_generate_buffer(MeshFactory& mesh_fac)
{

	unsigned int skeleton_vao, skeleton_vbo, skeleton_ibo;
	glGenVertexArrays(1, &skeleton_vao);
	glGenBuffers(1, &skeleton_vbo);
	glGenBuffers(1, &skeleton_ibo);


	glBindVertexArray(skeleton_vao);
	glBindBuffer(GL_ARRAY_BUFFER, skeleton_vbo);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, skeleton_ibo);


	mesh_fac.skeleton_VAO = skeleton_vao;
	glBindVertexArray(0);


}
void skeleton_buffer(const MeshFactory& mesh_fac)
{
	glBindVertexArray(mesh_fac.skeleton_VAO);
	glBufferData(GL_ARRAY_BUFFER, mesh_fac.mesh_skeleton_vec.skeleton_points.size() * sizeof(float), &mesh_fac.mesh_skeleton_vec.skeleton_points[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_fac.mesh_skeleton_vec.skeleton_indices.size() * sizeof(unsigned int), &mesh_fac.mesh_skeleton_vec.skeleton_indices[0], GL_STATIC_DRAW);
	glBindVertexArray(0);

}

void skeleton_get_dijkstra_endpoints(Skeleton& skeleton, int index1, std::vector<int>& vertex_list, std::vector<float>& dijkstra_distances)
{
	int N = skeleton.skeletonFormat.size();
	//run a dijkstra from index1 
	std::vector<float> distances(N, INFINITY);
	distances[index1] = 0.0f;

	std::vector<int> discovered_vertices;
	std::vector<bool> is_vertex_discovered(N, false);
	std::vector<int> predecessors(N, -1);


	// initila nei
	discovered_vertices.push_back(index1);
	is_vertex_discovered[index1] = true;
	int no_of_discovered = 1;
	while (no_of_discovered < N)
	{
		float minimum_distance = INFINITY;
		int minimum_distance_index = -1;
		int discovered_vertex = -1;
		//check adjacent non_discovered points
		for (size_t i = 0; i < discovered_vertices.size(); i++)
		{
			int discovered_vertex_adjacency_size = skeleton.adjacencies[discovered_vertices[i]].size();
			//check adjacencies
			for (size_t j = 0; j < discovered_vertex_adjacency_size; j++)
			{
				int discovered_vertex_adjaceny = skeleton.adjacencies[discovered_vertices[i]][j];
				if (!is_vertex_discovered[discovered_vertex_adjaceny])
				{
					//check distance and compare
					float dist = glm::distance(skeleton.skeletonFormat[discovered_vertices[i]].point, skeleton.skeletonFormat[discovered_vertex_adjaceny].point);
					if (dist < minimum_distance)
					{
						minimum_distance = dist;
						minimum_distance_index = discovered_vertex_adjaceny;
						discovered_vertex = discovered_vertices[i];
					}
				}
			}
		}

		distances[minimum_distance_index] = distances[discovered_vertex] + minimum_distance;;
		predecessors[minimum_distance_index] = discovered_vertex;

		is_vertex_discovered[minimum_distance_index] = true;
		discovered_vertices.push_back(minimum_distance_index);

		//change the weights
		for (size_t i = 0; i < skeleton.adjacencies[minimum_distance_index].size(); i++)
		{
			//if (is_vertex_discovered[skeleton.adjacencies[minimum_distance_index][i]])
			{
				glm::vec3 p_index1 = skeleton.skeletonFormat[minimum_distance_index].point;
				glm::vec3 p = skeleton.skeletonFormat[skeleton.adjacencies[minimum_distance_index][i]].point;
				float edge_dist = glm::distance(p_index1, p);
				if ((edge_dist + distances[minimum_distance_index]) < distances[skeleton.adjacencies[minimum_distance_index][i]])
				{
					distances[skeleton.adjacencies[minimum_distance_index][i]] = (edge_dist + distances[minimum_distance_index]);
					predecessors[skeleton.adjacencies[minimum_distance_index][i]] = minimum_distance_index;
				}
				else
				{
					int a = 1;
				}
			}

		}
		no_of_discovered += 1;
	}

	//now the only difference with the original is calculating the end points only
	std::vector<float> dijsktra_endpoints;
	std::vector<int> predecessor_endpoints;
	for (size_t i = 0; i <distances.size() ; i++)
	{
		if (skeleton.skeletonFormat[i].label == END)
		{
			dijsktra_endpoints.push_back(distances[i]);
			predecessor_endpoints.push_back(predecessors[i]);
		}
	}
	vertex_list = predecessor_endpoints;
	dijkstra_distances = dijsktra_endpoints;
	return;
}

unsigned int skeleton_calculate_closest_mesh_point(Skeleton& skeleton, Mesh* m, unsigned int skeleton_point_index)
{
	std::vector<unsigned int> end_points;
	skeleton_get_end_points(skeleton, end_points);
	glm::vec3 skeleton_point = skeleton.skeletonFormat[skeleton_point_index].point;
	int minimum_dist_index = -1;
	float minimum_distance = INFINITY;
	for (size_t j = 0; j < m->vertices.size(); j++)
	{
		float dist = glm::distance(skeleton_point, m->vertices[j]);
		if (dist < minimum_distance)
		{
			minimum_distance = dist;
			minimum_dist_index = j;
		}
	}
	return minimum_dist_index;
	
}

std::pair<std::vector<std::pair<int, int>>,float> do_unique_pairing(std::vector<std::pair<float, int>>& distances_and_hit_points, std::vector<unsigned int>& end_point_indices
,std::vector<std::vector<float>>& end_point_dijkstras,std::vector<std::vector<int>>& end_point_vertex_list , float backbone_length )
{
	std::vector<std::pair<float, std::pair<int, int>>> compareResults;
	int N = distances_and_hit_points.size();
	for (size_t i = 0; i < N; i++)
	{
		int hitpoint_i = distances_and_hit_points[i].second;
		float dist_to_backbone_i = distances_and_hit_points[i].first;
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				continue;
			}
			int hitpoint_j = distances_and_hit_points[j].second;
			float dist_to_backbone_j = distances_and_hit_points[j].first;
			
			float dist_between_hitpoints = 0;
			float distances_to_backbone = 0;

			//search for hitpoint k in vertex list
			for (size_t t = 0; t < end_point_vertex_list[j].size(); t++)
			{
				if (hitpoint_j == end_point_vertex_list[j][t])
				{
					dist_between_hitpoints = end_point_dijkstras[j][t];
					break;
				}
			}

			//normalize dist_between 
			dist_between_hitpoints = dist_between_hitpoints / backbone_length;
			//normalize dist to backbone
			distances_to_backbone = std::abs(dist_to_backbone_i - dist_to_backbone_j) / std::max(dist_to_backbone_i, dist_to_backbone_j);
			
			//glm::vec2 res_vec(dist_between_hitpoints, distances_to_backbone);
			compareResults.push_back({ dist_between_hitpoints +  distances_to_backbone, {i,j}});
		}
	}
	std::vector<std::pair<int, int>> selectedPairs;
	float skeleton_resemblance_error = 0;
	std::vector<bool> used(N, false);

	std::sort(compareResults.begin(), compareResults.end());

	// Greedily select pairs with smallest compare() result
	for (const auto& entry : compareResults) {
		int i = entry.second.first;
		int j = entry.second.second;
		if (!used[i] && !used[j]) {
			selectedPairs.push_back({ end_point_indices[i], end_point_indices[j] });
			skeleton_resemblance_error += entry.first;
			used[i] = used[j] = true;  // Mark these objects as used
		}
	}
	
	std::pair<std::vector<std::pair<int, int>> , float> return_val;
	
	return_val.first = selectedPairs;
	return_val.second = skeleton_resemblance_error;
	
	return return_val;
}
void skeleton_generate_backbone_w_midpoint(MeshFactory& meshFac, Skeleton skeleton, unsigned int mesh_index,
	BackBone& best_backbone, std::vector<unsigned int>& right_points, std::vector<unsigned int>& left_points)
{
	Mesh* mesh = &meshFac.mesh_vec[mesh_index];
	int skeleton_mid_point_index = skeleton.mid_point_index;
	int N = skeleton.skeletonFormat.size();


	std::vector<unsigned int> end_point_indices;
	std::vector<BackBone> candidate_backbones; 
	for (size_t i = 0; i < N; i++)
	{
		if (skeleton.skeletonFormat[i].label == END)
		{
			end_point_indices.push_back(i);
		}
	}
	
	int N_end = end_point_indices.size();
	std::vector<std::vector<float>> end_point_dijkstras(N_end);
	std::vector<std::vector<int>> end_point_vertex_list(N_end);
	for (size_t i = 0; i < N_end; i++)
	{
		skeleton_calculate_dijkstra(skeleton, end_point_indices[i], end_point_vertex_list[i], end_point_dijkstras[i]);
	}

	std::vector<std::vector<int>> predecessor_list_for_end_points(N_end);
	std::vector<std::vector<float>> distance_matrix(N_end);
	for (size_t i = 0; i < N_end; i++)
	{
		int index1 = end_point_indices[i];
		std::vector<int> point_list;
		std::vector<float> distances_from_index1;

		skeleton_calculate_dijkstra(skeleton, index1, point_list, distances_from_index1);

		for (size_t j = i; j < N_end; j++)
		{
			//assume i and j are the end points of backbone
			//generate the path
			int index2 = end_point_indices[j];
			float backbone_length;
			std::vector<int> vertex_list_index1_index2;
			if (i != j)
			{
				skeleton_get_distance_and_vertex_list(skeleton, index1, index2, point_list,
					vertex_list_index1_index2, backbone_length);

				BackBone backbone;
				backbone.start_index = index1;
				backbone.end_index = index2;
				backbone.vertex_list = vertex_list_index1_index2;
				/*bool is_mid_point_present = false; // mid point deactive !!!!!!!!!!!!
				// check if mid point index included
				for (size_t k = 0; k < backbone.vertex_list.size(); k++)
				{
					if (backbone.vertex_list[k] == skeleton_mid_point_index)
					{
						is_mid_point_present = true; 
						break;
					}
				}*/
				if (true)
				{
					candidate_backbones.push_back(backbone);
				}

			}
			predecessor_list_for_end_points[i] = point_list;

			distance_matrix[i] = distances_from_index1;

		}
	}
	//end results are obtained 
	std::vector<std::vector<std::pair< int, int>>> skeleton_resemblances(candidate_backbones.size());
	std::vector<float> skeleton_resemblances_score(candidate_backbones.size() , 0);
	for (size_t i = 0; i < candidate_backbones.size(); i++)
	{

		// 1 - get 	backbone total length
		float backbone_i_length = skeleton_get_backbone_length(mesh, &candidate_backbones[i]);
		std::vector<std::pair<float, int> > distances_and_hit_points; // first hitpoint distance from point to backbone
		// second is the skeleton index
		// 2 - divide end points of backbone
	/*	skeleton_divide_end_points_left_right(mesh, skeleton, end_point_indices, candidate_backbones[i],
		end_point_vertex_list,right_points , left_points);*/
		for (size_t j = 0; j < N_end; j++)
		{
			float dist; 
			int hit_index;
			std::vector<int> indices;
			std::vector<float> distance_matrix ;
			std::vector<int> predecessor_list;
			skeleton_point_to_backbone(skeleton, candidate_backbones[i], end_point_indices[j],
			hit_index,dist,indices,end_point_dijkstras[j], end_point_vertex_list[j]);

			distances_and_hit_points.push_back(std::pair<float,int>(dist,hit_index));
		}
		std::vector<std::pair< int, int>> skeleton_i_resemblances;
		// do a N^2 comaprison of each distance and hit point  


		 
		std::pair<std::vector<std::pair<int,int>>, float > skeleton_pair_vals  = do_unique_pairing(distances_and_hit_points , end_point_indices , end_point_dijkstras , end_point_vertex_list , backbone_i_length );
		skeleton_resemblances[i] = skeleton_pair_vals.first;
		skeleton_resemblances_score[i] = skeleton_pair_vals.second;
		/*for (size_t j = 0; j < N_end; j++)
		{
			float min_dist = INFINITY;
			int min_index = -1;
 
			int hitpoint_j = distances_and_hit_points[j].second;
			float dist_to_backbone_j = distances_and_hit_points[j].first;
			// do a skeleton dijkstra from hitpoint

			for (size_t k = 0; k < N_end; k++)
			{
				if (j == k)
				{
					continue; 
				}
				// 1 - compare distance to backbone
				float dist_to_backbone_k = distances_and_hit_points[k].second;
				int  hit_point_k = distances_and_hit_points[k].second;
				float dist_between_hitpoints = INFINITY;
				float distances_to_backbone;
				
				//search for hitpoint k in vertex list
				for (size_t t = 0; t < end_point_vertex_list[j].size(); t++)
				{
					if (hit_point_k == end_point_vertex_list[j][t])
					{
						dist_between_hitpoints = end_point_dijkstras[j][t];
						break; 
					}
				}

				//normalize dist_between 
				dist_between_hitpoints = dist_between_hitpoints / backbone_i_length;
				//normalize dist to backbone
				distances_to_backbone = std::abs(dist_to_backbone_j - dist_to_backbone_k) / std::max(dist_to_backbone_j, dist_to_backbone_k);
				
				// voting scheme
				float total_dist = 0.5f * dist_between_hitpoints + distances_to_backbone * 0.5;
				if (total_dist < min_dist)
				{
					min_dist = total_dist;
					min_index = k;
				}
			}
			skeleton_i_resemblances.push_back(std::pair<int,int>(end_point_indices[j], end_point_indices[min_index]));
			//skeleton_resemblances_score[i] += min_dist;
		}
		skeleton_resemblances[i] = skeleton_i_resemblances;*/


	}

	//select with the lowest score
	int min_index = -1;
	float min_val = INFINITY;
	for (size_t i = 0; i < skeleton_resemblances_score.size(); i++)
	{
		if (skeleton_resemblances_score[i] < min_val)
		{
			min_val = skeleton_resemblances_score[i];
			min_index = i; 
		}
	}

	for (size_t i = 0; i < N_end/2; i++)
	{
		right_points.push_back(skeleton_resemblances[min_index][i].first);
		left_points.push_back(skeleton_resemblances[min_index][i].second);
	}
	best_backbone = candidate_backbones[min_index];
	glBindVertexArray(meshFac.skeleton_VAO);
	for (size_t i = 0; i < skeleton.skeletonFormat.size(); i++)
	{
		for (size_t j = 0; j < best_backbone.vertex_list.size(); j++)
		{
			if (i == best_backbone.vertex_list[j])
			{
				meshFac.mesh_skeleton_vec.skeleton_points[i * 6 + 3] = 0;
				meshFac.mesh_skeleton_vec.skeleton_points[i * 6 + 4] = 0;
				meshFac.mesh_skeleton_vec.skeleton_points[i * 6 + 5] = 255;
			}
		}
	}
	glBufferData(GL_ARRAY_BUFFER, meshFac.mesh_skeleton_vec.skeleton_points.size() * sizeof(float), &meshFac.mesh_skeleton_vec.skeleton_points[0], GL_STATIC_DRAW);
}

float skeleton_get_backbone_length(Mesh* m, BackBone* backBone)
{
	float total_length = 0;
	for (size_t i = 0; i < backBone->vertex_list.size()-1; i++)
	{
		int mesh_index_1 = backBone->vertex_list[i];
		int mesh_index_2 = backBone->vertex_list[i+1];
		glm::vec3 p1 = m->vertices[mesh_index_1];
		glm::vec3 p2 = m->vertices[mesh_index_2];
		total_length = total_length + glm::distance(p1, p2);
	}
	return total_length;
}

void skeleton_divide_end_points_left_right(Mesh* mesh, Skeleton& skeleton, std::vector<int> end_point_indices, BackBone candidate_backbone,
 std::vector<std::vector<int>>& end_point_vertex_list,std::vector<std::vector<float>>& end_point_distances,
std::vector<int>& left_points , std::vector<int>& right_points)
{
	int N_end = end_point_indices.size();
	for (size_t i = 0; i < N_end; i++)
	{
		int hit_index = -1;
		float dist;
		std::vector<int> indices;
		skeleton_point_to_backbone(skeleton, candidate_backbone, end_point_indices[i],
		hit_index, dist, indices, end_point_distances[i], end_point_vertex_list[i]);



	}
}

//used for ladnmarked meshes 
void skeleton_left_right_test_for_endpoint(Mesh* m , Skeleton* skeleton,
std::vector<int>& right, std::vector<int>& left)
{
	
}
