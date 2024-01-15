#include "../Include/Skeleton.h"
#include "../Include/Mesh.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <src/Application/Include/MeshFactory.h>
#include <GL/glew.h>

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
		if (line_count > 16)
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
	it = keypoints.begin();
	while (it != keypoints.end())
	{
		glm::vec3 temp = it->second;

		float temp_f = temp.y;
		temp.y = temp.z;
		temp.z = temp_f;
		it->second = temp;
		it++;
	}

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
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LEye"].x);
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

	// 3 -between Shoulders 
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);

	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["RShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);


	// 4 - Full arm 
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].x);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].y);
	meshFactory.mesh_skeleton_vec.push_back(keypoints["LShoulder"].z);

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
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

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
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

	meshFactory.mesh_skeleton_vec.push_back(0.0f);
	meshFactory.mesh_skeleton_vec.push_back(0.0f);
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
	meshFactory.mesh_skeleton_vec.push_back(255.0f);


	unsigned int skeleton_vao, skeleton_vbo;
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

	glBindVertexArray(0);

}
