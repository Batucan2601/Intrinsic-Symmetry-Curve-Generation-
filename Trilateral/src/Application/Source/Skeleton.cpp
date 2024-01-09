#include "../Include/Skeleton.h"
#include "../Include/Mesh.h"
#include <iostream>
#include <sstream>
#include <fstream>
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

void match_skeleton_keypoints(Mesh* m , std::vector<float>& skeleton_bounding_box, std::map<std::string, glm::vec3>& keypoints)
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
	glm::vec3 skeleton_to_world_space_axis((bb_mesh[3] - bb_mesh[0]) / (skeleton_bounding_box[3] - skeleton_bounding_box[0]),
		(bb_mesh[4] - bb_mesh[1]) / (skeleton_bounding_box[4] - skeleton_bounding_box[1]),
		(bb_mesh[5] - bb_mesh[2]) / (skeleton_bounding_box[5] - skeleton_bounding_box[2]));



}
