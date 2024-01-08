#include "../Include/Skeleton.h"
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
	
	indata.open("../../Trilateral/Mesh/off/KIDS_skeleton" + file_name);
	if (!indata)
	{
		return bounding_box;
	}
	std::string line;
	while (std::getline(indata, line))
	{
		if (line_count > 16)
		{
			line_count += 1;
		}
		else
		{
			//read the values
			std::stringstream ss(line);
			bounding_box.push_back(0); // empty val
			ss >> bounding_box[bounding_box_index++];
		}

	}
	indata.close();

	return bounding_box;
}