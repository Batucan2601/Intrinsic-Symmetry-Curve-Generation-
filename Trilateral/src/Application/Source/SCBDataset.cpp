#include "../Include/SCBDataset.h"
#include <string>
#include "../Include/HeatKernelSignature.h"
std::vector<TrilateralMesh> SCB_mesh; 

#define SCAPE_DATA_SIZE 20
void SCB_read_SCAPE()
{
	for (size_t i = 0; i < SCAPE_DATA_SIZE; i++)
	{
		std::string mesh_name = "/mesh";
		std::string i_name = std::to_string(i);
		std::string zero;
		if (i < 10)
		{
			zero = "00";
		}
		else
		{
			zero = "0";
		}
		mesh_name = SCAPE_PATH +  mesh_name + zero + i_name + ".off";
		TrilateralMesh m((char*)mesh_name.c_str());
		SCB_mesh.push_back(m);
		
	}
}
void SCB_select_mesh(TrilateralMesh& m , int meshNo, Skeleton& skeleton) 
{
	std::string hks_path = "../../Trilateral/Mesh/SCB/Data/SCAPE/HKS/";
	std::string skeleton_path = "../../Trilateral/Mesh/SCB/Data/SCAPE/skeleton/";
	skeleton_path = skeleton_path + m.file_name.substr(0,m.file_name.size() - 3) + "swc";
	m = SCB_mesh[meshNo];
	HKS_read_kernel_signature(&m , hks_path );
	skeleton = skeleton_read_swc_file(&m, skeleton_path);
}