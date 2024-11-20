#include "../Include/SCBDataset.h"
#include <string>
#include "../Include/HeatKernelSignature.h"
std::vector<TrilateralMesh> SCB_mesh; 
typedef enum
{
	SCAPE,
	TOSCA,
	WATERTIGHT
}DATASET;

#define SCAPE_DATA_SIZE 20
#define TOSCA_CAT_SIZE 11
static DATASET dataset; 
void SCB_read_SCAPE()
{
	dataset = SCAPE;
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
void SCB_read_TOSCA()
{
	dataset = TOSCA;
	for (size_t i = 0; i < TOSCA_CAT_SIZE; i++)
	{
		std::string mesh_name = "/cat";
		std::string i_name = std::to_string(i);
	
		mesh_name = TOSCA_PATH + mesh_name + i_name + ".off";
		TrilateralMesh m((char*)mesh_name.c_str());
		SCB_mesh.push_back(m);

	}
}
void SCB_select_mesh(TrilateralMesh& m , int meshNo, Skeleton& skeleton) 
{
	std::string hks_path; 
	std::string skeleton_path;
	if (dataset == SCAPE)
	{
		hks_path = "../../Trilateral/Mesh/SCB/Data/SCAPE/HKS/";
		skeleton_path = "../../Trilateral/Mesh/SCB/Data/SCAPE/skeleton/";
	}
	else if (dataset == TOSCA)
	{
		hks_path = "../../Trilateral/Mesh/SCB/Data/TOSCA/HKS/";
		skeleton_path = "../../Trilateral/Mesh/SCB/Data/TOSCA/skeleton/";
	}
	skeleton_path = skeleton_path + m.file_name.substr(0,m.file_name.size() - 3) + "swc";
	m = SCB_mesh[meshNo];
	HKS_read_kernel_signature(&m , hks_path );
	skeleton = skeleton_read_swc_file(&m, skeleton_path);
}