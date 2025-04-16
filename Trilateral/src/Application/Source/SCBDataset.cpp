#include "../Include/SCBDataset.h"
#include <string>
#include "../Include/HeatKernelSignature.h"
std::vector<TrilateralMesh> SCB_mesh; 

#define SCAPE_DATA_SIZE 74
#define TOSCA_CAT_SIZE 13
#define TOSCA_CENTAUR_SIZE 7
#define TOSCA_DOG_SIZE 8
#define PRINCETON_DATASET_SIZE 400
static DATASET dataset; 

void SCB_select_dataset(DATASET cur_dataset)
{
	dataset = cur_dataset; 
}
void SCB_read_SCAPE()
{
	SCB_mesh.clear();
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
	SCB_mesh.clear();
	dataset = TOSCA;
	for (size_t i = 0; i < TOSCA_CAT_SIZE; i++)
	{
		std::string mesh_name = "/cat";
		std::string i_name = std::to_string(i);
	
		mesh_name = TOSCA_PATH + mesh_name + i_name + ".off";
		TrilateralMesh m((char*)mesh_name.c_str());
		SCB_mesh.push_back(m);

	}
	for (size_t i = 0; i < TOSCA_CENTAUR_SIZE; i++)
	{
		std::string mesh_name = "/centaur";
		std::string i_name = std::to_string(i);

		mesh_name = TOSCA_PATH + mesh_name + i_name + ".off";
		TrilateralMesh m((char*)mesh_name.c_str());
		SCB_mesh.push_back(m);
	}
	for (size_t i = 0; i < TOSCA_DOG_SIZE; i++)
	{
		std::string mesh_name = "/dog";
		std::string i_name = std::to_string(i);

		mesh_name = TOSCA_PATH + mesh_name + i_name + ".off";
		TrilateralMesh m((char*)mesh_name.c_str());
		SCB_mesh.push_back(m);
	}
}
void SCB_read_Princeton()
{
	SCB_mesh.clear();
	dataset = PRINCETON;
	for (size_t i = 1; i < PRINCETON_DATASET_SIZE; i++)
	{
		if (i == 262)
		{
			i = 281;
		}
		std::string mesh_name;
		std::string i_name = std::to_string(i);
		std::string tab = "/";
		mesh_name = PRINCETON_PATH  +  tab + i_name + ".off";
		std::cout << mesh_name << std::endl; 
		TrilateralMesh m((char*)mesh_name.c_str());
		SCB_mesh.push_back(m);

	}
}
void SCB_read_KIDS()
{
	dataset = KIDS; 
	SCB_mesh.clear();
	// total of 15 + 15 meshes 
	for (size_t i = 0; i < 10; i++)
	{
		std::string path("../../Trilateral/Mesh/off/");
		std::string isometry_batch_no("000" + std::to_string((i / 15) + 1));
		std::string isometry_no(std::to_string(i % 15 + 1));
		// read the meshes.
		path = path + isometry_batch_no + ".isometry." + isometry_no + ".off";
		TrilateralMesh m((char*)path.c_str());
		//read the symmetry format
		read_symmetry_format((char*)"../../Trilateral/Mesh/off/sym.txt", &m);
		std::string hks_path = "../../Trilateral/Mesh/off/HKS/";
		HKS_read_kernel_signature(&m, hks_path);
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
	else if (dataset == KIDS)
	{
		hks_path = "../../Trilateral/Mesh/off/HKS/";
		skeleton_path = "../../Trilateral/Mesh/off/KIDS_skeleton/";
		read_symmetry_format((char*)"../../Trilateral/TrilateralMesh/off/sym.txt", &m);

	}
	else if (dataset == PRINCETON)
	{
		hks_path = "../../Trilateral/Mesh/SCB/Data/Princeton/HKS/";
	}
	m = SCB_mesh[meshNo];
	skeleton_path = skeleton_path + m.file_name.substr(0,m.file_name.size() - 3) + "swc";
	HKS_read_kernel_signature(&m , hks_path );
	skeleton = skeleton_read_swc_file(&m, skeleton_path);
}
