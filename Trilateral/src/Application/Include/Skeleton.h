#pragma once 
#include <vector>
#include <string>
#include <map>
#include <glm/ext/vector_float3.hpp>
#include "../Include/Mesh.h"
#include "../Include/MeshFactory.h"

std::vector<float> generate_bounding_box(std::string file_name);
std::map<std::string, glm::vec3 > generate_skeleton_keypoints(std::string file_name);
void match_skeleton_keypoints( MeshFactory& meshFactory, Mesh* m, std::vector<float>& skeleton_bounding_box, std::map<std::string, glm::vec3>& keypoints);
void match_skeleton_lines(MeshFactory& meshFactory, Mesh* m, std::vector<float>& skeleton_bounding_box, std::vector<float> skeleton_lines);
std::vector<float> generate_skeleton_lines(std::string file_name);

//for cohen or 
enum POINT_LABEL{ 
	UNDEFINED = 0,
	SOMA = 1,
	FORK = 2, 
	END =6,
};
typedef struct
{
	int parent;
	POINT_LABEL label;
	glm::vec3 point;
}SkeletonFormat;
typedef struct
{
	std::vector<SkeletonFormat> skeletonFormat;
	std::vector<std::vector<unsigned int>> adjacencies;
}Skeleton;
typedef struct
{
	std::vector<int> vertex_list; 
	unsigned int start_index;
	unsigned int end_index;
}BackBone;
typedef struct
{
	float distance_to_backbone;
	unsigned int point_in_backbone;
}NodeAffinityParams;


Skeleton skeleton_read_swc_file(MeshFactory& meshFactory , std::string file_name);

void skeleton_calculate_dijkstra(Skeleton skeleton, int index1,std::vector<int>& predecessor_list, std::vector<float>& dijkstra_distances);
void skeleton_generate_backbone(MeshFactory& meshFac, Skeleton skeleton);
void skeleton_point_to_backbone(Skeleton skeleton, BackBone backbone, int index1, int& hitIndex, float& dist, std::vector<int>& indices,
	std::vector<float>& distance_matrix, std::vector<int>& predecessor_list_for_end_points);
void skeleton_get_distance_and_vertex_list(Skeleton& skeleton,
	int index1, int index2, std::vector<int>& predecessor_list, std::vector<int>& predecessor_index2_index1, float& geodesic_dist);
