#pragma once 
#include <vector>
#include <string>
#include <map>
#include <glm/ext/vector_float3.hpp>
#include "../Include/TrilateralMesh.h"
#include "../Include/DvorakEstimatingApprox.h"
#include "../Include/MeshFactory.h"

std::vector<float> generate_bounding_box(std::string file_name);
std::map<std::string, glm::vec3 > generate_skeleton_keypoints(std::string file_name);
void match_skeleton_keypoints( MeshFactory& meshFactory, TrilateralMesh* m, std::vector<float>& skeleton_bounding_box, std::map<std::string, glm::vec3>& keypoints);
void match_skeleton_lines(MeshFactory& meshFactory, TrilateralMesh* m, std::vector<float>& skeleton_bounding_box, std::vector<float> skeleton_lines);
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
	Color color;
}SkeletonFormat;
typedef struct
{
	unsigned int index;
}SkeletonEndPoint;
typedef struct
{
	std::vector<SkeletonFormat> skeletonFormat;
	std::vector<SkeletonEndPoint> endPoints;
	std::vector<std::vector<unsigned int>> adjacencies;
	glm::vec3 skeleton_mid_point;
	int mid_point_index; // index of skeleton point's who is closest to the mid point
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

struct SkeletonTreeNode
{
	int skeleton_index;
	int depth;
	SkeletonTreeNode* parent;
	std::vector<SkeletonTreeNode*> child;
};
typedef struct
{
	SkeletonTreeNode* head; // mid point
}SkeletonTree;
//buffer
void skeleton_generate_buffer( MeshFactory& mesh_fac);
void skeleton_buffer(const MeshFactory& mesh_fac);

Skeleton skeleton_read_swc_file(TrilateralMesh* m , std::string file_name);

void skeleton_calculate_dijkstra(Skeleton skeleton, int index1,std::vector<int>& vertex_list, std::vector<float>& dijkstra_distances);
void skeleton_generate_backbone(TrilateralMesh* m, Skeleton skeleton,
BackBone& best_backbone, std::vector<unsigned int>& right_points  , std::vector<unsigned int>& left_points);

void skeleton_generate_backbone_w_midpoint(MeshFactory& meshFac, Skeleton skeleton, unsigned int mesh_index,
	BackBone& best_backbone, std::vector<unsigned int>& right_points, std::vector<unsigned int>& left_points);

void skeleton_point_to_backbone(Skeleton skeleton, BackBone backbone, int index1, int& hitIndex, float& dist, std::vector<int>& indices,
	std::vector<float>& distance_matrix, std::vector<int>& predecessor_list_for_end_points);
void skeleton_get_distance_and_vertex_list(Skeleton& skeleton,
	int index1, int index2, std::vector<int>& predecessor_list, std::vector<int>& predecessor_index2_index1, float& geodesic_dist);
void skeleton_calculate_closest_mesh_points(Skeleton& skeleton, TrilateralMesh* m, std::vector<unsigned int >& mesh_vertex_indices);
unsigned int skeleton_calculate_closest_mesh_point(Skeleton& skeleton, TrilateralMesh* m, unsigned int skeleton_point_index);
void skeleton_get_end_points(Skeleton& skeleton, std::vector<unsigned int >& mesh_vertex_indices);
void skeleton_get_end_points_update_mesh(TrilateralMesh* m, Skeleton& skeleton, std::vector<unsigned int >& end_vertex_indices);
void skeleton_get_N_Lateral_points(MeshFactory& m_factory, Skeleton& skeleton, unsigned int selected_mesh, BackBone& best_backbone,
	std::vector<std::pair<unsigned int, unsigned int>> best_backbone_point_pairs, std::vector<unsigned int>& right_mesh_indices,
	std::vector<unsigned int>& left_mesh_indices);
void skeleton_get_dijkstra_endpoints(Skeleton& skeleton, int index1, std::vector<int>& vertex_list, std::vector<float>& dijkstra_distances);

float skeleton_get_backbone_length(TrilateralMesh* m,BackBone* backBone);

SkeletonTree skeleton_generate_skeleton_tree(TrilateralMesh* m, Skeleton& skeleton);
SkeletonTreeNode skeleton_get_skeleton_node(Skeleton& skeleton, SkeletonTree& skelTree, int skeletonIndex);
void skeleton_get_closest_skeleton_endpoints(TrilateralMesh* m, Skeleton& skeleton, std::vector<unsigned int>& mesh_points, 
std::vector<unsigned int >& skeleton_end_points);

void skeleton_left_right_test_for_endpoint(std::vector<int>& right , std::vector<int>& left );

void skeleton_generate_backbone_with_dvorak_pairs(TrilateralMesh* m, Skeleton& skeleton,BackBone& b, std::vector<DvorakPairs>& dvorakPairs  );
//TODO
void skeleton_generate_backbone_with_dominant_sym(MeshFactory& meshFactory , int selected_mesh , Skeleton& skeleton);

std::vector<float> skeleton_distance_to_midpoint(TrilateralMesh* m, Skeleton& skeleton, std::vector<unsigned int> indices );