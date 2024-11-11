#pragma once 
#include <vector>
#include <string>
#include <map>
#include <glm/ext/vector_float3.hpp>
#include "../Include/TrilateralMesh.h"
//for cohen or 
enum POINT_LABEL {
	UNDEFINED = 0,
	SOMA = 1,
	FORK = 2,
	END = 6,
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