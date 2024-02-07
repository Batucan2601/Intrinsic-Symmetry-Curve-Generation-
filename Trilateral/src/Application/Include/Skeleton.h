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

void skeleton_read_swc_file(MeshFactory& meshFactory , std::string file_name);
