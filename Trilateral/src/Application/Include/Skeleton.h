#include <vector>
#include <string>
#include <map>
#include <glm/ext/vector_float3.hpp>
#include "../Include/Mesh.h"
#include "../Include/MeshFactory.h"

std::vector<float> generate_bounding_box(std::string file_name);

std::map<std::string, glm::vec3 > generate_skeleton_keypoints(std::string file_name);

void match_skeleton_keypoints( MeshFactory& meshFactory, Mesh* m, std::vector<float>& skeleton_bounding_box, std::map<std::string, glm::vec3>& keypoints);

