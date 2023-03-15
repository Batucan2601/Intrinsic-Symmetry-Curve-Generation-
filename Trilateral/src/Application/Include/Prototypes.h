#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>

#pragma region prototypes
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
unsigned int set_shader(std::string path);
void poll_keys();
void imgui_new_frame();
void imgui_setup();
void imgui_render();
void imgui_close();
void imgui_histogram(std::vector<float>& histogram, int partition_size);
void imgui_singled_rocket_simulation_popup(bool is_hit, glm::vec3  lidar1_intersection, glm::vec3 lidar2_intersection, glm::vec3 lidar3_intersection);
void imgui_single_rocket_simulation_popup(int& lidar1_smallest, int& lidar2_smallest, int& lidar3_smallest, int& lidar1_biggest, int& lidar2_biggest, int& lidar3_biggest, float& total_area1, float& total_area2, float& total_area3, float& distance1, float& distance2, float& distance3, glm::vec3 vehicle_hit_highest_point_x, glm::vec3 vehicle_hit_lowest_point_x, float time_of_hit_lidar1, float  time_of_hit_lidar2, float  time_of_hit_lidar3);
void imgui_multiple_Rocket_simulation_popup();
void imgui_end_of_simulation_popup();
void init_hashtable();
void imgui_input_window();
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
glm::vec3 generate_random_point();
void generate_direction(glm::vec3& random_point);
#pragma region pix2hitarea.h
void rocket_hit_area(int  pix_min_1, int  pix_max_1, float  rocket_dist_1, int pix_min_2, int pix_max_2, float rocket_dist_2, int pix_min_3, int pix_max_3, float rocket_dist_3, glm::vec3* vehicle_hit_highest_point_x, glm::vec3* vehicle_hit_lowest_point_x);
void pixToPositionXYZ(int pix_min, int pix_max, int rocket_dist, double wall_angle, glm::vec3 lidar_position, glm::vec3* return_start_point, glm::vec3* return_end_point);
void calculate_hitpoint_on_trajectory(glm::vec3 p1, glm::vec3 p2, glm::vec3* return_hitpoint);
void vehicleHitRegionCalculation(glm::vec3 lidar_1_low, glm::vec3 lidar_1_high, glm::vec3 lidar_2_low, glm::vec3 lidar_2_high, glm::vec4 lidar_3_low, glm::vec3 lidar_3_high, glm::vec3* vehicle_region_x_high, glm::vec3* vehicle_region_x_low);
float time_to_hit_single_rocket(float z);
#pragma endregion pix2hitarea.h 
#pragma region simulation 
float time_passed_in_intersection(const glm::vec3& hit_point);
#pragma endregion 
//std::vector<glm::vec4> rocket_hit_area(int  pix_min_1, int  pix_max_1, float  rocket_dist_1, int pix_min_2, int pix_max_2, float rocket_dist_2, int pix_min_3, int pix_max_3, float rocket_dist_3);

void match_points_from2_mesh(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2,  int division_no); // sample size must be bigger than 3 
std::vector<float> match_points_from2_mesh_mock(MeshFactory& mesh_fac, int mesh_index1, int mesh_index2, int division_no); // sample size must be bigger than 3 
void brute_force_symmetry_extraction(MeshFactory& mesh_fac, int& selected_index);




#pragma endregion