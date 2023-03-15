#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include <vector>
void calculate_hitpoint_on_trajectory(glm::vec3 p1, glm::vec3 p2, glm::vec3* return_hitpoint)
{
    glm::vec3 hitpoint = glm::vec3();

    double vehicle_position_z = 0;

    double x_dif_p1_p2 = p2.x - p1.x;

    double y_dif_p1_p2 = p2.y - p1.y;

    double z_dif_p1_p2 = p2.z - p1.z;

    double z_dif_p2_vehicle = vehicle_position_z - p2.z;

    double move_ratio = z_dif_p2_vehicle / z_dif_p1_p2;


    hitpoint.x = p2.x + move_ratio * x_dif_p1_p2;


    hitpoint.y = p2.y + move_ratio * y_dif_p1_p2;


    hitpoint.z = p2.z + move_ratio * z_dif_p1_p2;

    *return_hitpoint = hitpoint;

}
float time_to_hit_single_rocket(float z)
{
    float z_distance_until_hit = glm::abs(z);
    float velocity_x_dir = glm::tan(glm::radians(global_variables["roket gelis acisi x"]));
    float velocity_y_dir = glm::tan(glm::radians(global_variables["roket gelis acisi y"]));
    float velocity_z_dir = 1.0f;
    float z_speed = global_variables["roket hizi"] * velocity_z_dir / (glm::sqrt(glm::pow(velocity_x_dir, 2) + glm::pow(velocity_y_dir, 2) + glm::pow(velocity_z_dir, 2)));
    float time_until_hit = z_distance_until_hit / z_speed; // m/(m/s)  in seconds
    return time_until_hit;
}
float time_passed_in_intersection(const glm::vec3 & hit_point )
{
    //get the total time to hit
   float total_time =  time_to_hit_single_rocket(global_variables["roket pozisyonu z"]);
   // the hit time =  (total_time - hit time calculated in perde )
   float hit_titme_after_intersection = time_to_hit_single_rocket(hit_point[2]);

   return total_time - hit_titme_after_intersection;
}