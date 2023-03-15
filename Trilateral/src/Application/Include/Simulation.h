#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include <vector>
void pixToPositionXYZ(int pix_min, int pix_max, double rocket_dist, double wall_angle_degree, glm::vec3 lidar_position, glm::vec3* return_start_point, glm::vec3* return_end_point)
{
    double rocket_diameter = global_variables["roket capi"];
    double wall_angle = wall_angle_degree / 180 * glm::pi<double>();
    double fov_rad = global_variables["FOV"] * glm::pi<double>() / 180;
    double ifov_rad = fov_rad / global_variables["sensor sayisi"];


    int active_pixel_count = (pix_max - pix_min) + 1;
    double pixel_width = ifov_rad * rocket_dist;

    double pixel_border_to_rpg_angle = glm::asin((rocket_diameter / 2) / rocket_dist);


    glm::vec3 rocket_center_point_lowest = glm::vec3();
    rocket_center_point_lowest = lidar_position;

    glm::vec3 rocket_center_point_highest = glm::vec3();
    rocket_center_point_highest = lidar_position;

    if (rocket_diameter < (active_pixel_count - 1) * pixel_width)
    {
        double min_pixel_right_angle = -1 * fov_rad / 2 + (pix_min)*ifov_rad;
        rocket_center_point_highest.x = lidar_position.x + rocket_dist * glm::sin(min_pixel_right_angle + pixel_border_to_rpg_angle);
        rocket_center_point_highest.y = lidar_position.y - rocket_dist * glm::cos(min_pixel_right_angle + pixel_border_to_rpg_angle) * glm::cos(wall_angle);
        rocket_center_point_highest.z = lidar_position.z - rocket_dist * glm::cos(min_pixel_right_angle + pixel_border_to_rpg_angle) * glm::sin(wall_angle);

        double max_pixel_left_angle = -1 * fov_rad / 2 + (pix_max - 1) * ifov_rad;
        rocket_center_point_lowest.x = lidar_position.x + rocket_dist * glm::sin(max_pixel_left_angle - pixel_border_to_rpg_angle);
        rocket_center_point_lowest.y = lidar_position.y - rocket_dist * glm::cos(max_pixel_left_angle - pixel_border_to_rpg_angle) * glm::cos(wall_angle);
        rocket_center_point_lowest.z = lidar_position.z - rocket_dist * glm::cos(max_pixel_left_angle - pixel_border_to_rpg_angle) * glm::sin(wall_angle);
    }
    else if (rocket_diameter > (active_pixel_count - 1) * pixel_width)
    {
        double min_pixel_left_angle = -1 * fov_rad / 2 + (pix_min - 1) * ifov_rad;
        rocket_center_point_lowest.x = lidar_position.x + rocket_dist * glm::sin(min_pixel_left_angle + pixel_border_to_rpg_angle);
        rocket_center_point_lowest.y = lidar_position.y - rocket_dist * glm::cos(min_pixel_left_angle + pixel_border_to_rpg_angle) * glm::cos(wall_angle);
        rocket_center_point_lowest.z = lidar_position.z - rocket_dist * glm::cos(min_pixel_left_angle + pixel_border_to_rpg_angle) * glm::sin(wall_angle);

        double max_pixel_right_angle = -1 * fov_rad / 2 + (pix_max)*ifov_rad;
        rocket_center_point_highest.x = lidar_position.x + rocket_dist * glm::sin(max_pixel_right_angle - pixel_border_to_rpg_angle);
        rocket_center_point_highest.y = lidar_position.y - rocket_dist * glm::cos(max_pixel_right_angle - pixel_border_to_rpg_angle) * glm::cos(wall_angle);
        rocket_center_point_highest.z = lidar_position.z - rocket_dist * glm::cos(max_pixel_right_angle - pixel_border_to_rpg_angle) * glm::sin(wall_angle);
    }

    /*std::cout << rocket_center_point_lowest.x << std::endl;
    std::cout << rocket_center_point_lowest.y << std::endl;
    std::cout << rocket_center_point_lowest.z << std::endl;

    std::cout << rocket_center_point_highest.x << std::endl;
    std::cout << rocket_center_point_highest.y << std::endl;
    std::cout << rocket_center_point_highest.z << std::endl;*/

    *return_start_point = rocket_center_point_lowest;
    *return_end_point = rocket_center_point_highest;
}

void vehicleHitRegionCalculation(glm::vec3 lidar_1_low, glm::vec3 lidar_1_high, glm::vec3 lidar_2_low, glm::vec3 lidar_2_high, glm::vec3 lidar_3_low, glm::vec3 lidar_3_high, glm::vec3* vehicle_region_x_high, glm::vec3* vehicle_region_x_low)
{
    glm::vec3 hitpoint_x_high = glm::vec3();//initialized as pos inf so that any point is better 
    hitpoint_x_high.x = std::numeric_limits<double>::infinity();
    hitpoint_x_high.y = std::numeric_limits<double>::infinity();
    hitpoint_x_high.z = std::numeric_limits<double>::infinity();

    glm::vec3 poss_hitpoint_high = glm::vec3();
    calculate_hitpoint_on_trajectory(lidar_1_low, lidar_2_high, &poss_hitpoint_high);
    if (poss_hitpoint_high.x < hitpoint_x_high.x)
    {
        hitpoint_x_high = poss_hitpoint_high;
    }
    calculate_hitpoint_on_trajectory(lidar_1_low, lidar_3_high, &poss_hitpoint_high);
    if (poss_hitpoint_high.x < hitpoint_x_high.x)
    {
        hitpoint_x_high = poss_hitpoint_high;
    }    
    calculate_hitpoint_on_trajectory(lidar_2_low, lidar_3_high, &poss_hitpoint_high);
    if (poss_hitpoint_high.x < hitpoint_x_high.x)
    {
        hitpoint_x_high = poss_hitpoint_high;
    }

    glm::vec3 hitpoint_x_low = glm::vec3(); //initialized as neg inf so that any point is better 
    hitpoint_x_low.x = -std::numeric_limits<double>::infinity();
    hitpoint_x_low.y = -std::numeric_limits<double>::infinity();
    hitpoint_x_low.z = -std::numeric_limits<double>::infinity();

    glm::vec3 poss_hitpoint_low = glm::vec3();
    calculate_hitpoint_on_trajectory(lidar_1_high, lidar_2_low, &poss_hitpoint_low);
    if (poss_hitpoint_low.x > hitpoint_x_low.x)
    {
        hitpoint_x_low = poss_hitpoint_low;
    }
    calculate_hitpoint_on_trajectory(lidar_1_high, lidar_3_low, &poss_hitpoint_low);
    if (poss_hitpoint_low.x > hitpoint_x_low.x)
    {
        hitpoint_x_low = poss_hitpoint_low;
    }    calculate_hitpoint_on_trajectory(lidar_2_high, lidar_3_low, &poss_hitpoint_low);
    if (poss_hitpoint_low.x > hitpoint_x_low.x)
    {
        hitpoint_x_low = poss_hitpoint_low;
    }

    *vehicle_region_x_high = hitpoint_x_high;
    *vehicle_region_x_low = hitpoint_x_low;
}
void rocket_hit_area(int  pix_min_1, int  pix_max_1, float  rocket_dist_1, int pix_min_2, int pix_max_2, float rocket_dist_2, int pix_min_3, int pix_max_3, float rocket_dist_3, glm::vec3* vehicle_hit_highest_point_x, glm::vec3* vehicle_hit_lowest_point_x)
{

    double wall_angle_1 = global_variables["lidar 1 perde acisi"];

    glm::vec3 lidar_position = glm::vec3();
    lidar_position.x = global_variables["arac uzunlugu"] / 2;
    lidar_position.y = global_variables["lidar 1 yuksekligi "];
    lidar_position.z = 0.0f;

    glm::vec3 lidar_1_point_low;
    glm::vec3 lidar_1_point_high;


    pixToPositionXYZ(pix_min_1, pix_max_1, rocket_dist_1, wall_angle_1, lidar_position, &lidar_1_point_low, &lidar_1_point_high);


    double wall_angle_2 = global_variables["lidar 2 perde acisi"];
    lidar_position = glm::vec3();
    lidar_position.x = global_variables["arac uzunlugu"] / 2;
    lidar_position.y = global_variables["lidar 2 yuksekligi "];
    lidar_position.z = 0.0f;

    glm::vec3 lidar_2_point_low;
    glm::vec3 lidar_2_point_high;


    pixToPositionXYZ(pix_min_2, pix_max_2, rocket_dist_2, wall_angle_2, lidar_position, &lidar_2_point_low, &lidar_2_point_high);

    double wall_angle_3 = global_variables["lidar 3 perde acisi"];
    lidar_position = glm::vec3();
    lidar_position.x = global_variables["arac uzunlugu"] / 2;
    lidar_position.y = global_variables["lidar 3 yuksekligi "];
    lidar_position.z = 0.0f;

    glm::vec3 lidar_3_point_low;
    glm::vec3 lidar_3_point_high;


    pixToPositionXYZ(pix_min_3, pix_max_3, rocket_dist_3, wall_angle_3, lidar_position, &lidar_3_point_low, &lidar_3_point_high);

    glm::vec3 vehicle_region_low = glm::vec3();
    glm::vec3 vehicle_region_high = glm::vec3();

    vehicleHitRegionCalculation(lidar_1_point_low, lidar_1_point_high, lidar_2_point_low, lidar_2_point_high, lidar_3_point_low, lidar_3_point_high, &vehicle_region_high, &vehicle_region_low);
    //std::vector<glm::vec3> vec;
    //vec.push_back(vehicle_region_high);
    //vec.push_back(vehicle_region_low);

    *vehicle_hit_highest_point_x = vehicle_region_high;
    *vehicle_hit_lowest_point_x = vehicle_region_low;
}

