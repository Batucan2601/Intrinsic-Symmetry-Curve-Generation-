#pragma once 
#include <vector>
#include "../Include/TrilateralMesh.h"
#include "../Include/FibonacciHeap.h"
#include "../Include/NLateralDescriptor.h"

bool Geodesic_proximity(TrilateralMesh& m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float proximity);
std::vector<float> Geodesic_dijkstra(TrilateralMesh& m, int point_index);
std::vector<int> Geodesic_dijkstra_predecessors(TrilateralMesh& m, int point_index);
std::vector<int> Geodesic_between_two_points(TrilateralMesh& m, int p1_index, int p2_index);
std::vector<unsigned int> Geodesic_avg_dijkstra(TrilateralMesh* m, int& number_of_points,float sweep_distance,int N_ring, bool is_color);
std::vector<unsigned int> Geodesic_min_dijkstra(TrilateralMesh* m,std::vector<unsigned int> average_geodesic_function,
float sweep_percentage,float tau, bool is_color);
std::vector<unsigned int> Geodesic_avg_dijkstra_modified(TrilateralMesh* m, float sweep_percentage, int N_ring, bool is_color, float& biggest_dijkstra);
std::vector<unsigned int> Geodesic_avg_dijkstra_modified_to_points(TrilateralMesh* m, std::vector<unsigned int> points,
int& no_of_points, float sweep_percentage, int N_ring, bool is_color);

std::vector<unsigned int> Geodesic_find_biggest_AGD(TrilateralMesh* m, float sweep_percentage, float stop_param);

void Geodesic_write_sampled_points(TrilateralMesh* m, std::vector<unsigned int>& agd_points);
void Geodesic_read_sampled_points(TrilateralMesh* m, std::vector<unsigned int>& sampled_points);
unsigned int Geodesic_find_midpoint(TrilateralMesh* m, unsigned int index1, unsigned int index2);
void Geodesic_mid_point_w_AGD(TrilateralMesh* m, unsigned int& p1, unsigned int& p2 , float& biggest_dijkstra);


std::vector<unsigned int> Geodesic_get_max_geodesics(TrilateralMesh* m,  float sweep_percentage  );

std::vector<unsigned int> conv_int_to_unsigned(std::vector<int> vec);
std::vector<unsigned int> Geodesic_generate_secondary_curve_w_midpoints(TrilateralMesh* m, unsigned int& midpoint1, unsigned int& midpoint2);
std::vector<unsigned int> Geodesic_generate_multiple_secondary_curve(TrilateralMesh* m, unsigned int& midpoint1, unsigned int& midpoint2);

std::vector<unsigned int> Geodesic_generate_secondary_curve(TrilateralMesh* m, unsigned int& midpoint1, unsigned int& midpoint2);
unsigned int Geodesic_send_ray_get_counterpart(TrilateralMesh* m, unsigned int& midpoint1);


bool Geodesic_path_intersection(TrilateralMesh* m, std::vector<unsigned int>& path1, std::vector<unsigned int>& path2, unsigned int& no_of_times);
void Geodesic_color_path(TrilateralMesh* m, unsigned int p1, unsigned int p2);
void Geodesic_color_midpoints(TrilateralMesh* m);
void Geodesic_color_according_to_midpoints(TrilateralMesh* m);
unsigned int Geodesic_get_midpoint_from_path(TrilateralMesh* m, unsigned int p1, unsigned int p2);
