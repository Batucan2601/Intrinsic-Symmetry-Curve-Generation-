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
std::vector<unsigned int> Geodesic_min_dijkstra(TrilateralMesh* m, int& number_of_points, std::vector<unsigned int> average_geodesic_function,
float sweep_percentage,float tau, bool is_color);
std::vector<unsigned int> Geodesic_avg_dijkstra_modified(TrilateralMesh* m, int& no_of_points, float sweep_percentage, int N_ring, bool is_color);