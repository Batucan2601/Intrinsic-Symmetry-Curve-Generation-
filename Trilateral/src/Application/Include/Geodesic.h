#pragma once 
#include <vector>
#include "../Include/Mesh.h"
#include "../Include/FibonacciHeap.h"

std::vector<float> Geodesic_dijkstra(Mesh& m, int point_index);
std::vector<int> Geodesic_dijkstra_predecessors(Mesh& m, int point_index);
std::vector<int> Geodesic_between_two_points(Mesh& m, int p1_index, int p2_index);