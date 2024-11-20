#pragma once 
#include <vector>
#include "../Include/TrilateralMesh.h"
#include "../Include/FibonacciHeap.h"
#include "../Include/NLateralDescriptor.h"

bool Geodesic_proximity(TrilateralMesh& m, NLateralDescriptor& desc1, NLateralDescriptor& desc2, float proximity);
std::vector<float> Geodesic_dijkstra(TrilateralMesh& m, int point_index);
std::vector<int> Geodesic_dijkstra_predecessors(TrilateralMesh& m, int point_index);
std::vector<int> Geodesic_between_two_points(TrilateralMesh& m, int p1_index, int p2_index);