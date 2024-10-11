#pragma once 
#include <vector>
#include "Histogram.h"

enum PointStatus
{
	INSIDE,
	OUTSIDE,
	EDGE,
	MIDPOINT,
};
enum TrilateralGeometry
{
	ENCLOSED, // expected case
	LINEAR, // no closed space
	ENCLOSED_TWO // two points enclose 
};
struct TrilateralDescriptor
{
	double area; // ROI
	double total_length;
	double geodesic_lenght_1_2;//  geodesic length between 1 - 2
	double geodesic_lenght_1_3;//  length between 1 - 3
	double geodesic_lenght_2_3;//  length between 2 - 3
	double euclidian_lenght_1_2;//  euclidian length between 1 - 2
	double euclidian_lenght_1_3;//  euclidian length between 1 - 2
	double euclidian_lenght_2_3;//  euclidian length between 1 - 2
	double curvature_1_2;
	double curvature_1_3;
	double curvature_2_3;
	unsigned int p1; //point indices
	unsigned int p2;
	unsigned int p3;
	// extras
	float n_ring_area_p1;
	float n_ring_area_p2;
	float n_ring_area_p3;
	// histogram
	Histogram histogram;
	std::vector<unsigned int> visited_indices;
	std::vector<int> path_1_2;
	std::vector<int> path_1_3;
	std::vector<int> path_2_3;
	TrilateralMesh m_inside;
	TrilateralGeometry geo;
	TrilateralDescriptor();
	bool check_colinearity();
};

void TrilateralDescriptor_generate_mesh_inside(TrilateralMesh* m, TrilateralDescriptor& desc);
void TrilateralDescriptor_generate_descriptor_with_resolution(TrilateralMesh* m_inside, TrilateralDescriptor& desc);
void TrilateralDescriptor_generate_mesh_with_resolution(TrilateralMesh* m, TrilateralDescriptor& desc, int res);