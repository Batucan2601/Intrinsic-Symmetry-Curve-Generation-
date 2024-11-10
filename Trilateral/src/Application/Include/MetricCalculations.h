#include "TrilateralMesh.h"
float Metric_get_geodesic_cost(TrilateralMesh* m, unsigned int point_index1 , unsigned int calculated_index1_correspondence , bool isNormalized );
float Metric_get_geodesic_cost_with_list(TrilateralMesh* m, std::vector<unsigned int> point_indices, std::vector<unsigned int> calculated_index_correspondence_list );
void Metric_set_gaussian(TrilateralMesh* m, int gaussian_point, float gaussian_dist);
float Metric_get_gaussian_dist(TrilateralMesh* m);
int Metric_get_gaussian_point_no(TrilateralMesh* m);
void Metric_set_N(int N);
void Metric_write_to_file(TrilateralMesh* m, const std::string& file_name);