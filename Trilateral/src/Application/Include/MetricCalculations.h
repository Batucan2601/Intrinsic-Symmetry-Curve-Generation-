#include "TrilateralMesh.h"
float Metric_get_geodesic_cost(TrilateralMesh* m, unsigned int point_index1 , unsigned int calculated_index1_correspondence , bool isNormalized );
float Metric_get_geodesic_cost_with_list(TrilateralMesh* m, std::vector<unsigned int> point_indices, std::vector<unsigned int> calculated_index_correspondence_list );
void Metric_write_to_file(TrilateralMesh* m, const std::string& file_name);