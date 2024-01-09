#include "Mesh.h"
float get_geodesic_cost(Mesh* m, unsigned int point_index1 , unsigned int calculated_index1_correspondence );
float get_geodesic_cost_with_list(Mesh* m, std::vector<unsigned int> point_indices, std::vector<unsigned int> calculated_index_correspondence_list );