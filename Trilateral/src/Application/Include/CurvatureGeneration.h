#include <vector>
#include "glm/glm.hpp"
#include "TrilateralMesh.h"
struct Curvature
{
	std::vector<unsigned int> curvature_main_points; //for points founed with geodesic
	std::vector<glm::vec3> points; 
	float curvature_quality; 
};
struct Curve
{
	std::vector<unsigned int> curve_path; 
	float length; 
};

Curvature CurvatureGeneration_generate(TrilateralMesh* m, float merge_distance_param);
float CurvatureGeneration_curvature_quality(TrilateralMesh* m, Curvature& curv);
std::vector<Curve> CurvatureGeneration_generate_curve_paths(TrilateralMesh* m);
float  CurvatureGeneration_get_curve_length(TrilateralMesh* m, Curve& curv  );
