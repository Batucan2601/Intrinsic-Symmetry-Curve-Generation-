#pragma once
#include "MeshFactory.h"
#include "TrilateralMesh.h"

typedef struct {
	float fuzziness;
	int pointIndex;
}FuzzyGeodesic;

typedef struct
{
	int startIndex;
	int endIndex;
	std::vector<FuzzyGeodesic> fuzzyGeodesicList;
}FuzzyGeodesicList;



FuzzyGeodesicList FuzzyGeodesic_calculateFuzzyGedoesic(TrilateralMesh* m, int startIndex, int endIndex, float sigma);
float FuzzyGeodesic_FuzzyArea(TrilateralMesh* m , const FuzzyGeodesicList& fuzzyList , bool color = false );
float FuzzyGeodesic_fuzzyGeo(TrilateralMesh* m, int start_index, int end_index, int index_x);