#pragma once
#include "MeshFactory.h"
#include "Mesh.h"

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

FuzzyGeodesicList FuzzyGeodesic_calculateFuzzyGedoesic(Mesh* m, int startIndex, int endIndex, float sigma);
void FuzzyGeodesic_colorMesh(Mesh* m , const FuzzyGeodesicList& fuzzyList , );