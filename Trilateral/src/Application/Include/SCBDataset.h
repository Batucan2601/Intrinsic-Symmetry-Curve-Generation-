#pragma once 
#include "TrilateralMesh.h"
#include "Skeleton.h"
typedef enum
{
	SCAPE,
	TOSCA,
	WATERTIGHT,
	KIDS,
	PRINCETON,
}DATASET;

void SCB_select_dataset(DATASET cur_dataset);
void SCB_read_SCAPE();
void SCB_read_WaterTight();
void SCB_read_Princeton();
void SCB_read_TOSCA();
void SCB_read_KIDS();
void SCB_select_mesh(TrilateralMesh& m, int meshNo, Skeleton& skeleton);


