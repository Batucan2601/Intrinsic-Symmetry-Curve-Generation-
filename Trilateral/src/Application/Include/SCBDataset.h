#pragma once 
#include "TrilateralMesh.h"
#include "Skeleton.h"

void SCB_read_SCAPE();
void SCB_read_WaterTight();
void SCB_read_TOSCA();
void SCB_select_mesh(TrilateralMesh& m, int meshNo, Skeleton& skeleton);


