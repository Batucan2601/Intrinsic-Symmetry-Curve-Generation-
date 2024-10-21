#include "../Include/KIDS.h"
#include "../Include/DominantSymmetry.h"
#include "../Include/DvorakEstimatingApprox.h"

std::vector<TrilateralMesh> Kids_dataset;
std::vector<std::vector<unsigned int>>  KIDS_dataset_gaussian_indices;

void KIDS_read_meshes()
{
    // total of 15 + 15 meshes 
    for (size_t i = 0; i < 30; i++)
    {
        std::string path("../../Trilateral/Mesh/off/");
        std::string isometry_batch_no("000" + std::to_string((i / 15) + 1));
        std::string isometry_no(std::to_string(i % 15 + 1));
        // read the meshes.
        path = path + isometry_batch_no + ".isometry." + isometry_no + ".off";
        TrilateralMesh m((char*)path.c_str());
        //read the symmetry format
        read_symmetry_format((char*)"../../Trilateral/TrilateralMesh/off/sym.txt", &m);
        Kids_dataset.push_back(m);
    }
}

void KIDS_dom_sym_generate_or_read_planes(float convergence_ratio)
{
	for (size_t i = 0; i < Kids_dataset.size(); i++)
	{
		TrilateralMesh* m = &Kids_dataset[i];
		snprintf(NULL, 0, "%f", convergence_ratio);
		std::string name = m->file_name + std::to_string(convergence_ratio) + ".pln";
		std::string path("../../Trilateral/Mesh/off/DomSym/");
		path = path + name;
		Plane plane;
		// try to read the plane
		bool is_already = dom_sym_read_plane(m, plane, path);
		if (!is_already )
		{
			plane = generate_dominant_symmetry_plane(m, convergence_ratio);
			dom_sym_write_plane(m, plane, path);
		}
	}
}

void KIDS_generate_gaussians(int no_of_gaussian , float sweep_distance )
{
	for (size_t i = 0; i < Kids_dataset.size(); i++)
	{
		TrilateralMesh* m = &Kids_dataset[i];
		std::vector<DvorakPairs> dvorak_pairs = dvorak_extraction_of_significant_points( m , no_of_gaussian);
		dvorak_pairs = dvorak_distance_sweep(m, dvorak_pairs, sweep_distance);
		
		std::vector<unsigned int> gaussian_indices_i;
		KIDS_dataset_gaussian_indices.push_back(gaussian_indices_i);
		
		for (size_t j = 0; j < dvorak_pairs.size(); j++)
		{
			int index = dvorak_pairs[j].p_index;
			KIDS_dataset_gaussian_indices[i].push_back(index);
		}
	}
}


TrilateralMesh* KIDS_select_mesh(  int mesh_no)
{
	return &Kids_dataset[mesh_no];
}