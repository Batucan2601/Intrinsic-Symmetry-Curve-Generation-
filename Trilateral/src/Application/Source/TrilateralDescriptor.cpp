#include "../Include/TrilateralDescriptor.h"
#include "../Include/Geodesic.h"
#include "../Include/glm/glm.hpp"
#include "../Include/ROI.h"
#include "../Include/HistogramFunctions.h"


bool TrilateralDescriptor::check_colinearity()
{
    //select the biggest path
    if (this->path_1_2.size() > this->path_1_3.size() && this->path_1_2.size() >= this->path_2_3.size())
    {
        for (size_t i = 0; i < this->path_1_2.size(); i++)
        {
            if (this->path_1_2[i] == this->p3)
            {
                return true; 
            }
        }
    }
    if (this->path_1_3.size() > this->path_1_2.size() && this->path_1_3.size() >= this->path_2_3.size())
    {
        for (size_t i = 0; i < this->path_1_3.size(); i++)
        {
            if (this->path_1_3[i] == this->p2)
            {
                return true;
            }
        }
    }
    if (this->path_2_3.size() > this->path_1_2.size() && this->path_2_3.size() >= this->path_1_3.size())
    {
        for (size_t i = 0; i < this->path_2_3.size(); i++)
        {
            if (this->path_2_3[i] == this->p1)
            {
                return true;
            }
        }
    }
    return false;
}

void TrilateralDescriptor_generate_mesh_inside(TrilateralMesh* m, TrilateralDescriptor& desc)
{
    //we need to create a new mesh and adjacencies
    TrilateralMesh m_inside;
    std::vector< int > desc_to_mesh_map(m->vertices.size(), -1);
    std::vector< int > mesh_to_desc_map;
    //add the points
    for (size_t i = 0; i < desc.visited_indices.size(); i++)
    {
        if (desc_to_mesh_map[desc.visited_indices[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.visited_indices[i]]);
            desc_to_mesh_map[desc.visited_indices[i]] = m_inside.vertices.size() - 1;
            mesh_to_desc_map.push_back(m_inside.vertices.size() - 1);
            m_inside.areas.push_back(desc.visited_indices[i]);
        }
    }
    for (size_t i = 0; i < desc.path_1_2.size(); i++)
    {
        if (desc_to_mesh_map[desc.path_1_2[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.path_1_2[i]]);
            desc_to_mesh_map[desc.path_1_2[i]] = m_inside.vertices.size() - 1;
            mesh_to_desc_map.push_back(m_inside.vertices.size() - 1);
            m_inside.areas.push_back(desc.path_1_2[i]);
        }

    }
    for (size_t i = 0; i < desc.path_1_3.size(); i++)
    {
        if (desc_to_mesh_map[desc.path_1_3[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.path_1_3[i]]);
            desc_to_mesh_map[desc.path_1_3[i]] = m_inside.vertices.size() - 1;
            mesh_to_desc_map.push_back(desc.path_1_3[i]);
            m_inside.areas.push_back(m->areas[desc.path_1_3[i]]);
        }
       

    }
    for (size_t i = 0; i < desc.path_2_3.size(); i++)
    {
        if (desc_to_mesh_map[desc.path_2_3[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.path_2_3[i]]);
            desc_to_mesh_map[desc.path_2_3[i]] = m_inside.vertices.size() - 1;
            mesh_to_desc_map.push_back(desc.path_2_3[i]);
            m_inside.areas.push_back(m->areas[desc.path_2_3[i]]);
        }
    }

    //get triangles 
    for (size_t i = 0; i < m->triangles.size(); i += 3)
    {
        int index1 = m->triangles[i];
        int index2 = m->triangles[i + 1];
        int index3 = m->triangles[i + 2];

        int new_index1 = desc_to_mesh_map[index1];
        int new_index2 = desc_to_mesh_map[index2];
        int new_index3 = desc_to_mesh_map[index3];
        if (new_index1 != -1 && new_index2 != -1 && new_index3 != -1)
        {
            m_inside.triangles.push_back(new_index1);
            m_inside.triangles.push_back(new_index2);
            m_inside.triangles.push_back(new_index3);
        }
    }

    //get adjacencies

    m_inside.adjacenies = std::vector<std::vector<std::pair<int, float>>>(m_inside.vertices.size());
    for (size_t i = 0; i < m->adjacenies.size(); i ++)
    {
        int index = desc_to_mesh_map[i];
        if (index == -1)
        {
            continue; 
        }
        for (size_t j = 0; j <m->adjacenies[i].size(); j++)
        {
            int index_adj = desc_to_mesh_map[m->adjacenies[i][j].first];
            if (index_adj == -1)
            {
                continue;
            }
            std::pair<int, float> pair;
            pair.first = index_adj;
            pair.second = m->adjacenies[i][j].second;
            m_inside.adjacenies[index].push_back(pair); 
        }
    }
    
    desc.m_inside = m_inside;
    desc.m_inside.desc_to_mesh_map = desc_to_mesh_map;
    desc.m_inside.mesh_to_desc_map = mesh_to_desc_map;

}

void TrilateralDescriptor_generate_mesh_with_resolution(TrilateralMesh* m, TrilateralDescriptor& desc,  int res )
{
    TrilateralDescriptor_generate_mesh_inside(m ,desc);
    for (size_t i = 0; i < res; i++)
    {
        TrilateralDescriptor_generate_descriptor_with_resolution(&desc.m_inside, desc);
    }
}

void TrilateralDescriptor_generate_descriptor_with_resolution(TrilateralMesh* m_inside, TrilateralDescriptor& desc)
{
    std::vector<unsigned int> triangles_new = m_inside->triangles;
    for (size_t i = 0; i < m_inside->triangles.size(); i += 3)
    {
        int index1 = m_inside->triangles[i];
        int index2 = m_inside->triangles[i+1];
        int index3 = m_inside->triangles[i+2];
        glm::vec3 p1 = m_inside->vertices[index1];
        glm::vec3 p2 = m_inside->vertices[index2];
        glm::vec3 p3 = m_inside->vertices[index3];
        
        glm::vec3 p_new = (p1 + p2 + p3) / 3.0f;
        int p_new_index = m_inside->vertices.size();
        //add p_new as a point
        m_inside->vertices.push_back(p_new);

        float len_1 = glm::length(p1 - p_new);
        std::pair<int, float>pair(p_new_index, len_1);
        m_inside->adjacenies[index1].push_back(pair);
        
        float len_2 = glm::length(p2 - p_new);
        pair.second = len_2;
        m_inside->adjacenies[index2].push_back(pair);

        float len_3 = glm::length(p3 - p_new);
        pair.second = len_3;
        m_inside->adjacenies[index3].push_back(pair);

        std::vector<std::pair<int, float>> pairs;
        m_inside->adjacenies.push_back(pairs);
        
        pair.first = index1;
        pair.second = len_1;
        m_inside->adjacenies[p_new_index].push_back(pair);
        
        pair.first = index2;
        pair.second = len_2;
        m_inside->adjacenies[p_new_index].push_back(pair);

        pair.first = index3;
        pair.second = len_3;
        m_inside->adjacenies[p_new_index].push_back(pair);


        triangles_new.push_back(p_new_index);
        triangles_new.push_back(index1);
        triangles_new.push_back(index2);

        triangles_new.push_back(p_new_index);
        triangles_new.push_back(index1);
        triangles_new.push_back(index3);

        triangles_new.push_back(p_new_index);
        triangles_new.push_back(index2);
        triangles_new.push_back(index3);

    }
    m_inside->triangles = triangles_new;

}
TrilateralDescriptor  TrilateralDescriptor_create(TrilateralMesh* m, int point_index1, int point_index2, int point_index3, bool is_simplified)
{
    TrilateralDescriptor trilateral_descriptor;//trialteral descriptor 
    //init descriptor
    trilateral_descriptor.p1 = point_index1;
    trilateral_descriptor.p2 = point_index2;
    trilateral_descriptor.p3 = point_index3;

    ROI_trilateral(m, trilateral_descriptor, 10, false);

    trilateral_descriptor.curvature_1_2 = trilateral_descriptor.geodesic_lenght_1_2 / glm::distance(m->vertices[point_index1], m->vertices[point_index2]);
    trilateral_descriptor.curvature_1_3 = trilateral_descriptor.geodesic_lenght_1_3 / glm::distance(m->vertices[point_index1], m->vertices[point_index3]);
    trilateral_descriptor.curvature_2_3 = trilateral_descriptor.geodesic_lenght_2_3 / glm::distance(m->vertices[point_index2], m->vertices[point_index3]);

    trilateral_descriptor.euclidian_lenght_1_2 = glm::distance(m->vertices[point_index1], m->vertices[point_index2]);
    trilateral_descriptor.euclidian_lenght_1_3 = glm::distance(m->vertices[point_index1], m->vertices[point_index3]);
    trilateral_descriptor.euclidian_lenght_2_3 = glm::distance(m->vertices[point_index2], m->vertices[point_index3]);

    trilateral_descriptor.curvature_1_2 = trilateral_descriptor.geodesic_lenght_1_2 / trilateral_descriptor.euclidian_lenght_1_2;
    trilateral_descriptor.curvature_1_3 = trilateral_descriptor.geodesic_lenght_1_3 / trilateral_descriptor.euclidian_lenght_1_3;
    trilateral_descriptor.curvature_2_3 = trilateral_descriptor.geodesic_lenght_2_3 / trilateral_descriptor.euclidian_lenght_2_3;

    trilateral_descriptor.n_ring_area_p1 = get_N_ring_area(m, trilateral_descriptor.p1, 1);
    trilateral_descriptor.n_ring_area_p2 = get_N_ring_area(m, trilateral_descriptor.p2, 1);
    trilateral_descriptor.n_ring_area_p3 = get_N_ring_area(m, trilateral_descriptor.p3, 1);


    //for only brute force research 
    return trilateral_descriptor;
}

void TrilateralDescriptor_write(std::string filename, std::vector<TrilateralDescriptor>& positive_desc, std::vector<TrilateralDescriptor>& negative_desc)
{
    std::ofstream file;                // Create an ofstream object for file output

    // Open the file in write mode
    file.open(filename);

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write some data to the file
    //desc format 1 - 
    for (size_t i = 0; i < positive_desc.size(); i++)
    {
        file << " desc ";
        file << positive_desc[i].p1 << " " << positive_desc[i].p2 << " " << positive_desc[i].p3 << std::endl;
        file << " vertices inside ";
        for (size_t j = 0; j < positive_desc[i].visited_indices.size(); j++)
        {
            file << positive_desc[i].visited_indices[j] << " ";
        }
        file << std::endl;
        file << " path12";
        for (size_t j = 0; j < positive_desc[i].path_1_2.size(); j++)
        {
            file << positive_desc[i].path_1_2[j] << " ";
        }
        file << std::endl;
        file << " path13";
        for (size_t j = 0; j < positive_desc[i].path_1_3.size(); j++)
        {
            file << positive_desc[i].path_1_3[j] << " ";
        }
        file << std::endl;
        file << " path23";
        for (size_t j = 0; j < positive_desc[i].path_2_3.size(); j++)
        {
            file << positive_desc[i].path_2_3[j] << " ";
        }
        file << std::endl;
    }
    file << "negative" << std::endl;
    for (size_t i = 0; i < negative_desc.size(); i++)
    {
        file << " desc ";
        file << negative_desc[i].p1 << " " << negative_desc[i].p2 << " " << negative_desc[i].p3 << std::endl;
        file << " vertices inside ";
        for (size_t j = 0; j < negative_desc[i].visited_indices.size(); j++)
        {
            file << negative_desc[i].visited_indices[j] << " ";
        }
        file << std::endl;
        file << " path12";
        for (size_t j = 0; j < negative_desc[i].path_1_2.size(); j++)
        {
            file << negative_desc[i].path_1_2[j] << " ";
        }
        file << std::endl;
        file << " path13";
        for (size_t j = 0; j < negative_desc[i].path_1_3.size(); j++)
        {
            file << negative_desc[i].path_1_3[j] << " ";
        }
        file << std::endl;
        file << " path23";
        for (size_t j = 0; j < negative_desc[i].path_2_3.size(); j++)
        {
            file << negative_desc[i].path_2_3[j] << " ";
        }
    }
    // Close the file
    file.close();

}

void TrilateralDescriptor_read(std::string filename, std::vector<TrilateralDescriptor>& positive_desc, std::vector<TrilateralDescriptor>& negative_desc)
{
    std::ifstream file(filename);                // Create an ofstream object for file output

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write some data to the file
    //desc format 1 - 
    int pos_size = positive_desc.size();
    std::vector<TrilateralDescriptor> descriptors;
    int negative_start_index = -1;
    bool is_pos_desc = true;
    std::string line;

    TrilateralDescriptor desc;
    while (std::getline(file, line))
    {

        if (line.find("negative") != std::string::npos)
        {
            is_pos_desc = false;
            negative_start_index = descriptors.size();
        }
        if (line.find("desc") != std::string::npos)
        {
            line = line.substr(5);
            std::stringstream ss(line);
            std::vector<int> nums;
            int num;
            while (ss >> num)
            {
                nums.push_back(num);
            }
            desc.p1 = nums[0];
            desc.p2 = nums[1];
            desc.p3 = nums[2];
        }
        if (line.find("vertices inside") != std::string::npos)
        {
            line = line.substr(16);
            std::stringstream ss(line);
            std::vector<unsigned int> visited;
            unsigned int num;
            while (ss >> num)
            {
                visited.push_back(num);
            }
            desc.visited_indices = visited;
        }
        if (line.find("path12") != std::string::npos)
        {
            line = line.substr(7);
            std::stringstream ss(line);
            std::vector<int> visited;
            int num;
            while (ss >> num)
            {
                visited.push_back(num);
            }
            desc.path_1_2 = visited;
        }
        if (line.find("path13") != std::string::npos)
        {
            line = line.substr(7);
            std::stringstream ss(line);
            std::vector<int> visited;
            int num;
            while (ss >> num)
            {
                visited.push_back(num);
            }
            desc.path_1_3 = visited;
        }
        if (line.find("path23") != std::string::npos)
        {
            line = line.substr(7);
            std::stringstream ss(line);
            std::vector<int> visited;
            int num;
            while (ss >> num)
            {
                visited.push_back(num);
            }
            desc.path_2_3 = visited;
            descriptors.push_back(desc);
        }
    }
    // Close the file
    file.close();

    for (size_t i = 0; i < negative_start_index; i++)
    {
        positive_desc.push_back(descriptors[i]);
    }
    for (size_t i = negative_start_index; i < descriptors.size(); i++)
    {
        negative_desc.push_back(descriptors[i]);
    }
    return;
}

//returns cdf normalized
Histogram TrilateralDescriptor_generate_cdf_of_areas(TrilateralMesh* m, TrilateralDescriptor& desc1 , int division_no)
{
    Histogram h = Histogram_triangle_area_w_res( m,desc1,division_no,6);
    Histogram cdf(division_no); 
    float sum = 0;
    for (size_t i = 0; i < division_no; i++)
    {
        sum += h[i];
        cdf[i] = sum; 
    }
    return cdf;
}