#include "../Include/TrilateralDescriptor.h"
#include "../Include/Geodesic.h"
#include "../Include/glm/glm.hpp"
#include "../Include/ROI.h"


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
    std::vector< int > mesh_to_mesh_map(m->vertices.size(), -1);
    //add the points
    for (size_t i = 0; i < desc.visited_indices.size(); i++)
    {
        if (mesh_to_mesh_map[desc.visited_indices[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.visited_indices[i]]);
            mesh_to_mesh_map[desc.visited_indices[i]] = m_inside.vertices.size() - 1;
        }
    }
    for (size_t i = 0; i < desc.path_1_2.size(); i++)
    {
        if (mesh_to_mesh_map[desc.path_1_2[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.path_1_2[i]]);
            mesh_to_mesh_map[desc.path_1_2[i]] = m_inside.vertices.size() - 1;
        }

    }
    for (size_t i = 0; i < desc.path_1_3.size(); i++)
    {
        if (mesh_to_mesh_map[desc.path_1_3[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.path_1_3[i]]);
            mesh_to_mesh_map[desc.path_1_3[i]] = m_inside.vertices.size() - 1;
        }
       

    }
    for (size_t i = 0; i < desc.path_2_3.size(); i++)
    {
        if (mesh_to_mesh_map[desc.path_2_3[i]] == -1)
        {
            m_inside.vertices.push_back(m->vertices[desc.path_2_3[i]]);
            mesh_to_mesh_map[desc.path_2_3[i]] = m_inside.vertices.size() - 1;
        }
    }

    //get triangles 
    for (size_t i = 0; i < m->triangles.size(); i += 3)
    {
        int index1 = m->triangles[i];
        int index2 = m->triangles[i + 1];
        int index3 = m->triangles[i + 2];

        int new_index1 = mesh_to_mesh_map[index1];
        int new_index2 = mesh_to_mesh_map[index2];
        int new_index3 = mesh_to_mesh_map[index3];
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
        int index = mesh_to_mesh_map[i];
        if (index == -1)
        {
            continue; 
        }
        for (size_t j = 0; j <m->adjacenies[i].size(); j++)
        {
            int index_adj = mesh_to_mesh_map[m->adjacenies[i][j].first];
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