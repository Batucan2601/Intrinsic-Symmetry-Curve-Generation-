#include "../Include/TrilateralDescriptor.h"
#include "../Include/glm/glm.hpp"

void TrilateralDescriptor_generate_mesh_inside(TrilateralMesh* m, TrilateralDescriptor& desc)
{
    //we need to create a new mesh and adjacencies
    TrilateralMesh m_inside;
    std::vector<unsigned int > mesh_to_mesh_map(m->vertices.size(), -1);
    //add the points
    for (size_t i = 0; i < desc.visited_indices.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.visited_indices[i]]);
        mesh_to_mesh_map[desc.visited_indices[i]] = m_inside.vertices.size() - 1;
    }
    for (size_t i = 0; i < desc.path_1_2.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.path_1_2[i]]);
        mesh_to_mesh_map[desc.path_1_2[i]] = m_inside.vertices.size() - 1;

    }
    for (size_t i = 0; i < desc.path_1_3.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.path_1_3[i]]);
        mesh_to_mesh_map[desc.path_1_3[i]] = m_inside.vertices.size() - 1;

    }
    for (size_t i = 0; i < desc.path_2_3.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.path_2_3[i]]);
        mesh_to_mesh_map[desc.path_2_3[i]] = m_inside.vertices.size() - 1;
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
    for (size_t i = 0; i < m_inside.triangles.size(); i += 3)
    {
        int index1 = m_inside.triangles[i];
        int index2 = m_inside.triangles[i + 1];
        int index3 = m_inside.triangles[i + 2];
        float len12 = glm::length(m_inside.vertices[index1] -  m_inside.vertices[index2]);
        float len13 = glm::length(m_inside.vertices[index1] -  m_inside.vertices[index3]);
        float len23 = glm::length(m_inside.vertices[index2] -  m_inside.vertices[index3]);
        std::pair<int, float> pair;
        pair.first = index2;
        pair.second = len12;
        m_inside.adjacenies[index1].push_back(pair);
        pair.first = index1;
        m_inside.adjacenies[index2].push_back(pair);
        pair.second = len13;
        m_inside.adjacenies[index1].push_back(pair);
        pair.first = index3;
        m_inside.adjacenies[index3].push_back(pair);
        pair.second = len23;
        m_inside.adjacenies[index2].push_back(pair);
        pair.first = index2;
        m_inside.adjacenies[index3].push_back(pair);
    }
    desc.m_inside = m_inside;

}

TrilateralMesh TrilateralDescriptor_generate_mesh_with_resolution(TrilateralMesh* m, TrilateralDescriptor& desc,  int res )
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

        std::pair<int, float>pair(p_new_index,glm::length(p1-p_new));
        m_inside->adjacenies[index1].push_back(pair);
        
        pair.second = glm::length(p2 - p_new);
        m_inside->adjacenies[index2].push_back(pair);

        pair.second = glm::length(p3 - p_new);
        m_inside->adjacenies[index3].push_back(pair);

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