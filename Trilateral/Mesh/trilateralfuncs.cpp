Histogram Trilateral_triangle_area_with_resolution(TrilateralMesh* m , TrilateralDesc& desc, int res = 1  )
{
    TrilateralMesh m_inside = TrilateralDescriptor_generate_mesh_inside(m ,  desc );
    for (size_t i = 0; i < m_inside.size(); i+=3)
    {
        
    }

    



}

TrilateralMesh TrilateralDescriptor_generate_mesh_inside(TrilateralMesh*m , TrilateralDescriptor& desc )
{
    //we need to create a new mesh and adjacencies
    TrilateralMesh m_inside;
    std::vector<unsigned int > mesh_to_mesh_map(  m->.vertices.size(),-1 );
    //add the points
    for (size_t i = 0; i < desc.visited_vertices.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.visited_vertices[i]]);
        mesh_to_mesh_mapping[desc.visited_vertices[i]] =m_inside_vertices.size()-1; 
    }
    for (size_t i = 0; i < desc.path_1_2.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.path_1_2[i]]);
        mesh_to_mesh_mapping[desc.path_1_2[i]] =m_inside_vertices.size()-1; 

    }
    for (size_t i = 0; i < desc.path_1_3.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.path_1_3[i]]);
        mesh_to_mesh_mapping[desc.path_1_3[i]] =m_inside_vertices.size()-1; 

    }
    for (size_t i = 0; i < desc.path_2_3.size(); i++)
    {
        m_inside.vertices.push_back(m->vertices[desc.path_2_3[i]]);
        mesh_to_mesh_mapping[desc.path_2_3[i]] =m_inside_vertices.size()-1; 
    }

    //get triangles 
    for (size_t i = 0; i < m->triangles.size(); i+=3)
    {
        int index1 = m->triangles[i];
        int index2 = m->triangles[i+1];
        int index3 = m->triangles[i+2];

        int new_index1 = mesh_to_mesh_index[index1];
        int new_index2 = mesh_to_mesh_index[index2];
        int new_index3 = mesh_to_mesh_index[index3];
        if( new_index1 != -1 && new_index2 != -1  && new_index3 != -1  )
        {
            m_inside.triangles.push_back( new_index1);
            m_inside.triangles.push_back( new_index2);
            m_inside.triangles.push_back( new_index3);
        }
    }
    
    //get adjacencies
    m_inside.adjacencies = std::vector<std::pair<int,float>>( m_inside.vertices.size());
    for (size_t i = 0; i < m_inside.triangles.size(); i+=3)
    {
        int index1 = m_inside.triangles[i];
        int index2 = m_inside.triangles[i+1];
        int index3 = m_inside.triangles[i+2];
        float len12 = glm::length(m_inside.vertices[index1] ,m_inside.vertices[index2] );  
        float len13 = glm::length(m_inside.vertices[index1] ,m_inside.vertices[index3] );
        float len23 = glm::length(m_inside.vertices[index2] ,m_inside.vertices[index3] );
        std::pair<int,float> pair; 
        pair.first = index2;
        pair.second = len12;
        m_inside.adjacencies[index1].push_back(pair);
        pair.first = index1;
        m_inside.adjacencies[index2].push_back(pair);
        pair.second = len13;
        m_inside.adjacencies[index1].push_back(pair);
        pair.first = index3;
        m_inside.adjacencies[index3].push_back(pair);
        pair.second = len23;
        m_inside.adjacencies[index2].push_back(pair);
        pair.first = index2;
        m_inside.adjacencies[index3].push_back(pair);

    }

    return m_inside; 
    

}
