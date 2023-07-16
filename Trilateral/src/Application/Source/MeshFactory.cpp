#pragma once 
#include "../Include/MeshFactory.h"
#include <iostream>
#include "glm/glm.hpp"
#include <glm/ext/matrix_transform.hpp>
MeshFactory::MeshFactory()
{
	std::cout << "hello" << std::endl;
}

MeshFactory::~MeshFactory()
{
	std::cout << "bye" <<std::endl;

}

void MeshFactory::add_mesh(Mesh& m)
{
	mesh_vec.push_back(m);
	
	//add indices
	//get the size of the points so far  and add it throughout the indices
	int size_so_far = points.size();
	for (size_t i = 0; i < m.triangles.size(); i++)
	{
		point_indices.push_back( (size_so_far + m.triangles[i]));
	}
	//add points
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		points.push_back(m.vertices[i]);
		colors.push_back(m.colors[i]);

	}
	
}
void MeshFactory::buffer_meshes()
{
	// 1 mesh part 
	std::vector<float> float_points_vec; 
	for (int i = 0; i < points.size(); i++ )
	{
		float_points_vec.push_back(points[i].x);
		float_points_vec.push_back(points[i].y);
		float_points_vec.push_back(points[i].z);

		float_points_vec.push_back(colors[i].x );
		float_points_vec.push_back(colors[i].y);
		float_points_vec.push_back(colors[i].z);
	}
		
	glBindBuffer(GL_ARRAY_BUFFER, 1); // VBO will lok into that later 
	glBufferData(GL_ARRAY_BUFFER, float_points_vec.size()    * sizeof(float), &float_points_vec[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, point_indices.size()  * sizeof(int), &point_indices[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 3); // VBO will lok into that later 
	float_points_vec.clear();
	for (int i = 0; i < mesh_point_pairs.size(); i++ )
	{
		for (int j = 0; j < mesh_point_pairs[i].point_pairs.size(); j++)
		{
			float_points_vec.push_back(mesh_point_pairs[i].point_pairs[j]);

		}
		
	}
	if(mesh_point_pairs.size( )> 0 )
	glBufferData(GL_ARRAY_BUFFER, float_points_vec.size() * sizeof(float), &float_points_vec[0], GL_STATIC_DRAW);


}

void MeshFactory::draw_meshes()
{
	glDrawElements(GL_TRIANGLES , point_indices.size() - correspondence_lines.size(), GL_UNSIGNED_INT ,0 );
	

}



void MeshFactory::remove_all()
{
	point_indices.clear();
	points.clear();
	colors.clear();
	
}
//void MeshFactory::create_lines(int mesh_no1, int mesh_no2 , int  no_of_lines_to_draw)
//{
//	srand(time(0));
//	no_of_lines = no_of_lines_to_draw;
//	Mesh* m1 = &mesh_vec[mesh_no1];
//	Mesh* m2 = &mesh_vec[mesh_no2];
//	glm::mat4 m1_model = m1->model_mat;
//	glm::mat4 m2_model = m2->model_mat;
//	
//	std::vector<glm::vec4> mesh_1_random_vertices;
//	std::vector<glm::vec4> mesh_2_random_vertices;
//	//generate random points for both of the meshes
//	for (size_t i = 0; i < no_of_lines; i++)
//	{
//		int random_index = rand() % m1->vertices.size();
//		glm::vec3 temp_point = m1->vertices[random_index];
//		//model mat
//		glm::vec4 temp_point_new =  projection * view * m1_model * glm::vec4(temp_point.x, temp_point.y, temp_point.z, 1.0f);
//		temp_point_new /= temp_point_new.w;
//		points.push_back(glm::vec3( temp_point_new.x , temp_point_new.y , temp_point_new.z) );
//
//
//		int random_index_2 = rand() % m2->vertices.size() ;
//		glm::vec3 temp_point_2 = m2->vertices[random_index_2];
//		//model mat
//		glm::vec4 temp_point_new_2 =  projection * view * m2_model * glm::vec4(temp_point_2.x, temp_point_2.y, temp_point_2.z, 1.0f);
//		temp_point_new_2 /= temp_point_new_2.w;
//		points.push_back(glm::vec3(temp_point_new_2.x, temp_point_new_2.y, temp_point_new_2.z));
//	}
//
//	// draw lines 
//}
//void MeshFactory::paint_points(int mesh_no1, int mesh_no2, int  no_of_lines_to_draw)
void MeshFactory::add_lines_between_meshes(int mesh_index1, int mesh_index2, int p1, int p2)
{
	correspondence_lines.push_back( std::make_pair( mesh_index1 , p1 ));
	correspondence_lines.push_back(std::make_pair( mesh_index2, p2));
}
void MeshFactory::get_camera_and_projection(glm::mat4 view_ , glm::mat4  projection_)
{
	view = view_;
	projection = projection_;
}
void MeshFactory::add_all()
{
	std::vector<Mesh> temp_mesh_vector = mesh_vec;
	mesh_vec.clear();
	for (size_t i = 0; i < temp_mesh_vector.size(); i++) // for every mesh  
	{
		add_mesh(temp_mesh_vector[i]);
	}
	buffer_meshes();
}

void MeshFactory::remove_mesh(int mesh_no)
{
	int mesh_size = mesh_vec[mesh_no].vertices.size();
	for (size_t i = 0; i < mesh_no; i++)
	{

	}
	mesh_vec.erase(mesh_vec.begin() + mesh_no);

}

void MeshFactory::draw_mesh(int mesh_no)
{
	int point_size_down_limit = 0;
	for (size_t i = 0; i < mesh_no; i++)
	{
		point_size_down_limit  += mesh_vec[i].triangles.size();
	}
	int point_size_up_limit = point_size_down_limit + mesh_vec[mesh_no].triangles.size();
	glDrawElements(GL_TRIANGLES, mesh_vec[mesh_no].triangles.size()   , GL_UNSIGNED_INT,  (const void * ) ( point_size_down_limit  * (sizeof(unsigned int) ))  );

}


