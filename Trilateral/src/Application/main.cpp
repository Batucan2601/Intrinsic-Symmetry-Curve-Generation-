/*#include <GL/glew.h>
#include <GLFW/glfw3.h>*/
#include <iostream>
#include "raylib.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "imgui/imgui_stdlib.h"
#include "imgui/imgui.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "Include/Sampling.h"
#include "Include/Mesh_imgui.h"
#include "Include/ImguiMenuBar.h"
#include "Include/Ray.h"
#include "rlImGui/rlImGui.h"


//include prototypes
#include "Include/Prototypes.h"
//#include "TrilateralMesh.h"
#include "Include/MeshFactory.h"
#include "Include/Shader.h"
#include "Include/TrilateralMap.h"
#pragma region  glcall 
#define ASSERT(x) if (!(x)) __debugbreak();
#define GLCall(x) GLClearError();x;ASSERT(GLLogCall(#x, __FILE__, __LINE__));
/*static void GLClearError()
{
    while (glGetError() != GL_NO_ERROR);
}

static bool GLLogCall(const char* function, const char* file, int line)
{
    while (GLenum error = glGetError())
    {
        std::cout << "[OpenGL Error] (" << error << "):" << function <<
            " " << file << ":" << line << std::endl;
        return false;
    }
    return true;
}
#pragma endregion
#pragma region shader compilation 
struct ShaderProgramSource
{
    std::string VertexSource;
    std::string FragmentSource;
};
static ShaderProgramSource ParseShader(const std::string& filepath)
{
    std::ifstream stream(filepath);

    enum class ShaderType {
        NONE = -1, VERTEX = 0, FRAGMENT = 1
    };

    ShaderType type = ShaderType::NONE;
    std::string line;
    std::stringstream ss[2];
    while (getline(stream, line))
    {
        if (line.find("#shader") != std::string::npos)
        {
            if (line.find("vertex") != std::string::npos)
                //set mode to vertex
            {
                type = ShaderType::VERTEX;
            }
            else if (line.find("fragment") != std::string::npos)
                //set mode to fragment
            {
                type = ShaderType::FRAGMENT;
            }
        }
        else
        {
            ss[(int)type] << line << "\n";
        }
    }

    return { ss[0].str(), ss[1].str() };

}
static int CompileShader(unsigned int type, const std::string& source)
{
    unsigned int id = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(id, 1, &src, nullptr);
    glCompileShader(id);

    // TODO: Error handling
    int result;
    glGetShaderiv(id, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE)
    {
        int length;
        glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
        char* message = (char*)alloca(length * sizeof(char));
        glGetShaderInfoLog(id, length, &length, message);
        std::cout << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << " shader!" << std::endl << message << std::endl;
        glDeleteShader(id);
        return 0;
    }
    return id;
}

static unsigned int CreateShader(const std::string& vertexShader, const std::string& fragmentShader)
{
    GLCall(unsigned int program_index = glCreateProgram());
    GLCall(unsigned int vertex_shader = CompileShader(GL_VERTEX_SHADER, vertexShader));
    GLCall(unsigned int fragment_shader = CompileShader(GL_FRAGMENT_SHADER, fragmentShader));

    GLCall(glAttachShader(program_index, vertex_shader));
    GLCall(glAttachShader(program_index, fragment_shader));
    GLCall(glLinkProgram(program_index));
    GLCall(glValidateProgram(program_index));

    GLCall(glDeleteShader(vertex_shader));
    GLCall(glDeleteShader(fragment_shader));

    return program_index;
}
unsigned int set_shader(std::string path)
{
    ShaderProgramSource vertex_shader_source = ParseShader(path);
    unsigned int vertex_shader_program = CreateShader(vertex_shader_source.VertexSource, vertex_shader_source.FragmentSource);
    GLCall(glUseProgram(vertex_shader_program));
    // enable increasing point size
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    glUseProgram(vertex_shader_program);
    return vertex_shader_program;
}*/
#pragma endregion

#pragma region GLFW variables
#pragma endregion 
#pragma region mvp variables
//code of the program
unsigned int vertex_shader_program; 
//view  projection direction
glm::mat4 proj;
glm::mat4 view;
glm::mat4 model;
glm::vec3 direction;
const float cameraSpeed = 10.0f;// speed of the camera
bool camera_mouse_lock = true; // mouse callback lock for up and down movement
// mvp matrix
glm::mat4 MVP;
//camera variables
glm::vec3 cameraPos = glm::vec3(-10.0f, 1.25f, 0.0f); // where camera starts
glm::vec3 cameraFront = glm::vec3(1.0f, 0.0f, 0.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
#pragma endregion 
GLFWwindow* window;
bool is_e_pressed = false; // edit mode enter
std::map<std::string, float> global_variables;


// simulation elements
bool is_single_simulation_started = false;
bool is_multiple_random_simulation_started = false;
bool is_multiple_rocket_popup_open = false; 
bool is_multiple_simulation_popup_started = false; 

//multiple simulation deprecated 
float no_of_rockets = 1000;
float min_x_range = 0.0f, min_y_range = 0.0f, min_z_range = -500.0f;
float max_x_range = 5.0f, max_y_range = 2.5f, max_z_range = 0.0f; 
// multiple simulation anew
float grid_x_partition_size = 10.0f; float  tank_x_partition_size = 10.0f;
float grid_y_partiion_size = 10.0f; float tank_y_partition_size = 10.0f;
float distance_from_tank = -30.0f; 
float degree_from_left_right = 5.0f; //tank degree coming from left and right
float degree_from_up_down = 45.0f; 
int sensor_size_min = 1; 
int sensor_size_max = 10;
int sensor_size_increment = 1; 

float radius;
std::string* record_text_name; // name of the text to be recorded with  
std::ofstream* recorded_file;
//iterations
float iteration_speed = 0.001f;
float iteration_no = 1;
int total_hit = 0;
int total_iteration = 0;


int main(void) 
{

    // first of all initialise the hashtable values to default
    // glfw: initialize and configure
    // ------------------------------
    /*glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);*/



    // glfw window creation
    // --------------------
    /*const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    window = glfwCreateWindow(mode->width, mode->height, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);*/

    InitWindow(500,500 , " Trialteral");
    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    /*if (glewInit() != GLEW_OK)
    {
        std::cout << "problems" << std::endl;
    }
    // variables for the camera movement
    glEnable(GL_DEPTH_TEST);*/
    proj = glm::perspective(glm::radians(90.0f), (float)mode->width / (float)mode->height, 0.1f, 1000.0f);
    view = glm::mat4(1.0f);
    model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, 0.0));
    // build and compile our shader program
   // ------------------------------------
#pragma region shaders
    //TrilateralShader default_shader("../../Trilateral/shaders/MeshVertex.shader" , "../../Trilateral/shaders/MeshFragment.shader");

    //default_shader.use();
#pragma endregion 
    /*unsigned int VBO, VAO, IBO;
    unsigned int VBO_matching_points , VAO_matching_points , IBO_matchig_points; 
    glGenVertexArrays(1, &VAO);
    glGenVertexArrays(1, &VAO_matching_points);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &IBO);
    glGenBuffers(1, &VBO_matching_points);
    glGenBuffers(1, &IBO_matchig_points);
    unsigned int VBO_pairs; 
    glGenBuffers(1, &VBO_pairs);*/

#pragma region cs589 before loop init

    TrilateralMesh m1((char*)"../../Trilateral/Mesh/off/0001.isometry.8.off");


    
    
    
    MeshFactory mesh_fac;
    mesh_fac.add_mesh(m1);



   
    // END OF SUGGESTED MESHES
    /*glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    int selected_mesh = 0;
    mesh_fac.buffer_meshes();*/
    
     


    /*glBindVertexArray(VAO_matching_points);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_matching_points);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO_matchig_points);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(VAO);*/

    
#pragma endregion cs589 before loop init
    // render loop
    // -----------
    //imgui_setup();
    rlImGuiSetup(true);
    while (!WindowShouldClose())
    {


        // input
        // -----
        //poll_keys();

        //imgui_new_frame();
        // render
        // ------

        /*glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);*/

        // calculate movement
        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);


#pragma region  cs589e
        //mesh shader 
        //default_shader.use();
        //classic VAO 

        
        BeginDrawing();
        rlImGuiBegin();

        // show ImGui Content
        bool open = true;
        ImGui::ShowDemoWindow(&open);

        // end ImGui Content
        rlImGuiEnd();
        DrawLine(0, 0, 100, 100, RED);
        EndDrawing();
        
        //glBindVertexArray(VAO);

        /*for (size_t i = 0; i < mesh_fac.mesh_vec.size(); i++)
        {
            glm::mat4 model = mesh_fac.mesh_vec[i].model_mat;
            MVP = proj * view * model;
            mesh_fac.mesh_vec[i].MVP = MVP;
            glUniformMatrix4fv(glGetUniformLocation(default_shader.ID, "u_MVP"), 1, GL_FALSE, &MVP[0][0]);
            mesh_fac.mesh_vec[0].model_mat = model; 
            mesh_fac.draw_mesh(i);
        }

        glBindVertexArray(VAO_matching_points);
        for (size_t i = 0; i < mesh_fac.mesh_vec.size(); i++)
        {
            if (mesh_fac.mesh_vec[i].calculated_symmetry_pairs.size() > 0)
            {
                glm::mat4 model = mesh_fac.mesh_vec[i].model_mat;
                MVP = proj * view * model;
                mesh_fac.mesh_vec[i].MVP = MVP;
                glUniformMatrix4fv(glGetUniformLocation(default_shader.ID, "u_MVP"), 1, GL_FALSE, &MVP[0][0]);
                glLineWidth(5.0f);
                glDrawArrays(GL_LINES, 0, mesh_fac.mesh_vec[i].calculated_symmetry_pairs.size() * 2);
                glLineWidth(1.0f);
            }
            
        }
       glBindVertexArray(mesh_fac.skeleton_VAO);
        if (mesh_fac.mesh_skeleton_vec.skeleton_points.size() > 0)
        {
            glm::mat4 model = mesh_fac.mesh_vec[0].model_mat;
            MVP = proj * view * model;
            mesh_fac.mesh_vec[0].MVP = MVP;
            glUniformMatrix4fv(glGetUniformLocation(default_shader.ID, "u_MVP"), 1, GL_FALSE, &MVP[0][0]);

            glLineWidth(5.0f);
            //glDrawArrays(GL_LINES, 0, mesh_fac.mesh_skeleton_vec.skeleton_points.size()/2 );
            glDrawElements(GL_LINES, mesh_fac.mesh_skeleton_vec.skeleton_indices.size(), GL_UNSIGNED_INT, 0);
            //glDraw(GL_LINES, 0, mesh_fac.mesh_skeleton_vec.skeleton_points.size()/2 );
            glLineWidth(1.0f);

        }
        glBindVertexArray(mesh_fac.mesh_vec[0].vao_normals);
        if (is_normals_shown)
        {
            glm::mat4 model = mesh_fac.mesh_vec[0].model_mat;
            MVP = proj * view * model;
            mesh_fac.mesh_vec[0].MVP = MVP;
            glUniformMatrix4fv(glGetUniformLocation(default_shader.ID, "u_MVP"), 1, GL_FALSE, &MVP[0][0]);

            glLineWidth(5.0f);
            glDrawArrays(GL_LINES, 0, mesh_fac.mesh_vec[0].normals_display.size() / 6);
            glLineWidth(1.0f);
        }
        if (is_skeletalNLateral_created)
        {
            SkeletalNLateral_draw( mesh_fac, default_shader.ID);
        }

        glBindVertexArray(0);

        mesh_fac.get_camera_and_projection(view, proj);
       
        if (activate_histogram)
        {
            imgui_histogram(histogram, partition_no);
        }*/
#pragma endregion 
        // ui 
        //imgui_input_window();
                //imgui new frame
        /*imgui_mesh_window(selected_mesh, mesh_fac);
        imgui_selected_mesh_properties_window(selected_mesh, mesh_fac);
        imgui_KIDS_skeleton(selected_mesh, mesh_fac);
        imgui_N_Lateral_Parameters(selected_mesh, mesh_fac);
        imgui_debug_layer(selected_mesh, mesh_fac,cameraPos, cameraFront, cameraUp);
        imgui_menu_bar(selected_mesh, mesh_fac);
        //imgui_N_Lateral_Parameters();
        if (is_trilateral_generated)
        {
            imgui_trilateralConfiguration(selected_mesh , mesh_fac);
        }
        
        imgui_render();*/

        /*glfwSwapBuffers(window);
        glfwPollEvents();*/
    }


   
    /*glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &IBO);

    glDeleteVertexArrays(1, &VAO_matching_points);
    glDeleteBuffers(1, &VBO_matching_points);
    glDeleteBuffers(1, &IBO_matchig_points);*/
    

    //imgui_close();

    rlImGuiShutdown();
    CloseWindow();
    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    return 0;
}

#include "Include/callback.h"
#include "Include/imgui_stuff.h"

