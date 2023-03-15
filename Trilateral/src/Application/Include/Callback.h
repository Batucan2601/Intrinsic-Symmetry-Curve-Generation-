#pragma once
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{


    if (yoffset > 0.0f)
    {
        cameraPos += cameraSpeed * cameraFront * 100000.0f;
    }
    else
    {
        cameraPos -= cameraSpeed * cameraFront * 100000.0f;

    }
}
void poll_keys()
{
        /* else if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS) {
             cameraPos = glm::vec3(0.0f, 0.0f, 5.0f);
             cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
             cameraUp = glm::vec3(1.0f, 0.0f, 0.0f);

             camera_mouse_lock = true;
         }*/

         // w a s d movement keys ( for data testing will be removed later)
        float movement_speed = cameraSpeed / 50;
        if (!is_e_pressed)
        {
            if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
                cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * movement_speed;
            else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
                cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * movement_speed;
            else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            {
                cameraPos += movement_speed * glm::normalize(cameraFront);
                //cameraPos.z += movement_speed * glm::normalize(cameraFront).z;
            }
            else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            {
                cameraPos -= movement_speed * glm::normalize(cameraFront);
            }
            // return to original location
            else if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
            {
                cameraPos = glm::vec3(0.0f, 0.0f, -5.0f); // where camera starts
                cameraFront = glm::vec3(0.0f, 0.0f, 0.0f);
                cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
            }
            else if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
            {
                 cameraPos = glm::vec3(20.0f, 1.25f, 0.0f); // where camera starts
                 cameraFront = glm::vec3(-1.0f, 0.0f, 0.0f);
                 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
            }
        }
        
        
  

}
bool firstMouse = true ; //holds first mouse movement
float lastX = 0, lastY = 0; // holds the global camera variables
float yaw = -90.0f;
float pitch = 45.0f;
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (!is_e_pressed)
    {
        if (firstMouse)
        {
            lastX = xpos;
            lastY = ypos;
            firstMouse = false;
        }

        float xoffset = xpos - lastX;
        float yoffset = lastY - ypos;
        lastX = xpos;
        lastY = ypos;

        float sensitivity = 0.1f;
        xoffset *= sensitivity;
        yoffset *= sensitivity;

        yaw += xoffset;
        pitch += yoffset;

        if (pitch > 89.0f)
            pitch = 89.0f;
        if (pitch < -89.0f)
            pitch = -89.0f;

        glm::vec3 direction;
        direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        direction.y = sin(glm::radians(pitch));
        direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraFront = glm::normalize(direction);
    }
}
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_E && action == GLFW_RELEASE)
    {
        is_e_pressed = !is_e_pressed;
        if (is_e_pressed)
        {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }
        else
        {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

        }
        std::cout << "e released" << std::endl;
        
    }
}