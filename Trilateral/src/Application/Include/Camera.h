#include <raylib.h>
#include <raylib/examples/shaders/rlights.h>
struct ModifiedCamera
{
	Camera camera; 
	Light light;
	Shader shader;
	void set_camera();
	void set_shader(Shader& shader);
	void set_light(); // do this after setting the shader pls
	void update();
};
