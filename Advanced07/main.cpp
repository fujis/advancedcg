/*!
  @file main.cpp

  @brief メインファイル(GLUTコントローラクラスを呼び出す)

  @author Makoto Fujisawa
  @date   2020-06
*/


#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "glfw3dll.lib")

// コマンドプロンプトを出したくない場合はここをコメントアウト
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <iostream>
#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
//#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>

#include "scene.h"


#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif


using namespace std;

// scene entries
struct SceneEntry
{
	std::string name;
	void(*Init)(int argc, char* argv[]);
	void(*Draw)();
	void(*Timer)();
	void(*Cursor)(GLFWwindow* window, double xpos, double ypos);
	void(*Mouse)(GLFWwindow* window, int button, int action, int mods);
	void(*Keyboard)(GLFWwindow* window, int key, int mods);
	void(*Resize)(GLFWwindow* window, int w, int h);
	void(*ImGui)(GLFWwindow* window);
	void(*Destroy)();
};

#define SceneRegister(name, className) { name, (className::Init), (className::Draw), (className::Timer), (className::Cursor), (className::Mouse), \
										 (className::Keyboard), (className::Resize), (className::ImGui), (className::Destroy) }

SceneEntry g_SceneEntries[] = {
	SceneRegister("Scene 7: Character Animation", SceneCA)
};

#undef SceneRegister

const int g_NumSceneEntries = sizeof(g_SceneEntries) / sizeof(SceneEntry);
int g_SceneIndex = g_NumSceneEntries - 1;



static void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
	g_SceneEntries[g_SceneIndex].Cursor(window, xpos, ypos);
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	g_SceneEntries[g_SceneIndex].Mouse(window, button, action, mods);
}

void keyboardCallback(GLFWwindow* window, int key, int scanCode, int action, int mods)
{
	if(action == GLFW_PRESS){
		switch(key)
		{
		case GLFW_KEY_ESCAPE:
		case GLFW_KEY_Q:
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;
		default:
			g_SceneEntries[g_SceneIndex].Keyboard(window, key, mods);
		}
	}
}

void resizeCallback(GLFWwindow* window, int w, int h)
{
	g_SceneEntries[g_SceneIndex].Resize(window, w, h);
}


//-----------------------------------------------------------------------------
// メイン関数
//-----------------------------------------------------------------------------
/*!
 * メインルーチン
 * @param[in] argc コマンドライン引数の数
 * @param[in] argv コマンドライン引数
 */
int main(int argc, char *argv[])
{
	int winw = 1024, winh = 1024, winx = 100, winy = 100;

	if(!glfwInit()) return 1;
	glfwSetErrorCallback(glfw_error_callback);

#ifdef __APPLE__
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
#endif

	// Create window with graphics context
	GLFWwindow* window = glfwCreateWindow(winw, winh, "Advanced CG 07: Assignments", NULL, NULL);
	if(window == NULL) return 1;

	glfwMakeContextCurrent(window);
	glewExperimental = GL_TRUE;
	glfwSwapInterval(0); // Disable vsync

	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	cout << "Vendor: " << glGetString(GL_VENDOR) << endl;
	cout << "Renderer: " << glGetString(GL_RENDERER) << endl;

	for(int i = 0; i < g_NumSceneEntries; i++)
		g_SceneEntries[i].Init(argc, argv);

	// Setup callback functions
	glfwSetCursorPosCallback(window, cursorPosCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetKeyCallback(window, keyboardCallback);
	glfwSetFramebufferSizeCallback(window, resizeCallback);
	resizeCallback(window, winw, winh);

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsClassic();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL2_Init();
	//ImGui_ImplOpenGL3_Init();

	// Settings for timer
	float dt = DT;	// 1.0/FPS
	float cur_time = 0.0f, last_time = 0.0f, elapsed_time = 0.0f;
	glfwSetTime(0.0);	// Initialize the glfw timer

	// Main loop
	while(!glfwWindowShouldClose(window))
	{
		// Poll and handle events (inputs, window resize, etc.)
		glfwPollEvents();

		// OpenGL Rendering & Animation function
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		g_SceneEntries[g_SceneIndex].Draw();

		// Timer
		cur_time = glfwGetTime();
		elapsed_time = cur_time-last_time;
		if(elapsed_time >= dt){
			g_SceneEntries[g_SceneIndex].Timer();
			last_time = glfwGetTime();
		}

		// Start the ImGui frame
		ImGui_ImplOpenGL2_NewFrame();
		//ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// GUI
		ImGui::Begin("ImGui Window");
		ImGui::Text("Framerate: %.3f ms/frame (%.1f fps)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

		if(ImGui::BeginListBox("Scenes", ImVec2(0, ImGui::GetTextLineHeightWithSpacing() * g_NumSceneEntries * 1.05))){
			for(int i = 0; i < g_NumSceneEntries; i++){
				const bool is_selected = (g_SceneIndex == i);
				if(ImGui::Selectable(g_SceneEntries[i].name.c_str(), is_selected)) g_SceneIndex = i;
				if(is_selected) ImGui::SetItemDefaultFocus();
			}
			ImGui::EndListBox();
		}
		ImGui::Separator();
		g_SceneEntries[g_SceneIndex].ImGui(window);
		ImGui::End();

		// Rendering of the ImGUI frame in opengl canvas
		ImGui::Render();
		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
		//ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	// Cleanup
	for(int i = 0; i < g_NumSceneEntries; i++) g_SceneEntries[i].Destroy();
	ImGui_ImplOpenGL2_Shutdown();
	//ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();


	return 0;
}


