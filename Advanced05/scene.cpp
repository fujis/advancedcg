/*!
  @file rx_controller.cpp

  @brief GLUTによるOpenGLウィンドウクラス

  @author Makoto Fujisawa
  @date   2020-06
*/
// FILE --rx_controller.cpp--

#pragma warning (disable: 4996)
#pragma warning (disable: 4819)


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "scene.h"


//-----------------------------------------------------------------------------
// 定数・変数
//-----------------------------------------------------------------------------
const GLfloat RX_LIGHT0_POS[4] = { 2.0f, 4.0f, 2.0f, 0.0f };
const GLfloat RX_LIGHT1_POS[4] = { -1.0f, -10.0f, -1.0f, 0.0f };
const GLfloat RX_LIGHT_AMBI[4] = { 0.3f, 0.3f, 0.3f, 1.0f };
const GLfloat RX_LIGHT_DIFF[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat RX_LIGHT_SPEC[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat RX_FOV = 45.0f;


//-----------------------------------------------------------------------------
// SceneMLSクラスのstatic変数の定義と初期化
//-----------------------------------------------------------------------------
int SceneMLS::m_winw = 1000;					//!< 描画ウィンドウの幅
int SceneMLS::m_winh = 1000;					//!< 描画ウィンドウの高さ
double SceneMLS::m_bgcolor[3] = { 1, 1, 1 };	//!< 背景色
bool SceneMLS::m_animation_on = false;			//!< アニメーションON/OFF

int SceneMLS::m_draw = 0;						//!< 描画フラグ
int SceneMLS::m_currentstep = 0;				//!< 現在のステップ数
int SceneMLS::m_simg_spacing = -1;				//!< 画像保存間隔(=-1なら保存しない)

rxMeshDeform2D* SceneMLS::m_md = 0;				//!< メッシュ変形
glm::vec2 SceneMLS::m_envmin = glm::vec2(-1.4);	//!< メッシュ描画範囲(最小座標)
glm::vec2 SceneMLS::m_envmax = glm::vec2(1.4);	//!< メッシュ描画範囲(最大座標)
GLuint SceneMLS::m_tex = 0;						//!< テクスチャ
int SceneMLS::m_picked = -1;					//!< マウスピックされた頂点番号

rxGLSL SceneMLS::m_shader;						//!< GLSLシェーダ


/*!
 * 初期化関数
 */
void SceneMLS::Init(int argc, char* argv[])
{
	// GLEWの初期化
	GLenum err = glewInit();
	if(err != GLEW_OK) cout << "GLEW Error : " << glewGetErrorString(err) << endl;

	// マルチサンプル設定
	//GLint buf, sbuf;
	//glGetIntegerv(GL_SAMPLE_BUFFERS, &buf);
	//cout << "number of sample buffers is " << buf << endl;
	//glGetIntegerv(GL_SAMPLES, &sbuf);
	//cout << "number of samples is " << sbuf << endl;
	glEnable(GL_MULTISAMPLE);

	// 背景色の設定
	glClearColor((GLfloat)m_bgcolor[0], (GLfloat)m_bgcolor[1], (GLfloat)m_bgcolor[2], 1.0f);
	glClearDepth(1.0f);

	// 描画設定
	glEnable(GL_POINT_SMOOTH);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// 描画フラグ初期化
	m_draw = 0;
	m_draw |= RXD_EDGE;
	m_draw |= RXD_FACE;
	m_draw |= RXD_FIX;
	m_draw |= RXD_TEXTURE;

	// ファイルパス検索
	PathFinder p;
	p.addSearchPath(".");
	p.addSearchPath("bin");
	p.addSearchPath("../bin");
	p.addSearchPath("../../bin");

	// シェーダの初期化
	m_shader = CreateGLSLFromFile(p.find("shaders/shading.vp"), p.find("shaders/shading.fp"), "shading");

	// テクスチャ
	std::string filename = p.find("sample.bmp");
	//std::string filename = "sample.bmp";
	glActiveTexture(GL_TEXTURE0);
	if(!loadTexture(filename, m_tex, false, false)){
		cout << "Error : could not find the image file - " << filename << endl;
	}

	// メッシュファイル初期化
	m_md = new rxMeshDeform2D();
	m_md->InitVAO();

	// アニメーションON
	switchanimation(1);
}


/*!
 * 再描画イベント処理関数
 */
void SceneMLS::Draw(void)
{
	// プロジェクション変換行列の設定
	//glm::mat4 mp = glm::perspective(glm::radians(RX_FOV), (float)m_winw/m_winh, 0.03f, 1000.0f);
	glm::mat4 mp = glm::ortho(m_envmin[0], m_envmax[0], m_envmin[1], m_envmax[1], -1.0f, 1.0f);

	// 描画バッファのクリア
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_tex);

	glDisable(GL_LIGHTING);

	// GLSLシェーダをセット
	glUseProgram(m_shader.Prog);

	// uniform
	int loc_color = glGetUniformLocation(m_shader.Prog, "color");
	int loc_mproj = glGetUniformLocation(m_shader.Prog, "mat_proj");
	int loc_mv = glGetUniformLocation(m_shader.Prog, "mat_mv");
	int loc_tex = glGetUniformLocation(m_shader.Prog, "tex");
	int loc_alpha = glGetUniformLocation(m_shader.Prog, "alpha");

	glUniformMatrix4fv(loc_mproj, 1, GL_FALSE, glm::value_ptr(mp));	// 変換行列
	glUniform1i(loc_tex, 0);	// テクスチャユニット番号(glGenTextureで得られたテクスチャIDとは異なるので注意)

	// 面描画(テクスチャON/OFFあり)
	if(m_draw & RXD_FACE){
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glm::mat4 m = glm::translate(glm::vec3(0.0, 0.0, 0.0));
		glUniformMatrix4fv(loc_mv, 1, GL_FALSE, glm::value_ptr(m));	// 変換行列
		glUniform3f(loc_color, 0.0, 0.5, 1.0);	// 描画色
		glUniform1f(loc_alpha, ((m_draw&RXD_TEXTURE) ? 1.0f : 0.0f));
		if(m_md) m_md->DrawMesh();
	}

	// エッジ描画
	if(m_draw & RXD_EDGE){
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glm::mat4 m = glm::translate(glm::vec3(0.0, 0.0, 0.01));
		glUniformMatrix4fv(loc_mv, 1, GL_FALSE, glm::value_ptr(m));	// 変換行列
		glUniform3f(loc_color, 0.8f, 0.8f, 0.8f);	// 描画色
		glUniform1f(loc_alpha, 0.0f);
		if(m_md) m_md->DrawMesh();
	}

	// 頂点描画
	if(m_draw & RXD_VERTEX){
		glm::mat4 m = glm::translate(glm::vec3(0.0, 0.0, 0.02));
		glUniformMatrix4fv(loc_mv, 1, GL_FALSE, glm::value_ptr(m));	// 変換行列
		glUniform3f(loc_color, 1.0f, 1.0f, 0.0f);	// 描画色
		glUniform1f(loc_alpha, 0.0f);
		glPointSize(5.0);
		if(m_md) m_md->DrawPoints();
	}

	// 固定頂点描画
	if(m_draw & RXD_FIX){
		glm::mat4 m = glm::translate(glm::vec3(0.0, 0.0, 0.03));
		glUniformMatrix4fv(loc_mv, 1, GL_FALSE, glm::value_ptr(m));	// 変換行列
		glUniform3f(loc_color, 1.0f, 0.0f, 0.0f);	// 描画色
		glUniform1f(loc_alpha, 0.0f);
		glPointSize(24.0);
		if(m_md) m_md->DrawFixPoints();
	}

	glBindTexture(GL_TEXTURE_2D, 0);
}

/*!
 * タイマーコールバック関数
 */
void SceneMLS::Timer(void)
{
	if(m_animation_on){
		// 描画を画像ファイルとして保存
		if(m_simg_spacing > 0 && m_currentstep%m_simg_spacing == 0) savedisplay(-1);

		if(m_md){
			m_md->Update(0.01);
		}

		if(m_currentstep > MAX_STEPS) m_currentstep = 0;
		m_currentstep++;
	}
}


/*!
 * マウスイベント処理関数
 * @param[in] button マウスボタン(GLFW_MOUSE_BUTTON_LEFT,GLFW_MOUSE_BUTTON_MIDDLE,GLFW_MOUSE_BUTTON_RIGHT)
 * @param[in] action マウスボタンの状態(GLFW_PRESS, GLFW_RELEASE)
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void SceneMLS::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	glm::vec2 mpos(x/(double)m_winw, (m_winh-y-1.0)/(double)m_winh);
	mpos *= m_envmax-m_envmin;
	mpos += m_envmin;

	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS){
			int picked = -1;
			picked = m_md->SearchCP(mpos, 0.05);
			if(picked == -1) picked = m_md->Search(mpos, 0.05);
			if(picked != -1){
				m_picked = picked;
				m_md->SetCP(m_picked, mpos, false);
			}
			else if(m_picked){
				m_picked = -1;
			}
		}
		else if(action == GLFW_RELEASE){
			if(m_picked != -1){
				m_picked = -1;
			}
		}
	}
	else if(button == GLFW_MOUSE_BUTTON_RIGHT){
		if(action == GLFW_PRESS){
			int picked = -1;
			picked = m_md->SearchCP(mpos, 0.05);
			if(picked != -1){
				m_md->UnsetCP(picked);
				m_picked = -1;
			}
		}
	}
}
/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void SceneMLS::Cursor(GLFWwindow* window, double x, double y)
{
	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE){
		return;
	}

	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
		if(m_picked != -1){
			glm::vec2 mpos(x/(double)m_winw, (m_winh-y-1.0)/(double)m_winh);
			mpos *= m_envmax-m_envmin;
			mpos += m_envmin;
			m_md->SetCP(m_picked, mpos, true);
			//m_md->Update(0.01);
		}
	}
}

/*!
 * キーボードイベント処理関数
 * @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void SceneMLS::Keyboard(GLFWwindow* window, int key, int mods)
{
	switch(key){
	case GLFW_KEY_R: // リセット
		if(m_md) m_md->Init();
		break;

	case GLFW_KEY_X: // デバッグ用
		break;

	default:
		break;
	}
}


/*!
 * リサイズイベント処理関数
 * @param[in] w キャンバス幅(ピクセル数)
 * @param[in] h キャンバス高さ(ピクセル数)
 */
void SceneMLS::Resize(GLFWwindow* window, int w, int h)
{
	m_winw = w; m_winh = h;
	glViewport(0, 0, m_winw, m_winh);
}

/*!
 * ImGUIのウィジット配置
 */
void SceneMLS::ImGui(GLFWwindow* window)
{
	ImGui::Text("menu:");

	if(ImGui::Button("start/stop")){ switchanimation(-1); }
	if(ImGui::Button("save screenshot")){ savedisplay(-1); }
	ImGui::Separator();
	ImGui::CheckboxFlags("vertex", &m_draw, RXD_VERTEX);
	ImGui::CheckboxFlags("edge", &m_draw, RXD_EDGE);
	ImGui::CheckboxFlags("face", &m_draw, RXD_FACE);
	ImGui::CheckboxFlags("fix points", &m_draw, RXD_FIX);
	ImGui::CheckboxFlags("texture", &m_draw, RXD_TEXTURE);
	ImGui::Separator();
	if(ImGui::Button("reset(grid)")){ if(m_md) m_md->Init(0); }
	if(ImGui::Button("reset(random)")){ if(m_md) m_md->Init(1); }
	ImGui::Separator();
	ImGui::InputDouble("alpha", &(m_md->m_alpha), 0.05, 3.0, "%.2f");
	if(ImGui::RadioButton("affine", m_md->m_deform_type == 0)){ m_md->m_deform_type = 0; }
	if(ImGui::RadioButton("similarity", m_md->m_deform_type == 1)){ m_md->m_deform_type = 1; }
	if(ImGui::RadioButton("rigid", m_md->m_deform_type == 2)){ m_md->m_deform_type = 2; }
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

/*!
 * 終了処理関数
 */
void SceneMLS::Destroy()
{
	if(m_md) delete m_md;
}




/*!
 * アニメーションN/OFF
 * @param[in] on trueでON, falseでOFF
 */
bool SceneMLS::switchanimation(int on)
{
	m_animation_on = (on == -1) ? !m_animation_on : (on ? true : false);
	return m_animation_on;
}

/*!
 * 現在の画面描画を画像ファイルとして保存(連番)
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 */
void SceneMLS::savedisplay(const int &stp)
{
	static int nsave = 1;
	string fn = CreateFileName("img_", ".bmp", (stp == -1 ? nsave++ : stp), 5);
	saveFrameBuffer(fn, m_winw, m_winh);
	std::cout << "saved the screen image to " << fn << std::endl;
}



