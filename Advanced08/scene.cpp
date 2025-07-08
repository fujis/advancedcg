/*!
  @file scene.cpp

  @brief GLFWによるOpenGL描画

  @author Makoto Fujisawa
  @date   2021-04
*/


#pragma warning (disable: 4996)
#pragma warning (disable: 4819)


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "scene.h"

// テクスチャ・画像
#include "rx_texture.h"
#include "rx_bitmap.h"

// OpenGL描画関連
#include "rx_trackball.h"	// 視点変更用トラックボールクラス
#include "rx_shaders.h"		// GLSL関数



//-----------------------------------------------------------------------------
// 定数・変数
//-----------------------------------------------------------------------------
const GLfloat RX_LIGHT0_POS[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
const GLfloat RX_LIGHT1_POS[4] = { -1.0f, -10.0f, -1.0f, 0.0f };
const GLfloat RX_LIGHT_AMBI[4] = { 0.3f, 0.3f, 0.3f, 1.0f };
const GLfloat RX_LIGHT_DIFF[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat RX_LIGHT_SPEC[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat RX_FOV = 45.0f;


//-----------------------------------------------------------------------------
// SceneSWEクラスのstatic変数の定義と初期化
//-----------------------------------------------------------------------------
int SceneSWE::m_winw = 1000;					//!< 描画ウィンドウの幅
int SceneSWE::m_winh = 1000;					//!< 描画ウィンドウの高さ
rxTrackball SceneSWE::m_view;					//!< トラックボール
float SceneSWE::m_bgcolor[3] = { 1, 1, 1 };	//!< 背景色
bool SceneSWE::m_animation_on = false;			//!< アニメーションON/OFF

int SceneSWE::m_draw = 0;						//!< 描画フラグ
int SceneSWE::m_currentstep = 0;				//!< 現在のステップ数
int SceneSWE::m_simg_spacing = -1;				//!< 画像保存間隔(=-1なら保存しない)

int SceneSWE::m_picked = -1;
float SceneSWE::m_pickdist = 1.0;	//!< ピックされた点までの距離

WaveSWE* SceneSWE::m_wave = 0;
float SceneSWE::m_scale = 4.0;
int SceneSWE::m_res = 64;



void SceneSWE::Init(int argc, char* argv[])
{
	// GLEWの初期化
	GLenum err = glewInit();
	if(err != GLEW_OK) cout << "GLEW Error : " << glewGetErrorString(err) << endl;

	// マルチサンプル設定の確認
	//GLint buf, sbuf;
	//glGetIntegerv(GL_SAMPLE_BUFFERS, &buf);
	//cout << "number of sample buffers is " << buf << endl;
	//glGetIntegerv(GL_SAMPLES, &sbuf);
	//cout << "number of samples is " << sbuf << endl;
	glEnable(GL_MULTISAMPLE);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glEnable(GL_POINT_SMOOTH);

	// 光源&材質描画設定
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, RX_LIGHT0_POS);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, RX_LIGHT_DIFF);
	glLightfv(GL_LIGHT0, GL_SPECULAR, RX_LIGHT_SPEC);
	glLightfv(GL_LIGHT0, GL_AMBIENT, RX_LIGHT_AMBI);

	glShadeModel(GL_SMOOTH);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// 視点初期化
	resetview();

	// 描画フラグ初期化
	m_draw = 0;
	m_draw |= RXD_EDGE;
	m_draw |= RXD_FACE;
	m_draw |= RXD_BBOX;

	// 波の初期化
	initWave();

	// アニメーションON
	switchanimation(1);
}


/*!
 * 再描画イベント処理関数
 */
void SceneSWE::Draw(void)
{
	// ビューポート,透視変換行列,モデルビュー変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glm::mat4 mp = glm::perspective(RX_FOV, (float)m_winw/m_winh, 0.2f, 1000.0f);
	glMultMatrixf(glm::value_ptr(mp));
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// 描画バッファのクリア
	glClearColor((GLfloat)m_bgcolor[0], (GLfloat)m_bgcolor[1], (GLfloat)m_bgcolor[2], 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

	// マウスによる回転・平行移動の適用
	m_view.Apply();

	// シミュレーション空間情報
	glm::vec3 cen(0.0), dim(2.0);
	if(m_wave) m_wave->GetDim(cen, dim);

	// 床描画
	glEnable(GL_LIGHTING);
	glm::vec3 lightpos(RX_LIGHT0_POS[0], RX_LIGHT0_POS[1], RX_LIGHT0_POS[2]);
	glm::vec3 lightcol(RX_LIGHT_DIFF[0], RX_LIGHT_DIFF[1], RX_LIGHT_DIFF[2]);
	if(m_draw & RXD_FLOOR) drawFloor(lightpos, lightcol, cen[1]-0.5*dim[1], 30.0);

	// 境界,軸描画
	glDisable(GL_LIGHTING);
	if(m_draw & RXD_BBOX) drawWireCuboid(cen-0.5f*dim, cen+0.5f*dim, glm::vec3(0.5, 1.0, 0.5), 2.0);
	if(m_draw & RXD_AXIS) drawAxis(0.6*dim[0], 3.0);

	// オブジェクト描画
	glEnable(GL_LIGHTING);
	glPushMatrix();

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glColor3d(0.0, 0.5, 1.0);
	if(m_wave) m_wave->Draw(m_draw);

	glPopMatrix();

	glPopMatrix();
}

/*!
 * タイマーコールバック関数
 */
void SceneSWE::Timer(void)
{
	if(m_animation_on){
		// 描画を画像ファイルとして保存
		if(m_simg_spacing > 0 && m_currentstep%m_simg_spacing == 0) savedisplay(-1);

		// シミュレーションステップを進める
		if(m_wave) m_wave->Update(m_currentstep, g_dt);

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
void SceneSWE::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	glm::vec2 mpos(x/(float)m_winw, (m_winh-y-1.0)/(float)m_winh);

	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS){
			// マウスピック
			if(m_wave && !(mods & 0x02)){
				// マウスクリックした方向のレイ
				glm::vec3 ray_from, ray_to;
				glm::vec3 init_pos = glm::vec3(0.0);
				m_view.CalLocalPos(ray_from, init_pos);
				m_view.GetRayTo(x, y, RX_FOV, ray_to);
				glm::vec3 dir = glm::normalize(ray_to-ray_from);	// 視点からマウス位置へのベクトル

				// レイと各頂点(球)の交差判定
				float d;
				int v = m_wave->IntersectRay(ray_from, dir, d);
				if(v >= 0){
					//cout << "vertex " << v << " is selected - " << d << endl;
					m_picked = v;
					m_pickdist = d;
					m_wave->AddHeightArea(m_picked, 0.2);
				}
			}

			if(m_picked == -1){
				// マウスドラッグによる視点移動
				m_view.Start(x, y, mods);
			}

		}
		else if(action == GLFW_RELEASE){
			m_view.Stop(x, y);
			m_picked = -1;
		}
	}
}
/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void SceneSWE::Cursor(GLFWwindow* window, double x, double y)
{
	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE && 
	   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
	   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE){
		return;
	}

	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
		if(m_picked == -1){
			m_view.Motion(x, y);
		}
		else{
			if(m_wave){
				glm::vec3 ray_from, ray_to;
				glm::vec3 init_pos = glm::vec3(0.0);
				m_view.CalLocalPos(ray_from, init_pos);
				m_view.GetRayTo(x, y, RX_FOV, ray_to);
				glm::vec3 dir = glm::normalize(ray_to-ray_from);	// 視点からマウス位置へのベクトル

				float d;
				int v = m_wave->IntersectRay(ray_from, dir, d);
				if(v >= 0){
					m_picked = v;
					m_pickdist = d;
					m_wave->AddHeightArea(m_picked, 0.2);
				}
			}
		}
	}
}

/*!
 * キーボードイベント処理関数
 * @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void SceneSWE::Keyboard(GLFWwindow* window, int key, int mods)
{
	switch(key){
	case GLFW_KEY_S: // デバッグ用
		switchanimation(-1);
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
void SceneSWE::Resize(GLFWwindow* window, int w, int h)
{
	m_winw = w; m_winh = h;
	m_view.SetRegion(w, h);
	glViewport(0, 0, m_winw, m_winh);
}

/*!
 * ImGUIのウィジット配置
 */
void SceneSWE::ImGui(GLFWwindow* window)
{
	ImGui::Text("menu:");

	if(ImGui::Button("start/stop")){ switchanimation(-1); }
	if(ImGui::Button("run a step")){ m_animation_on = true; Timer(); m_animation_on = false; }
	if(ImGui::Button("reset viewpos")){ resetview(); }
	if(ImGui::Button("save screenshot")){ savedisplay(-1); }
	ImGui::Separator();
	ImGui::CheckboxFlags("vertex", &m_draw, RXD_VERTEX);
	ImGui::CheckboxFlags("edge", &m_draw, RXD_EDGE);
	ImGui::CheckboxFlags("face", &m_draw, RXD_FACE);
	ImGui::CheckboxFlags("aabb", &m_draw, RXD_BBOX);
	ImGui::CheckboxFlags("axis", &m_draw, RXD_AXIS);
	ImGui::CheckboxFlags("floor", &m_draw, RXD_FLOOR);
	ImGui::Separator();
	if(ImGui::Button("flat")){ initWave(); }
	if(ImGui::Button("slope")){ initWaveSlope(); }
	if(ImGui::Button("mountain")){ initWaveMountain(); }
	if(m_wave){
		ImGui::Separator();
		ImGui::Checkbox("sine-wave", &(m_wave->m_surf));
		ImGui::InputFloat("viscosity", &(m_wave->m_nu), 0.001f, 0.01f, "%.4f");
	}
	ImGui::Separator();
	ImGui::Text("scene parameters:");
	ImGui::InputInt("resolution", &(m_res), 4, 16);
	ImGui::InputFloat("scale", &(m_scale), 0.1f, 1.0f, "%.1f");
	ImGui::InputFloat("dt", &(g_dt), 0.001f, 0.01f, "%.3f");
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

void SceneSWE::Destroy()
{
	if(m_wave) delete m_wave;
}


/*!
 * アニメーションN/OFF
 * @param[in] on trueでON, falseでOFF
 */
bool SceneSWE::switchanimation(int on)
{
	m_animation_on = (on == -1) ? !m_animation_on : (on ? true : false);
	return m_animation_on;
}

/*!
 * 現在の画面描画を画像ファイルとして保存(連番)
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 */
void SceneSWE::savedisplay(const int &stp)
{
	static int nsave = 1;
	string fn = CreateFileName("img_", ".bmp", (stp == -1 ? nsave++ : stp), 5);
	saveFrameBuffer(fn, m_winw, m_winh);
	std::cout << "saved the screen image to " << fn << std::endl;
}

/*!
 * 視点の初期化
 */
void SceneSWE::resetview(void)
{
	double q[4] = { 1, 0, 0, 0 };
	m_view.SetQuaternion(q);
	m_view.SetRotation(20.0, 1.0, 0.0, 0.0);
	m_view.SetScaling(-5.0);
	m_view.SetTranslation(0.0, 0.0);
}



//! 底面の高さを返す関数
//  - x,yの範囲は(0,0)-(m_scale,m_scale)
float SceneSWE::getHeightFlat(float x, float y){ return -0.8; }
float SceneSWE::getHeightSlope(float x, float y){ return (0.05*x-0.1)*m_scale-0.8; }
float SceneSWE::getHeightMountain(float x, float y)
{
	x /= m_scale; y /= m_scale;
	return 0.5*exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/(2.0*0.1*0.1))-1.0;
}

/*!
 * 底面がフラットな波のシーン生成
 */
void SceneSWE::initWave(void)
{
	// 波クラス生成
	if(m_wave) delete m_wave;
	m_wave = new WaveSWE;

	m_wave->Init(m_res, m_scale, getHeightFlat);
}

/*!
 * 底面がフラットな波のシーン生成
 */
void SceneSWE::initWaveSlope(void)
{
	// 波クラス生成
	if(m_wave) delete m_wave;
	m_wave = new WaveSWE;

	m_wave->Init(m_res, m_scale, getHeightSlope);

}

void SceneSWE::initWaveMountain(void)
{
	// 波クラス生成
	if(m_wave) delete m_wave;
	m_wave = new WaveSWE;

	m_wave->Init(m_res, m_scale, getHeightMountain);
}


