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

// テクスチャ・画像
#include "rx_texture.h"
#include "rx_bitmap.h"

// 3Dモデル
#include "rx_model_sa.h"

// マウスピック
#include "rx_pick.h"

// OpenGL描画関連
#include "rx_trackball.h"	// 視点変更用トラックボールクラス
#include "rx_gldraw.h"		// 各種OpenGL描画関数
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
// SceneMLSクラスのstatic変数の定義と初期化
//-----------------------------------------------------------------------------
int SceneMLS::m_winid = 0;						//!< GLUTのウインドウID
int SceneMLS::m_winw = 1000;					//!< 描画ウィンドウの幅
int SceneMLS::m_winh = 1000;					//!< 描画ウィンドウの高さ
int SceneMLS::m_mousebutton = 0;				//!< マウスボタンの状態
int SceneMLS::m_keymod;							//!< 修飾キーの状態
rxTrackball SceneMLS::m_view;					//!< トラックボール
double SceneMLS::m_bgcolor[3] = { 1, 1, 1 };	//!< 背景色
bool SceneMLS::m_animation_on = false;			//!< アニメーションON/OFF
bool SceneMLS::m_fullscreen_on = false;			//!< フルスクリーンON/OFF

int SceneMLS::m_currentstep = 0;				//!< 現在のステップ数
bool SceneMLS::m_pause = false;					//!< シミュレーションのポーズフラグ

int SceneMLS::m_draw = 0;						//!< 描画フラグ
		
int SceneMLS::m_simg_spacing = -1;				//!< 画像保存間隔(=-1なら保存しない)

rxWave* SceneMLS::m_wave = 0;
double SceneMLS::m_dt = 0.01;	//!< 時間ステップ幅

int SceneMLS::m_view_drag = 0;

rxPolygonsE SceneMLS::m_poly_org;
rxPolygonsE SceneMLS::m_poly;
GLuint SceneMLS::m_tex = 0;



void SceneMLS::Init(int argc, char* argv[])
{
	// OpenGLのバージョンとGLEWの初期化
	RXCOUT << "OpenGL Ver. " << glGetString(GL_VERSION) << endl;
	GLenum err = glewInit();
	if(err == GLEW_OK) RXCOUT << "GLEW OK : Glew Ver. " << glewGetString(GLEW_VERSION) << endl;
	else RXCOUT << "GLEW Error : " << glewGetErrorString(err) << endl;

	// マルチサンプル設定の確認
	GLint buf, sbuf;
	glGetIntegerv(GL_SAMPLE_BUFFERS, &buf);
	RXCOUT << "number of sample buffers is " << buf << endl;
	glGetIntegerv(GL_SAMPLES, &sbuf);
	RXCOUT << "number of samples is " << sbuf << endl;

	m_mousebutton = 0;

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glEnable(GL_MULTISAMPLE);

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

	// テクスチャ
	//loadTexture("sample.bmp", m_tex, false, false);

	// トラックボール初期姿勢
	m_view.SetTranslation(0, 0);
	m_view.SetScaling(-5);

	// 描画フラグ初期化
	m_draw = 0;
	m_draw |= RXD_EDGE;
	//m_draw |= RXD_BBOX;

	// 波クラス生成
	int n = 64;
	double scale = 4.0;
	double water_depth = 0.1;
	m_wave = new rxWave(n, scale, water_depth);

}


/*!
 * 再描画イベント処理関数
 */
void SceneMLS::Draw(void)
{
	// ビューポート,透視変換行列,モデルビュー変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(RX_FOV, (float)m_winw/(float)m_winh, 0.2f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// 描画バッファのクリア
	glClearColor((GLfloat)m_bgcolor[0], (GLfloat)m_bgcolor[1], (GLfloat)m_bgcolor[2], 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

	// マウスによる回転・平行移動の適用
	m_view.Apply();

	// 各種シミュレーションデータの取得
	Vec3 cen(0.0), dim(2.0);			// シミュレーション空間情報

	// 床描画
	glEnable(GL_LIGHTING);
	if(m_draw & RXD_FLOOR) drawFloor(Vec3f(RX_LIGHT0_POS), Vec3f(RX_LIGHT_DIFF), cen[1]-0.5*dim[1], 30.0);

	// 境界,軸描画
	glDisable(GL_LIGHTING);
	if(m_draw & RXD_BBOX) drawWireCuboid(cen-0.5*dim, cen+0.5*dim, Vec3(0.5, 1.0, 0.5), 2.0);
	if(m_draw & RXD_AXIS) drawAxis(0.6*dim[0], 3.0);

	// オブジェクト描画
	glEnable(GL_LIGHTING);
	glPushMatrix();

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glColor3d(0.0, 0.5, 1.0);
	if(m_wave) m_wave->Draw(m_draw);


	// メッシュ描画
	//if(m_draw & RXD_MESH) m_poly.Draw(4);


	glPopMatrix();

	glPopMatrix();
}

/*!
 * タイマーコールバック関数
 */
void SceneMLS::Timer(void)
{
	if(m_animation_on){
		// 描画を画像ファイルとして保存
		if(m_simg_spacing > 0 && m_currentstep%m_simg_spacing == 0) savedisplay(-1);

		// シミュレーションステップを進める
		if(m_wave) m_wave->Update(m_dt);

		if(m_currentstep > RX_MAX_STEPS) m_currentstep = 0;
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
	Vec2 mpos(x/(double)m_winw, (m_winh-y-1.0)/(double)m_winh);

	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS){
			// マウスピック
			rxGLPick pick;
			pick.Set(SceneMLS::RenderSceneForPick, 0, SceneMLS::ProjectionForPick, 0);
			int idx = pick.Pick(x, y);
			if(m_wave) m_wave->AddHeight(idx, 0.15);
			if(idx == -1){
				// マウスドラッグによる視点移動
				m_view.Start(x, y, mods);
				m_view_drag = 1;
			}
		}
		else if(action == GLFW_RELEASE){
			m_view.Stop(x, y);
			m_view_drag = 0;
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
		if(m_view_drag){
			m_view.Motion(x, y);
		}
		else{
			rxGLPick pick;
			pick.Set(SceneMLS::RenderSceneForPick, 0, SceneMLS::ProjectionForPick, 0);
			int idx = pick.Pick(x, y);
			if(m_wave) m_wave->AddHeight(idx, 0.15);
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
void SceneMLS::Resize(GLFWwindow* window, int w, int h)
{
	m_winw = w; m_winh = h;
	m_view.SetRegion(w, h);
	glViewport(0, 0, m_winw, m_winh);
}

/*!
 * ImGUIのウィジット配置
 */
void SceneMLS::ImGui(GLFWwindow* window)
{
	ImGui::Text("menu:");

	if(ImGui::Button("start/stop")){ switchanimation(-1); }
	if(ImGui::Button("run a step")){ Timer(); }
	if(ImGui::Button("reset viewpos")){ resetview(); }
	if(ImGui::Button("save screenshot")){ savedisplay(-1); }
	if(ImGui::Button("save mesh")){
		if(m_wave){
			rxPolygons &p = m_wave->GetMesh();
			RxModel::SaveOBJ("wave.obj", p);
		}
	}
	ImGui::Separator();
	ImGui::CheckboxFlags("vertex", &m_draw, RXD_VERTEX);
	ImGui::CheckboxFlags("edge", &m_draw, RXD_EDGE);
	ImGui::CheckboxFlags("face", &m_draw, RXD_FACE);
	ImGui::CheckboxFlags("aabb", &m_draw, RXD_BBOX);
	ImGui::CheckboxFlags("axis", &m_draw, RXD_AXIS);
	ImGui::CheckboxFlags("floor", &m_draw, RXD_FLOOR);
	ImGui::Separator();
	if(ImGui::Button("reset simulation")){ if(m_wave) m_wave->Init(); }
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

void SceneMLS::Destroy()
{
	if(m_wave) delete m_wave;
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
	string fn = CreateFileName(RX_DEFAULT_IMAGE_DIR+"pbf_", ".bmp", (stp == -1 ? nsave++ : stp), 5);
	saveFrameBuffer(fn, m_winw, m_winh);
	std::cout << "saved the screen image to " << fn << std::endl;
}

/*!
 * 視点の初期化
 */
void SceneMLS::resetview(void)
{
	double q[4] = { 1, 0, 0, 0 };
	m_view.SetQuaternion(q);
	m_view.SetScaling(-5.0);
	m_view.SetTranslation(0.0, 0.0);
}





/*!
 * 画面出力用の文字列の生成
 * @return 文字列
 */
vector<string> SceneMLS::setDrawString(void)
{
	vector<string> str;

	// 計測された計算時間
	str.push_back("time : ");
	string tstr;
	RXTIMER_STRING(tstr);
	divideString(tstr, str);

	return str;
}

/*!
 * マウスピック用の描画関数
 */
void SceneMLS::RenderSceneForPick(void* x)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	m_view.Apply();	// マウスによる回転・平行移動の適用
	if(m_wave) m_wave->DrawForPick();
	glPopMatrix();
}
void SceneMLS::ProjectionForPick(void* x)
{
	gluPerspective(RX_FOV, (float)m_winw/(float)m_winh, 0.2f, 1000.0f);
}

