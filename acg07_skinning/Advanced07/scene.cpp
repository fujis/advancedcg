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

// OpenGL描画関連
#include "rx_trackball.h"	// 視点変更用トラックボールクラス
#include "rx_shaders.h"		// GLSL関数



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
// SceneCAクラスのstatic変数の定義と初期化
//-----------------------------------------------------------------------------
int SceneCA::m_winw = 1000;					//!< 描画ウィンドウの幅
int SceneCA::m_winh = 1000;					//!< 描画ウィンドウの高さ
rxTrackball SceneCA::m_view;				//!< トラックボール
float SceneCA::m_bgcolor[3] = { 1, 1, 1 };	//!< 背景色
bool SceneCA::m_animation_on = false;		//!< アニメーションON/OFF

int SceneCA::m_draw = 0;					//!< 描画フラグ
int SceneCA::m_currentstep = 0;				//!< 現在のステップ数
int SceneCA::m_simg_spacing = -1;			//!< 画像保存間隔(=-1なら保存しない)

int SceneCA::m_picked = -1;
float SceneCA::m_pickdist = 1.0;	//!< ピックされた点までの距離

CharacterAnimation SceneCA::m_ca;
double SceneCA::m_scale = 1.0;
double SceneCA::m_yoffset = 0.0;

rxPolygonsE SceneCA::m_poly_org;
rxPolygonsE SceneCA::m_poly;
GLuint SceneCA::m_tex = 0;


/*!
 * 初期化
 */
void SceneCA::Init(int argc, char* argv[])
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
	m_draw |= RXD_BBOX;
	m_draw |= RXD_FLOOR;
	m_draw |= RXD_PARAMS;

	// キャラクターアニメーション初期設定
	initCA();
}

/*!
 * 再描画イベント処理関数
 */
void SceneCA::Draw(void)
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

	// 床描画
	glEnable(GL_LIGHTING);
	glm::vec3 lightpos(RX_LIGHT0_POS[0], RX_LIGHT0_POS[1], RX_LIGHT0_POS[2]);
	glm::vec3 lightcol(RX_LIGHT_DIFF[0], RX_LIGHT_DIFF[1], RX_LIGHT_DIFF[2]);
	if(m_draw & RXD_FLOOR) drawFloor(lightpos, lightcol, cen[1]-0.5*dim[1], 30.0);

	// 境界,軸描画
	glDisable(GL_LIGHTING);
	if(m_draw & RXD_BBOX) drawWireCuboid(cen-0.5f*dim, cen+0.5f*dim, glm::vec3(0.5, 1.0, 0.5), 2.0);
	if(m_draw & RXD_AXIS) drawAxis(0.6*dim[0], 3.0);

	// スケルトンの描画
	glEnable(GL_LIGHTING);
	glColor3d(0.0, 0.0, 1.0);
	glPushMatrix();
	glTranslated(0.0, m_yoffset, 0.0);

	m_ca.Draw(m_currentstep, m_scale);

	//for(int i = 0; i < m_ca.m_pos.size(); ++i){
	//	glm::vec3 p = m_ca.m_pos[i]*m_scale;
	//	glPushMatrix();
	//	glTranslatef(p[0], p[1], p[2]);
	//	glColor3f(1.0, 1.0, 0.0);
	//	glutSolidSphere(0.02, 20, 10);
	//	glPopMatrix();

	//}

	// メッシュ描画
	glScaled(m_scale, m_scale, m_scale);
	if(m_draw & RXD_MESH) m_poly.Draw(4);


	glPopMatrix();

	glPopMatrix();
}

/*!
 * タイマーコールバック関数
 */
void SceneCA::Timer(void)
{
	if(m_animation_on){
		// 描画を画像ファイルとして保存
		if(m_simg_spacing > 0 && m_currentstep%m_simg_spacing == 0) savedisplay(-1);

		m_poly.vertices = m_poly_org.vertices;
		m_ca.SkinningLBS(m_currentstep, m_poly.vertices, m_poly.weights);
		//m_ca.SkinningDQS(m_currentstep, m_poly.vertices, m_poly.weights);

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
void SceneCA::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	glm::vec2 mpos(x/(float)m_winw, (m_winh-y-1.0)/(float)m_winh);

	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS){
			// マウスドラッグによる視点移動
			m_view.Start(x, y, mods);
		}
		else if(action == GLFW_RELEASE){
			m_view.Stop(x, y);
		}
	}
}
/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void SceneCA::Cursor(GLFWwindow* window, double x, double y)
{
	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE && 
	   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
	   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE){
		return;
	}

	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
		m_view.Motion(x, y);
	}
}

/*!
 * キーボードイベント処理関数
 * @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void SceneCA::Keyboard(GLFWwindow* window, int key, int mods)
{
	switch(key){
	case GLFW_KEY_W: // weightのスムージング
		WeightFairingByUmbrella(m_poly);
		NormalizeWeights(m_poly);
		break;

	case GLFW_KEY_X: // デバッグ用
		DebugWeights(m_poly, "_weights.txt");
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
void SceneCA::Resize(GLFWwindow* window, int w, int h)
{
	m_winw = w; m_winh = h;
	m_view.SetRegion(w, h);
	glViewport(0, 0, m_winw, m_winh);
}

/*!
 * ImGUIのウィジット配置
 */
void SceneCA::ImGui(GLFWwindow* window)
{
	ImGui::Text("menu:");

	if(ImGui::Button("start/stop")){ switchanimation(-1); }
	if(ImGui::Button("run a step")){ Timer(); }
	if(ImGui::Button("reset viewpos")){ resetview(); }
	if(ImGui::Button("save screenshot")){ savedisplay(-1); }
	ImGui::Separator();
	ImGui::CheckboxFlags("aabb", &m_draw, RXD_BBOX);
	ImGui::CheckboxFlags("axis", &m_draw, RXD_AXIS);
	ImGui::CheckboxFlags("floor", &m_draw, RXD_FLOOR);
	ImGui::CheckboxFlags("params", &m_draw, RXD_PARAMS);
	ImGui::CheckboxFlags("mesh", &m_draw, RXD_MESH);
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

void SceneCA::Destroy()
{
}




/*!
 * アニメーションN/OFF
 * @param[in] on trueでON, falseでOFF
 */
bool SceneCA::switchanimation(int on)
{
	m_animation_on = (on == -1) ? !m_animation_on : (on ? true : false);
	return m_animation_on;
}

/*!
 * 現在の画面描画を画像ファイルとして保存(連番)
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 */
void SceneCA::savedisplay(const int &stp)
{
	static int nsave = 1;
	string fn = CreateFileName("img_", ".bmp", (stp == -1 ? nsave++ : stp), 5);
	saveFrameBuffer(fn, m_winw, m_winh);
	std::cout << "saved the screen image to " << fn << std::endl;
}

/*!
 * 視点の初期化
 */
void SceneCA::resetview(void)
{
	double q[4] = { 1, 0, 0, 0 };
	m_view.SetQuaternion(q);
	m_view.SetScaling(-5.0);
	m_view.SetTranslation(0.0, 0.0);
}





/*!
 * キャラクターアニメーションの初期化
 */
void SceneCA::initCA(void)
{
	// BVHファイル読み込み
	m_ca.Read("sample_arm2.bvh");
	//m_ca.Read("sample_walking1.bvh");

	// スケルトンモデルのAABB取得とスケール計算
	glm::vec3 minp, maxp;
	m_ca.AABB(minp, maxp);
	cout << "AABB of BVH : " << glm::to_string(minp) << " - " << glm::to_string(maxp) << endl;

	// BVHモデルに合わせて描画スケーリングとオフセットを設定
	glm::vec3 dim = maxp-minp;
	m_scale = 2.0/glm::max(dim[0], dim[1], dim[2]);
	m_yoffset = minp[1]*m_scale;
	cout << " scale : " << m_scale << endl;

	// メッシュモデル読み込み
	m_poly = rxPolygonsE();
	rxOBJ model;
	model.Read("simple_arm.obj", m_poly.vertices, m_poly.normals, m_poly.faces, m_poly.materials, true);
	if(!m_poly.vertices.empty()){
		cout << "[polygon]" << endl;
		cout << " num of vertices : " << m_poly.vertices.size() << endl;
		cout << " num of polygons : " << m_poly.faces.size() << endl;
		m_poly.open = 1;

		// メッシュモデルをBVHのAABBに合わせてスケーリング
		glm::vec3 poly_minp, poly_maxp;
		FindBBox(poly_minp, poly_maxp, m_poly.vertices); // 現在のBBoxの大きさを調べる

		glm::vec3 ctr = 0.5f*(maxp+minp), sl = 0.5f*(maxp-minp);
		glm::vec3 poly_ctr = 0.5f*(poly_maxp+poly_minp), poly_sl = 0.5f*(poly_maxp-poly_minp);
		glm::vec3 size_conv = 1.05f*sl/poly_sl;	// ボーンがはみ出ないように少しだけ大きくする

		// 全ての頂点をbboxにあわせて変換
		for(int i = 0; i < m_poly.vertices.size(); ++i){
			m_poly.vertices[i] = (m_poly.vertices[i]-poly_ctr)*size_conv+ctr;
		}

		// 元のメッシュモデルの縦横比を保存したままのスケーリングの場合こちらを使う
		//FitVertices(0.5*(maxp+minp), 0.5*(maxp-minp), m_poly.vertices);

		// Skinning weightの初期値
		m_poly.weights.resize(m_poly.vertices.size());
		m_ca.Weight(m_poly.vertices, m_poly.weights);

		// エッジデータの作成(weightのスムージングのため)
		SearchVertexFace(m_poly);
		SearchEdge(m_poly);

		// weightのスムージング
		for(int i = 0; i < 15; ++i) WeightFairingByUmbrella(m_poly);
		NormalizeWeights(m_poly);

		m_poly_org = m_poly;
	}
}


