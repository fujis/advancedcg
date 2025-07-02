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
// ScenePBDクラスのstatic変数の定義と初期化
//-----------------------------------------------------------------------------
int ScenePBD::m_winw = 1000;					//!< 描画ウィンドウの幅
int ScenePBD::m_winh = 1000;					//!< 描画ウィンドウの高さ
rxTrackball ScenePBD::m_view;					//!< トラックボール
float ScenePBD::m_bgcolor[3] = { 1, 1, 1 };	//!< 背景色
bool ScenePBD::m_animation_on = false;			//!< アニメーションON/OFF

int ScenePBD::m_draw = 0;						//!< 描画フラグ
int ScenePBD::m_currentstep = 0;				//!< 現在のステップ数
int ScenePBD::m_simg_spacing = -1;				//!< 画像保存間隔(=-1なら保存しない)

ElasticPBD* ScenePBD::m_elasticbody = 0;

int ScenePBD::m_picked = -1;
float ScenePBD::m_pickdist = 1.0;	//!< ピックされた点までの距離


void ScenePBD::Init(int argc, char* argv[])
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

	// PBD初期設定
	initRod();
	//initCloth();
	//initBall();
}


/*!
 * 再描画イベント処理関数
 */
void ScenePBD::Draw(void)
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

	// オブジェクト描画
	glEnable(GL_LIGHTING);
	glPushMatrix();
	if(m_elasticbody) m_elasticbody->Draw(m_draw);
	glPopMatrix();

	glPopMatrix();
}

/*!
 * タイマーコールバック関数
 */
void ScenePBD::Timer(void)
{
	if(m_animation_on){
		// 描画を画像ファイルとして保存
		if(m_simg_spacing > 0 && m_currentstep%m_simg_spacing == 0) savedisplay(-1);

		// シミュレーションを進める
		if(m_elasticbody) m_elasticbody->Update(DT);

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
void ScenePBD::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	glm::vec2 mpos(x/(float)m_winw, (m_winh-y-1.0)/(float)m_winh);

	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS){
			// マウスピック
			if(m_elasticbody && !(mods & 0x02)){
				// マウスクリックした方向のレイ
				glm::vec3 ray_from, ray_to;
				glm::vec3 init_pos = glm::vec3(0.0);
				m_view.CalLocalPos(ray_from, init_pos);
				m_view.GetRayTo(x, y, RX_FOV, ray_to);
				glm::vec3 dir = glm::normalize(ray_to-ray_from);	// 視点からマウス位置へのベクトル

				// レイと各頂点(球)の交差判定
				float d;
				int v = m_elasticbody->IntersectRay(ray_from, dir, d);
				if(v >= 0){
					//cout << "vertex " << v << " is selected - " << d << endl;
					m_picked = v;
					m_pickdist = d;
				}
			}

			if(m_picked == -1){
				// マウスドラッグによる視点移動
				m_view.Start(x, y, mods);
			}

		}
		else if(action == GLFW_RELEASE){
			m_view.Stop(x, y);
			if(m_picked != -1){
				if(mods & 0x01){	// SHIFTを押しながらボタンを離すことでその位置に固定
					m_elasticbody->FixVertex(m_picked);
				}
				else{
					m_elasticbody->UnFixVertex(m_picked);
					m_picked = -1;
				}
			}
		}
	}
}
/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void ScenePBD::Cursor(GLFWwindow* window, double x, double y)
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
			glm::vec3 ray_from, ray_to;
			glm::vec3 init_pos = glm::vec3(0.0);
			m_view.CalLocalPos(ray_from, init_pos);
			m_view.GetRayTo(x, y, RX_FOV, ray_to);
			glm::vec3 dir = glm::normalize(ray_to-ray_from);	// 視点からマウス位置へのベクトル

			glm::vec3 cur_pos = m_elasticbody->GetVertexPos(m_picked);
			glm::vec3 new_pos = ray_from+dir*static_cast<float>(m_pickdist);
			m_elasticbody->FixVertex(m_picked, new_pos);
		}
	}
}

/*!
 * キーボードイベント処理関数
 * @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void ScenePBD::Keyboard(GLFWwindow* window, int key, int mods)
{
	switch(key){
	case GLFW_KEY_S: // シミュレーションスタート/ストップ
		switchanimation(-1);
		break;

	case GLFW_KEY_C:	// 固定点の全解除
		m_elasticbody->UnFixAllVertex();
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
void ScenePBD::Resize(GLFWwindow* window, int w, int h)
{
	m_winw = w; m_winh = h;
	m_view.SetRegion(w, h);
	glViewport(0, 0, m_winw, m_winh);
}

/*!
 * ImGUIのウィジット配置
 */
void ScenePBD::ImGui(GLFWwindow* window)
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
	ImGui::CheckboxFlags("tetrahedral(wire)", &m_draw, RXD_TETRA);
	ImGui::CheckboxFlags("tetrahedral(cs)", &m_draw, RXD_TETCS);
	ImGui::CheckboxFlags("aabb", &m_draw, RXD_BBOX);
	ImGui::CheckboxFlags("axis", &m_draw, RXD_AXIS);
	ImGui::CheckboxFlags("floor", &m_draw, RXD_FLOOR);
	ImGui::Separator();
	if(ImGui::Button("rod (1D)")){ initRod(); }
	if(ImGui::Button("cloth (2D)")){ initCloth(); }
	if(ImGui::Button("ball (3D)")){ initBall(); }
	ImGui::Checkbox("use edges inside", &(m_elasticbody->m_bUseInEdge));
	ImGui::Separator();
	ImGui::InputInt("iterations", &(m_elasticbody->m_iNmax), 1, 20);
	ImGui::InputFloat("stiffness", &(m_elasticbody->m_fK), 0.01f, 1.0f, "%.2f");
	ImGui::InputFloat("wind", &(m_elasticbody->m_fWind), 0.01f, 0.5f, "%.2f");
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

void ScenePBD::Destroy()
{
	if(m_elasticbody) delete m_elasticbody;
}




/*!
 * アニメーションN/OFF
 * @param[in] on trueでON, falseでOFF
 */
bool ScenePBD::switchanimation(int on)
{
	m_animation_on = (on == -1) ? !m_animation_on : (on ? true : false);
	return m_animation_on;
}

/*!
 * 現在の画面描画を画像ファイルとして保存(連番)
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 */
void ScenePBD::savedisplay(const int &stp)
{
	static int nsave = 1;
	string fn = CreateFileName("img_", ".bmp", (stp == -1 ? nsave++ : stp), 5);
	saveFrameBuffer(fn, m_winw, m_winh);
	std::cout << "saved the screen image to " << fn << std::endl;
}

/*!
 * 視点の初期化
 */
void ScenePBD::resetview(void)
{
	double q[4] = { 1, 0, 0, 0 };
	m_view.SetQuaternion(q);
	m_view.SetScaling(-3.0);
	m_view.SetTranslation(0.0, 0.0);
}




/*!
 * 衝突処理関数
 * @param[in] p 現在の座標
 * @param[out] np 衝突後の座標
 * @param[in] v 速度
 * @param[in] obj オブジェクト番号
 */
void Collision(glm::vec3 &p, glm::vec3 &np, glm::vec3 &v, int obj)
{
	// 床以外の物体との衝突処理をしたい場合にここに書く

	//	// 球体との衝突判定
	//	glm::vec3 rpos = p-g_v3Cen;
	//	float d = norm(rpos)-g_fRad;
	//	if(d < 0.0){
	//		glm::vec3 n = Unit(rpos);
	//		np = g_v3Cen+n*g_fRad;
	//	}
}


/*!
 * ストランドで初期化
 */
void ScenePBD::initRod(void)
{
	m_draw |= RXD_VERTEX;
	m_draw |= RXD_EDGE;

	m_currentstep = 0;
	glm::vec3 env_min(-3.0, GOUND_HEIGHT, -3.0);
	glm::vec3 env_max(3.0, GOUND_HEIGHT+3.0, 3.0);

	if(m_elasticbody) delete m_elasticbody;
	m_elasticbody = new ElasticPBD(0);

	// PBDの設定
	m_elasticbody->SetSimulationSpace(env_min, env_max);
	m_elasticbody->SetCollisionFunc(Collision);

	// エッジ構造生成
	m_elasticbody->GenerateStrand(glm::vec3(0.0), glm::vec3(0.7, -0.7, 0.0), 30);

	// 固定点
	m_elasticbody->FixVertex(0);

	m_picked = -1;
}

/*!
 * 三角形メッシュで構成された布形状で初期化
 */
void ScenePBD::initCloth(void)
{
	m_draw |= RXD_EDGE;
	m_draw |= RXD_FACE;

	m_currentstep = 0;
	glm::vec3 env_min(-3.0, GOUND_HEIGHT, -3.0);
	glm::vec3 env_max(3.0, GOUND_HEIGHT+3.0, 3.0);

	if(m_elasticbody) delete m_elasticbody;
	m_elasticbody = new ElasticPBD(0);

	// PBDの設定
	m_elasticbody->SetSimulationSpace(env_min, env_max);
	m_elasticbody->SetCollisionFunc(Collision);

	// 三角形メッシュ生成
	int n = 24;	// 各軸の分割数:処理が思い場合はここを小さくする
	m_elasticbody->GenerateMesh(glm::vec2(-0.5, -0.5), glm::vec2(0.5, 0.5), n, n);

	// 固定点
	m_elasticbody->FixVertex(0);
	m_elasticbody->FixVertex(n*(n-1));

	m_picked = -1;
}


/*!
 * 四面体メッシュで構成されたボール形状で初期化
 */
void ScenePBD::initBall(void)
{
	m_draw &= ~RXD_VERTEX;
	m_draw |= RXD_EDGE;
	m_draw |= RXD_FACE;

	m_currentstep = 0;
	glm::vec3 env_min(-3.0, GOUND_HEIGHT, -3.0);
	glm::vec3 env_max(3.0, GOUND_HEIGHT+3.0, 3.0);

	if(m_elasticbody) delete m_elasticbody;
	m_elasticbody = new ElasticPBD(0);

	// PBDの設定
	m_elasticbody->SetSimulationSpace(env_min, env_max);
	m_elasticbody->SetCollisionFunc(Collision);

	// ファイルパス検索
	PathFinder p;
	p.addSearchPath(".");
	p.addSearchPath("bin");
	p.addSearchPath("../bin");
	p.addSearchPath("../../bin");

	// 四面体構造生成：処理が重い場合は"sphere_sparse.msh"の方を使う
	std::string filename = p.find("sphere.msh");
	//std::string filename = p.find("sphere_sparse.msh");
	m_elasticbody->GenerateTetrahedralMeshFromFile(filename);

	m_picked = -1;
}

/*!
 * マウスピックの解除
 */
void ScenePBD::clearPick(void)
{
	if(m_elasticbody && m_picked != -1) m_elasticbody->UnFixVertex(m_picked);
	m_picked = -1;
}
