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

CharacterAnimation* SceneCA::m_ca = 0;
float SceneCA::m_scale = 1.0;
glm::vec3 SceneCA::m_offset = glm::vec3(0.0);

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
	m_draw |= RXD_FACE;
	m_draw |= RXD_BBOX;
	m_draw |= RXD_FLOOR;
	m_draw |= RXD_SKELETON;

	// キャラクターアニメーション初期設定
	initCA("simple_arm.bvh", "simple_arm.obj");
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

	glEnable(GL_LIGHTING);
	glColor3d(0.0, 0.0, 1.0);
	glPushMatrix();
	glTranslatef(m_offset[0], m_offset[1], m_offset[2]);

	// スケルトンの描画
	if(m_ca && m_draw & RXD_SKELETON) m_ca->Draw(m_currentstep, m_scale);

	// メッシュ描画
	glScaled(m_scale, m_scale, m_scale);
	if(m_draw & RXD_MESH) drawPolygon(m_poly, m_draw);

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
		if(m_ca){
			m_ca->Skinning(m_currentstep, m_poly.vertices, m_poly.weights);

			// 頂点法線の更新
			CalVertexNormals(m_poly.vertices, m_poly.vertices.size(), m_poly.faces, m_poly.faces.size(), m_poly.normals);
		}

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
	case GLFW_KEY_S: // シミュレーションスタート/ストップ
		switchanimation(-1);
		break;

	case GLFW_KEY_W: // weightのスムージング
		WeightFairingByUmbrella(m_poly);
		NormalizeWeights(m_poly);
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
	ImGui::CheckboxFlags("skeleton", &m_draw, RXD_SKELETON);
	ImGui::CheckboxFlags("mesh", &m_draw, RXD_MESH);
	ImGui::CheckboxFlags("weights", &m_draw, RXD_WEIGHT);
	if(m_ca){
		ImGui::Checkbox("rest pose", &(m_ca->m_rest_pose));
		ImGui::Separator();
		if(ImGui::RadioButton("LBS", m_ca->m_skinning == 0)){
			m_ca->m_skinning = 0;
			m_poly.vertices = m_poly_org.vertices;
			m_ca->Skinning(m_currentstep, m_poly.vertices, m_poly.weights);
			CalVertexNormals(m_poly.vertices, m_poly.vertices.size(), m_poly.faces, m_poly.faces.size(), m_poly.normals);
		}
		if(ImGui::RadioButton("DQS", m_ca->m_skinning == 1)){
			m_ca->m_skinning = 1;
			m_poly.vertices = m_poly_org.vertices;
			m_ca->Skinning(m_currentstep, m_poly.vertices, m_poly.weights);
			CalVertexNormals(m_poly.vertices, m_poly.vertices.size(), m_poly.faces, m_poly.faces.size(), m_poly.normals);
		}
		ImGui::Separator();
		if(ImGui::Button("arm")){ initCA("simple_arm.bvh", "simple_arm.obj"); }
		if(ImGui::Button("arm (twist)")){ initCA("simple_arm_twist.bvh", "simple_arm.obj"); }
		if(ImGui::Button("walking")){ initCA("simple_walking.bvh", "simple_human.obj"); }
	}
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GL_TRUE); }

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
void SceneCA::initCA(const string &fn_bvh, const string &fn_obj)
{
	if(m_ca) delete m_ca;
	m_ca = new CharacterAnimation;

	// ファイルパス検索
	PathFinder p;
	p.addSearchPath(".");
	p.addSearchPath("bin");
	p.addSearchPath("../bin");
	p.addSearchPath("../../bin");

	// BVHファイル読み込み
	if(!m_ca->Read(p.find(fn_bvh))){
		return;
	}
	//m_ca->CheckData(); // デバッグ用

	// スケルトンモデルのAABB取得とスケール計算
	glm::vec3 minp, maxp;
	m_ca->AABB(minp, maxp, false);
	cout << "AABB of BVH : " << glm::to_string(minp) << " - " << glm::to_string(maxp) << endl;

	// BVHモデルに合わせて描画スケーリングとオフセットを設定
	glm::vec3 dim = maxp-minp;
	m_scale = 2.0/glm::max(dim[0], dim[1], dim[2]);
	m_offset = -0.5f*(minp+maxp)*m_scale;

	// メッシュモデル読み込み
	m_poly = rxPolygonsE();
	rxOBJ model;
	model.Read(p.find(fn_obj), m_poly.vertices, m_poly.normals, m_poly.faces, m_poly.materials, false);

	if(!m_poly.vertices.empty()){
		cout << "[polygon]" << endl;
		cout << " num of vertices : " << m_poly.vertices.size() << endl;
		cout << " num of polygons : " << m_poly.faces.size() << endl;
		m_poly.open = 1;

		// メッシュモデル用にスケルトンの初期ポーズ(キャラクターモデルならT-Pose)でのAABBを取得
		m_ca->AABB(minp, maxp, true);

		// メッシュモデルをBVHのAABBに合わせてスケーリング
		glm::vec3 poly_minp, poly_maxp;
		FindBBox(poly_minp, poly_maxp, m_poly.vertices); // 現在のBBoxの大きさを調べる

		glm::vec3 ctr = 0.5f*(maxp+minp), sl = 0.5f*(maxp-minp);
		glm::vec3 poly_ctr = 0.5f*(poly_maxp+poly_minp), poly_sl = 0.5f*(poly_maxp-poly_minp);
		glm::vec3 size_conv = 1.05f*sl/poly_sl;	// ボーンがはみ出ないように少しだけ大きくする

		// 元のメッシュモデルの縦横比を保存したままのスケーリング
		FitVertices(0.5f*(maxp+minp), 0.5f*(maxp-minp), m_poly.vertices);

		// Skinning weightの初期値
		m_poly.weights.resize(m_poly.vertices.size());
		m_ca->Weight(m_poly.vertices, m_poly.weights);

		// エッジデータの作成(weightのスムージングのため)
		SearchVertexFace(m_poly);
		SearchEdge(m_poly);
		cout << " num of edges : " << m_poly.edges.size() << endl;

		// 頂点法線の計算
		CalVertexNormals(m_poly.vertices, m_poly.vertices.size(), m_poly.faces, m_poly.faces.size(), m_poly.normals);

		// weightのスムージング
		for(int i = 0; i < 100; ++i) WeightFairingByUmbrella(m_poly);
		NormalizeWeights(m_poly);

		cout << " finish." << endl;

		m_poly_org = m_poly;
	}

	m_currentstep = 0;
}



/*!
 * ポリゴンの描画
 * @param[in] polys ポリゴンデータ
 * @param[in] draw 描画フラグ(下位ビットから頂点,エッジ,面,法線 - 1,2,4,8)
 */
void SceneCA::drawPolygon(rxPolygonsE &poly, int draw, float dn, bool col)
{
	// 頂点数とポリゴン数
	int vn = (int)poly.vertices.size();
	int pn = (int)poly.faces.size();

	if(draw & RXD_EDGE){
		// エッジ描画における"stitching"をなくすためのオフセットの設定
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}

	if(draw & RXD_FACE){
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_LIGHTING);
		
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnable(GL_COLOR_MATERIAL);
		glColor3d(0.2, 0.2, 0.7);

		bool tc = false;
		rxMaterialOBJ *cur_mat = NULL, *pre_mat = NULL;

		glm::vec3 cols[6] = { glm::vec3(1.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 1.0, 0.0),
							  glm::vec3(0.0, 1.0, 1.0), glm::vec3(1.0, 1.0, 0.0), glm::vec3(1.0, 0.0, 1.0) };

		// すべてのポリゴンを描画
		for(int i = 0; i < pn; ++i){
			const rxFace *face = &(poly.faces[i]);
			int n = (int)face->vert_idx.size();

			// ポリゴン描画
			glBegin(GL_POLYGON);
			for(int j = 0; j < n; ++j){
				int idx = face->vert_idx[j];
				if(idx >= 0 && idx < vn){
					glm::vec3 col(0.0);
					if(draw & RXD_WEIGHT && !poly.weights.empty() && !poly.weights[idx].empty()){
						int ncol = 0;
						map<int, double>::iterator iter;
						for(iter = poly.weights[idx].begin(); iter != poly.weights[idx].end(); ++iter){
							int joint_idx = iter->first;
							col += float(iter->second)*cols[joint_idx%6];
							ncol++;
						}
						glColor3fv(glm::value_ptr(col));
					}

					glNormal3fv(glm::value_ptr(poly.normals[idx]));
					if(!face->texcoords.empty()) glTexCoord2fv(glm::value_ptr(face->texcoords[j]));
					glVertex3fv(glm::value_ptr(poly.vertices[idx]));
				}
			}
			glEnd();
		}

		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_TEXTURE_2D);
	}


	// 頂点描画
	if(draw & RXD_VERTEX){
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glDisable(GL_LIGHTING);
		glPointSize(5.0);
		glColor3d(1.0, 0.3, 0.3);
		glBegin(GL_POINTS);
		for(int i = 0; i < vn; ++i){
			glVertex3fv(glm::value_ptr(poly.vertices[i]));
		}
		glEnd();
	}

	// エッジ描画
	if(draw & RXD_EDGE){
		glDisable(GL_LIGHTING);
		glColor3d(0.5, 0.5, 0.5);
		glLineWidth(1.0);
		for(int i = 0; i < pn; ++i){
			const rxFace *face = &poly.faces[i];
			int n = (int)face->vert_idx.size();

			glBegin(GL_LINE_LOOP);
			for(int j = 0; j < n; ++j){
				int idx = face->vert_idx[j];
				if(idx >= 0 && idx < vn){
					glVertex3fv(glm::value_ptr(poly.vertices[idx]));
				}
			}
			glEnd();
		}
	}
}
