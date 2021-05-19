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

GLuint SceneMLS::m_tex = 0;						//!< テクスチャ
int SceneMLS::m_picked = -1;					//!< マウスピックされた頂点番号

rxTrackball SceneMLS::m_view;					//!< トラックボール

rxGLSL SceneMLS::m_shader;						//!< GLSLシェーダ

Primitive SceneMLS::m_sphere;
Primitive SceneMLS::m_cylinder;
Primitive SceneMLS::m_capsule;



static GLuint makeSphereMesh(int &nvrts, int &ntris, float rad, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();
	vector<glm::vec3> vrts, nrms;
	vector<GLuint> tris;

	for(int j = 0; j <= stacks; ++j){
		float t = float(j)/float(stacks);
		float y = rad*cos(pi*t);
		float rj = rad*sin(pi*t);	// 高さyでの球の断面円半径
		for(int i = 0; i <= slices; ++i){
			float s = float(i)/float(slices);
			float x = rj*sin(2*pi*s);
			float z = rj*cos(2*pi*s);
			vrts.push_back(glm::vec3(x, y, z));
			nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
		}
	}
	nvrts = static_cast<int>(vrts.size());
	// メッシュ作成
	int nx = slices+1;
	for(int j = 0; j < stacks; ++j){
		for(int i = 0; i < slices; ++i){
			tris.push_back((i)+(j)*nx);
			tris.push_back((i+1)+(j+1)*nx);
			tris.push_back((i+1)+(j)*nx);

			tris.push_back((i)+(j)*nx);
			tris.push_back((i)+(j+1)*nx);
			tris.push_back((i+1)+(j+1)*nx);
		}
	}
	ntris = static_cast<int>(tris.size()/3);

	return CreateVAO((GLfloat*)(&vrts[0]), nvrts, 3, &tris[0], ntris, (GLfloat*)(&nrms[0]), nvrts);
}

static GLuint makeCylinderMesh(int &nvrts, int &ntris, float rad, float len, int slices = 16)
{
	const float pi = glm::pi<float>();
	vector<glm::vec3> vrts, nrms;
	vector<GLuint> tris;

	for(int i = 0; i <= slices; ++i){
		float t = float(i)/float(slices);
		float x = rad*cos(2*pi*t);
		float y = rad*sin(2*pi*t);
		vrts.push_back(glm::vec3(x, y, -0.5*len));
		nrms.push_back(glm::normalize(glm::vec3(x, y, 0.0)));
		vrts.push_back(glm::vec3(x, y,  0.5*len));
		nrms.push_back(glm::normalize(glm::vec3(x, y, 0.0)));

	}
	nvrts = static_cast<int>(vrts.size());

	// メッシュ作成
	for(int i = 0; i < 2*slices; i += 2){
		tris.push_back(i);
		tris.push_back((i+2 >= 2*slices ? 0 : i+2));
		tris.push_back(i+1);

		tris.push_back(i+1);
		tris.push_back((i+2 >= 2*slices ? 0 : i+2));
		tris.push_back((i+2 >= 2*slices ? 1 : i+3));
	}
	ntris = static_cast<int>(tris.size()/3);

	return CreateVAO((GLfloat*)(&vrts[0]), nvrts, 3, &tris[0], ntris, (GLfloat*)(&nrms[0]), nvrts);
}

static GLuint makeCapsuleMesh(int &nvrts, int &ntris, float rad, float len, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();
	vector<glm::vec3> vrts, nrms;
	vector<GLuint> tris;
	int offset = 0;

	//// 円筒部分頂点
	//for(int i = 0; i <= slices; ++i){
	//	float t = float(i)/float(slices);
	//	float x = rad*cos(2*pi*t);
	//	float y = rad*sin(2*pi*t);
	//	vrts.push_back(glm::vec3(x, y, -0.5*len));
	//	nrms.push_back(glm::normalize(glm::vec3(x, y, 0.0)));
	//	vrts.push_back(glm::vec3(x, y, 0.5*len));
	//	nrms.push_back(glm::normalize(glm::vec3(x, y, 0.0)));

	//}
	//offset = static_cast<int>(vrts.size());

	//// 円筒部分メッシュ作成
	//for(int i = 0; i < 2*slices; i += 2){
	//	tris.push_back(i);
	//	tris.push_back((i+2 >= 2*slices ? 0 : i+2));
	//	tris.push_back(i+1);

	//	tris.push_back(i+1);
	//	tris.push_back((i+2 >= 2*slices ? 0 : i+2));
	//	tris.push_back((i+2 >= 2*slices ? 1 : i+3));
	//}

	// 球体部分頂点
	for(int j = 0; j <= stacks; ++j){
		float t = float(j)/float(stacks);
		float z = rad*cos(pi*t);
		float rj = rad*sin(pi*t);	// 高さyでの球の断面円半径
		float zlen = (j < stacks/2 ? 0.5*len : -0.5*len);
		if(j == stacks/2){
			for(int i = 0; i <= slices; ++i){
				float s = float(i)/float(slices);
				float x = rj*sin(2*pi*s);
				float y = rj*cos(2*pi*s);
				vrts.push_back(glm::vec3(x, y, z-zlen));
				nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
			}
		}
		for(int i = 0; i <= slices; ++i){
			float s = float(i)/float(slices);
			float x = rj*sin(2*pi*s);
			float y = rj*cos(2*pi*s);
			vrts.push_back(glm::vec3(x, y, z+zlen));
			nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
		}
	}
	// メッシュ作成
	int nx = slices+1;
	for(int j = 0; j < stacks+1; ++j){
		for(int i = 0; i < slices; ++i){
			tris.push_back((i)+(j)*nx+offset);
			tris.push_back((i+1)+(j)*nx+offset);
			tris.push_back((i+1)+(j+1)*nx+offset);

			tris.push_back((i)+(j)*nx+offset);
			tris.push_back((i+1)+(j+1)*nx+offset);
			tris.push_back((i)+(j+1)*nx+offset);
		}
	}

	nvrts = static_cast<int>(vrts.size());
	ntris = static_cast<int>(tris.size()/3);

	return CreateVAO((GLfloat*)(&vrts[0]), nvrts, 3, &tris[0], ntris, (GLfloat*)(&nrms[0]), nvrts);
}



static inline void drawPrimitiveVAO(const Primitive &obj, int draw)
{
	// エッジ描画における"stitching"をなくすためのオフセットの設定
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	// 図形の描画
	glDisable(GL_LIGHTING);
	glBindVertexArray(obj.vao);
	if(draw & 0x01){
		glColor3d(1.0, 0.3, 0.3);
		glPointSize(5.0);
		glDrawArrays(GL_POINTS, 0, obj.nvrts);
	}
	if(draw & 0x02){
		glColor3d(0.5, 0.9, 0.9);
		glLineWidth(4.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, obj.ntris*3, GL_UNSIGNED_INT, 0);
	}
	if(draw & 0x04){
		glEnable(GL_LIGHTING);
		//glDisable(GL_CULL_FACE);
		glEnable(GL_AUTO_NORMAL);
		glEnable(GL_NORMALIZE);
		glColor3d(0.1, 0.5, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, obj.ntris*3, GL_UNSIGNED_INT, 0);
	}
	glBindVertexArray(0);
}


static inline void drawCapsule(const Primitive &cylinder, const Primitive &sphere, int draw)
{
}

/*!
 * 初期化関数
 */
void SceneMLS::Init(int argc, char* argv[])
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

	// 描画フラグ初期化
	m_draw = 0;
	m_draw |= RXD_EDGE;
	m_draw |= RXD_FACE;
	m_draw |= RXD_FLOOR;
	m_draw |= RXD_BBOX;

	// シェーダの初期化
	//m_shader = CreateGLSLFromFile("shaders/shading.vp", "shaders/shading.fp", "shading");

	//// テクスチャ
	//std::string filename = "sample.bmp";
	//glActiveTexture(GL_TEXTURE0);
	//if(!loadTexture(filename, m_tex, false, false)){
	//	loadTexture("./bin/"+filename, m_tex, false, false);
	//}

	// メッシュ初期化
	m_sphere.vao = makeSphereMesh(m_sphere.nvrts, m_sphere.ntris, 1.0, 16, 8);
	m_cylinder.vao = makeCylinderMesh(m_cylinder.nvrts, m_cylinder.ntris, 0.55, 2.2, 16);
	m_capsule.vao = makeCapsuleMesh(m_capsule.nvrts, m_capsule.ntris, 0.5, 2.2, 16);

	// トラックボール初期姿勢
	m_view.SetTranslation(0, 0);
	m_view.SetScaling(-5);

	// アニメーションON
	switchanimation(0);
}



/*!
 * 再描画イベント処理関数
 */
void SceneMLS::Draw(void)
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

	// 図形の描画
	if(m_draw & RXD_SPHERE) drawPrimitiveVAO(m_sphere, m_draw);
	if(m_draw & RXD_CYLINDER) drawPrimitiveVAO(m_cylinder, m_draw);
	if(m_draw & RXD_CAPSULE) drawPrimitiveVAO(m_capsule, m_draw);
	//if(m_draw & RXD_CAPSULE) drawCapsule(m_sphere, m_cylinder, m_draw);


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
void SceneMLS::Cursor(GLFWwindow* window, double x, double y)
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
void SceneMLS::Keyboard(GLFWwindow* window, int key, int mods)
{
	switch(key){
	case GLFW_KEY_R: // リセット
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
	if(ImGui::Button("save screenshot")){ savedisplay(-1); }
	ImGui::Separator();
	ImGui::CheckboxFlags("vertex", &m_draw, RXD_VERTEX);
	ImGui::CheckboxFlags("edge", &m_draw, RXD_EDGE);
	ImGui::CheckboxFlags("face", &m_draw, RXD_FACE);
	ImGui::CheckboxFlags("aabb", &m_draw, RXD_BBOX);
	ImGui::CheckboxFlags("axis", &m_draw, RXD_AXIS);
	ImGui::CheckboxFlags("floor", &m_draw, RXD_FLOOR);
	ImGui::Separator();
	ImGui::CheckboxFlags("sphere", &m_draw, RXD_SPHERE);
	ImGui::CheckboxFlags("cylinder", &m_draw, RXD_CYLINDER);
	ImGui::CheckboxFlags("capsule", &m_draw, RXD_CAPSULE);
	ImGui::Separator();
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

/*!
 * 終了処理関数
 */
void SceneMLS::Destroy()
{
	if(m_sphere.vao) glDeleteVertexArrays(1, &m_sphere.vao);
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



