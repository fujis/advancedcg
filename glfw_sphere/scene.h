/*!
  @file rx_controller.h
	
  @brief GLUTによるOpenGLウィンドウクラス
 
  @author Makoto Fujisawa 
  @date   2020-06
*/

#ifndef _RX_CONTROLLER_H_
#define _RX_CONTROLLER_H_



//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"
#include "imgui.h"

// テクスチャ・画像
#include "rx_texture.h"
#include "rx_bitmap.h"

// 視点移動
#include "rx_trackball.h"

// メッシュファイル読み込み関連
#include "rx_obj.h"


using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
const int MAX_STEPS = 1000000;
const float DT = 0.01;

// 描画フラグ
const int RXD_VERTEX   = 0x0001;	//!< 頂点描画
const int RXD_EDGE     = 0x0002;   //!< エッジ描画
const int RXD_FACE     = 0x0004;	//!< 面描画
const int RXD_SPHERE   = 0x0008;
const int RXD_CYLINDER = 0x0010;
const int RXD_CAPSULE  = 0x0020;

const int RXD_BBOX  = 0x0100;	//!< AABB(シミュレーション空間)
const int RXD_AXIS  = 0x0200;   //!< 軸
const int RXD_FLOOR = 0x0400;	//!< 床

struct Primitive
{
	int nvrts, ntris;
	GLuint vao;
};


//-----------------------------------------------------------------------------
//! rxControllerクラス
//   - GLUT(freeglut)によるアプリケーション全体のコントローラ
//   - GUI関係の処理全般を行う(イベントハンドラ全般)
//-----------------------------------------------------------------------------
class SceneMLS
{
protected:
	static int m_winw;						//!< 描画ウィンドウの幅
	static int m_winh;						//!< 描画ウィンドウの高さ
	static double m_bgcolor[3];				//!< 背景色

	static bool m_animation_on;				//!< アニメーションON/OFF

	static int m_draw;						//!< 描画フラグ
	static int m_currentstep;				//!< 現在のステップ数
	static int m_simg_spacing;				//!< 画像保存間隔(=-1なら保存しない)
	static int m_picked;					//!< マウスピックされた頂点番号

	static rxTrackball m_view;				//!< 視点移動トラックボール

	static GLuint m_tex;					//!< テクスチャ

	static rxGLSL m_shader;					//!< GLSLシェーダ

	static Primitive m_sphere;				//!< VAO:球
	static Primitive m_cylinder;			//!< VAO:円筒
	static Primitive m_capsule;				//!< VAO:カプセル

public:
	//! コンストラクタ
	SceneMLS(){}

	//! デストラクタ
	~SceneMLS(){}

public:
	// コールバック関数
	static void Init(int argc, char* argv[]);
	static void Draw();
	static void Timer();
	static void Cursor(GLFWwindow* window, double xpos, double ypos);
	static void Mouse(GLFWwindow* window, int button, int action, int mods);
	static void Keyboard(GLFWwindow* window, int key, int mods);
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui(GLFWwindow* window);
	static void Destroy();

private:
	// アニメーション切り替え
	static bool switchanimation(int on);

	// ファイル入出力
	static void savedisplay(const int &stp);
};




#endif // #ifdef _RX_CONTROLLER_H_
