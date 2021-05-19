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


// トラックボール＆テクスチャ
#include "rx_trackball.h"
#include "rx_texture.h"

// 流体シミュレーション
#include "wave.h"

// メッシュファイル読み込み関連
#include "rx_obj.h"


using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
#define RX_OUTPUT_TIME

const string RX_PROGRAM_NAME = "acg";
const int RX_MAX_STEPS = 1000000;

// 描画フラグ
static int RXD_VERTEX = 0x0001;	//!< 頂点描画
static int RXD_EDGE = 0x0002;	//!< エッジ描画
static int RXD_FACE = 0x0004;	//!< 面描画
static int RXD_NORMAL = 0x0008;	//!< 法線描画
static int RXD_BBOX = 0x0010;	//!< AABB(シミュレーション空間)
static int RXD_AXIS = 0x0020;   //!< 軸
static int RXD_FLOOR = 0x0040;	//!< 床



//-----------------------------------------------------------------------------
//! rxControllerクラス
//   - GLUT(freeglut)によるアプリケーション全体のコントローラ
//   - GUI関係の処理全般を行う(イベントハンドラ全般)
//-----------------------------------------------------------------------------
class SceneMLS
{
protected:
	static int m_winid;						//!< GLUTのウインドウID
	static int m_winw;						//!< 描画ウィンドウの幅
	static int m_winh;						//!< 描画ウィンドウの高さ
	static int m_mousebutton;				//!< マウスボタンの状態
	static int m_keymod;					//!< 修飾キーの状態
	static rxTrackball m_view;				//!< トラックボール

	static double m_bgcolor[3];				//!< 背景色
	static bool m_animation_on;				//!< アニメーションON/OFF
	static bool m_fullscreen_on;			//!< フルスクリーンON/OFF

	static int m_currentstep;				//!< 現在のステップ数
	static bool m_pause;					//!< シミュレーションのポーズフラグ

	// ハイトフィールド
	static rxWave *m_wave;
	static double m_dt;	//!< 時間ステップ幅

	static int m_view_drag;

	// メッシュ
	static rxPolygonsE m_poly_org;
	static rxPolygonsE m_poly;
	static GLuint m_tex;


public:
	// 描画フラグ
	static int m_draw;						//!< 描画フラグ

	// 画像出力
	static int m_simg_spacing;		//!< 画像保存間隔(=-1なら保存しない)

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

	// 視点
	static void resetview(void);

	// 描画関係
	static vector<string> setDrawString(void);

public:
	// マウスピック用描画
	static void RenderSceneForPick(void* x = 0);
	static void ProjectionForPick(void* x = 0);
};




#endif // #ifdef _RX_CONTROLLER_H_
