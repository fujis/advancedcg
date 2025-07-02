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

// BVHファイル
#include "characteranimation.h"

// メッシュファイル読み込み関連
#include "rx_obj.h"


using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
#define RX_OUTPUT_TIME

const int RX_MAX_STEPS = 1000000;
const float GOUND_HEIGHT = -1.0;

const float DT = 0.01;

// 描画フラグ
const int RXD_VERTEX = 0x0001;
const int RXD_EDGE = 0x0002;
const int RXD_FACE = 0x0004;
const int RXD_BBOX = 0x0010;		//!< AABB(シミュレーション空間)
const int RXD_AXIS = 0x0020;		//!< 軸
const int RXD_FLOOR = 0x0040;		//!< 床
const int RXD_SKELETON = 0x0100;	//!< スケルトン描画
const int RXD_MESH = 0x0200;		//!< メッシュ描画
const int RXD_WEIGHT = 0x0400;		//!< 頂点の重み描画



//-----------------------------------------------------------------------------
//! SceneCAクラス
//   - スケルトンを用いたキャラクターアニメーション
//   - GUI関係の処理全般を行う(イベントハンドラ全般)
//-----------------------------------------------------------------------------
class SceneCA
{
protected:
	static int m_winw;						//!< 描画ウィンドウの幅
	static int m_winh;						//!< 描画ウィンドウの高さ
	static rxTrackball m_view;				//!< トラックボール

	static float m_bgcolor[3];				//!< 背景色
	static bool m_animation_on;				//!< アニメーションON/OFF

	static int m_draw;						//!< 描画フラグ
	static int m_currentstep;				//!< 現在のステップ数
	static int m_simg_spacing;				//!< 画像保存間隔(=-1なら保存しない)

	// キャラクターアニメーション関連
	static CharacterAnimation *m_ca;
	static float m_scale;
	static glm::vec3 m_offset;

	//! マウスピック
	static int m_picked;
	static float m_pickdist;	//!< ピックされた点までの距離

	// メッシュ(Skinning用)
	static rxPolygonsE m_poly_org;
	static rxPolygonsE m_poly;
	static GLuint m_tex;

public:
	//! コンストラクタ
	SceneCA(){}

	//! デストラクタ
	~SceneCA(){}

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

	// キャラクターアニメーションの初期化
	static void initCA(const string &bvh, const string &obj);

	// ポリゴン描画(重みによる色分け切り替えあり)
	static void drawPolygon(rxPolygonsE &poly, int draw = 0x04, float dn = 0.02f, bool col = true);
};




#endif // #ifdef _RX_CONTROLLER_H_
