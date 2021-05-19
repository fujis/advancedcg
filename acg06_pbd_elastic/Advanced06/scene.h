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

// PBD
#include "pbd.h"

// メッシュファイル読み込み関連
#include "rx_obj.h"


using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
#define RX_OUTPUT_TIME

const int MAX_STEPS = 1000000;
const float GOUND_HEIGHT = -1.0;

const float DT = 0.01;


// 描画フラグ
const int RXD_VERTEX = 0x0001;
const int RXD_EDGE = 0x0002;
const int RXD_FACE = 0x0004;
const int RXD_TETRA = 0x0008;	//!< 四面体(ワイヤフレーム)
const int RXD_TETCS = 0x0010;	//!< 四面体(断面)
const int RXD_BBOX = 0x0020;	//!< AABB(シミュレーション空間)
const int RXD_AXIS = 0x0040;   //!< 軸
const int RXD_FLOOR = 0x0080;	//!< 床


//-----------------------------------------------------------------------------
//! ScenePBDクラス
//   - PBDによる弾性変形シミュレーション
//   - GUI関係の処理全般を行う(イベントハンドラ全般)
//-----------------------------------------------------------------------------
class ScenePBD
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

	// PBDによる弾性変形
	static float m_dt;
	static ElasticPBD* m_elasticbody;

	//! マウスピック
	static int m_picked;
	static float m_pickdist;	//!< ピックされた点までの距離

public:
	//! コンストラクタ
	ScenePBD(){}

	//! デストラクタ
	~ScenePBD(){}

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

	// PBD
	static void initRod(void);
	static void initCloth(void);
	static void initBall(void);

	// マウスピック
	static void clearPick(void);
};




#endif // #ifdef _RX_CONTROLLER_H_
