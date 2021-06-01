/*!
  @file scene.h

  @brief GLFWによるOpenGL描画

  @author Makoto Fujisawa
  @date   2021-04
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

const int MAX_STEPS = 1000000;

// 描画フラグ
const int RXD_VERTEX = 0x0001;	//!< 頂点描画
const int RXD_EDGE = 0x0002;	//!< エッジ描画
const int RXD_FACE = 0x0004;	//!< 面描画

const int RXD_BBOX = 0x0010;	//!< AABB(シミュレーション空間)
const int RXD_AXIS = 0x0020;   //!< 軸
const int RXD_FLOOR = 0x0040;	//!< 床

static float g_dt = 0.002;	//!< 時間ステップ幅


//-----------------------------------------------------------------------------
//! SceneSWEクラス
//   - SWEによる波シミュレーション
//   - GUI関係の処理全般を行う(イベントハンドラ全般)
//-----------------------------------------------------------------------------
class SceneSWE
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

	//! マウスピック
	static int m_picked;					//!< ピックされたグリッドセル番号
	static float m_pickdist;				//!< ピックされた点までの距離

	// ハイトフィールド
	static WaveSWE *m_wave;

	// シーンパラメータ
	static float m_scale;
	static int m_res;

public:
	//! コンストラクタ
	SceneSWE(){}

	//! デストラクタ
	~SceneSWE(){}

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
	// 波の初期化
	static void initWave(void);
	static void initWaveSlope(void);
	static void initWaveMountain(void);

	static float getHeightFlat(float x, float y);
	static float getHeightSlope(float x, float y);
	static float getHeightMountain(float x, float y);

private:
	// アニメーション切り替え
	static bool switchanimation(int on);

	// ファイル入出力
	static void savedisplay(const int &stp);

	// 視点
	static void resetview(void);
};




#endif // #ifdef _RX_CONTROLLER_H_
