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

// テクスチャ
#include "rx_texture.h"

// メッシュ変形
#include "deformation.h"



using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
const int MAX_STEPS = 1000000;

// 描画フラグ
const int RXD_VERTEX = 0x0001;	//!< 頂点描画
const int RXD_EDGE = 0x0002;   //!< エッジ描画
const int RXD_FACE = 0x0004;	//!< 面描画
const int RXD_FIX = 0x0008;	//!< 固定点描画
const int RXD_TEXTURE = 0x0010;	//!< テクスチャ描画


//-----------------------------------------------------------------------------
//! SceneMLSクラス
//   - MLS変形シーン
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

	static rxMeshDeform2D* m_md;			//!< メッシュ変形
	static glm::vec2 m_envmin, m_envmax;	//!< メッシュ描画範囲
	static GLuint m_tex;					//!< テクスチャ

	static int m_picked;					//!< マウスピックされた頂点番号

	static rxGLSL m_shader;					//!< GLSLシェーダ

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
