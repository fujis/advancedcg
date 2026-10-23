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
	static glm::vec2 m_envmin, m_envmax;	//!< メッシュ描画範囲(ウィンドウの縦横比に合わせて広げて表示)
	static glm::vec2 m_view_center;			//!< 描画範囲の中心(注視点)
	static bool m_panning;					//!< Ctrl+左ドラッグで注視点移動中
	static double m_pan_x, m_pan_y;			//!< 注視点移動中の直前のマウス座標
	static GLuint m_tex;					//!< テクスチャ
	static string m_imgfile;				//!< メッシュ生成・テクスチャに使う画像ファイル名
	static vector<string> m_imglist;		//!< 実行時フォルダにある画像ファイルのリスト
	static rxImageMeshParams m_img_prm;		//!< 画像からのメッシュ生成パラメータ

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
	// 描画範囲(ウィンドウの縦横比に合わせる)とマウス座標->ワールド座標変換
	static void viewRange(glm::vec2 &vmin, glm::vec2 &vmax);
	static glm::vec2 screenToWorld(GLFWwindow* window, double x, double y);

	// 実行時フォルダの画像ファイルをリストアップ
	static void scanImages(void);

	// 画像を読み込んでテクスチャとメッシュを作成
	static bool loadImageMesh(void);

	// アニメーション切り替え
	static bool switchanimation(int on);

	// ファイル入出力
	static void savedisplay(const int &stp);
};




#endif // #ifdef _RX_CONTROLLER_H_
