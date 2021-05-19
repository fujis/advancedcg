/*! 
  @file rx_wave.h
	
  @brief SWEを使った波のシミュレーション
 
  @author Makoto Fujisawa
  @date 2018-12
*/

#ifndef _RX_WAVE_H_
#define _RX_WAVE_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <ctime>

// STL
#include <vector>
#include <string>

// OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

// 各種ユーティリティ
#include "rx_utility.h"

// ポリゴンメッシュ
#include "rx_mesh.h"


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
using namespace std;


//-----------------------------------------------------------------------------
// rxWave
//  - SWEによる波のシミュレーション
//  - A. T. Layton and M. Van De Panne, “A Numerically Efficient and Stable Algorithm for Animating Water Waves”,
//    The Visual Computer, Vol. 18, No. 1, pp. 41-53, 2002.
//-----------------------------------------------------------------------------
class rxWave
{
	vector<double> m_vH;			//!< ハイトフィールド
	vector<double> m_vHprev;		//!< 前ステップのハイトフィールド
	vector<double> m_vU, m_vV;		//!< 速度場(u,v)
	vector<double> m_vUprev, m_vVprev;	//!< 前ステップの速度場(u,v)
	rxPolygons m_Mesh;				//!< ハイトフィールドメッシュ

	int m_iNx, m_iNy;				//!< グリッド分割数
	double m_fScale;				//!< 全体のスケール

	double m_fGravity;				//!< 重力
	double m_fDx, m_fDy;			//!< グリッド幅
	double m_fWaterDepth;			//!< 水深

	int m_iStep;					//!< 計算ステップ

public:
	//! コンストラクタ
	rxWave(int n, double scale, double depth);
	
	//! デストラクタ
	~rxWave();

public:
	//! 初期化
	void Init(void);
	//! 更新
	int Update(double dt);
	//! 描画
	void Draw(int drw = 2+4);
	void DrawForPick(void);

	//! アクセスメソッド
	rxPolygons& GetMesh(void){ return m_Mesh; }
	int& GetStep(void){ return m_iStep; }

	//! グリッドインデックスの計算
	inline int IDX(int i, int j){ return i+j*(m_iNx); }
	inline int IDX(int i, int j, int n){ return i+j*n; }

	//! 波の初期化(式(25)の計算)
	void InitWave(void);

	//! 高さ値の直接変更
	void AddHeight(int idx, double h);
	void AddHeight(int i, int j, double h);

	//! ハイトフィールドデータの入出力
	void InputHeightField(const string &fn);
	void OutputHeightField(const string &fn);

protected:
	//! n×nの頂点を持つメッシュ生成(x-z平面)
	void generateMesh(Vec3 c1, Vec3 c2);

	//! ハイトフィールドメッシュの更新
	void updateMesh(const vector<double> &h);

	//! SWEによるハイトフィールドの更新
	void updateHeightBySWE(double *h_new, double *h, double *u_new, double *v_new, double *u, double *v, double dt);

	//! ハイトフィールドにガウスフィルタを適用
	void gaussianHeight(const vector<double> &h, vector<double> &hs, int size = 9, int sigma = 1.0);

	// ハイトフィールドデータの入出力
	void inputHeightField(const string &fn, vector<double> &hf, int &nx, int &ny);
	void outputHeightField(const string &fn, const vector<double> &hf, int nx, int ny);

};




#endif // #ifdef _RX_WAVE_H_


