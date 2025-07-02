/*! 
  @file wave.h
	
  @brief SWEを使った波のシミュレーション
 
  @author Makoto Fujisawa
  @date 2021-05
*/

#ifndef _RX_WAVE_H_
#define _RX_WAVE_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"
#include "rx_mesh.h"


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
using namespace std;


//-----------------------------------------------------------------------------
//  SWEによる波のシミュレーション
//  - 参考: A. T. Layton and M. Van De Panne, “A Numerically Efficient and Stable Algorithm for Animating Water Waves”,
//    The Visual Computer, Vol. 18, No. 1, pp. 41-53, 2002.
//-----------------------------------------------------------------------------
class WaveSWE
{
	vector<float> m_h;				//!< ハイトフィールド
	vector<float> m_d;				//!< ハイトフィールド
	vector<float> m_dprev;			//!< 前ステップのハイトフィールド
	vector<float> m_u, m_v;		//!< 速度場(u,v)
	vector<float> m_uprev, m_vprev;//!< 前ステップの速度場(u,v)

	rxPolygons m_mesh;				//!< 水面ハイトフィールドメッシュ(描画用)
	rxPolygons m_mesh_ground;		//!< 水底ハイトフィールドメッシュ(描画用)
	glm::vec3 m_min, m_max;			//!< シミュレーション空間の描画上の大きさ

	int m_nx, m_ny;					//!< グリッド分割数
	float m_scale;					//!< 全体のスケール

	float m_gravity;				//!< 重力
	float m_dx, m_dy;				//!< グリッド幅
	float m_avg_h;					//!< 平均高さ

	float(*m_ground)(float, float);	//!< 水底の高さを返す関数
public:
	bool m_surf;					//!< 打ち寄せる波の生成フラグ
	bool m_gs;						//!< 地形の高さの影響を陰解法(ガウス・ザイデル)で解くフラグ
	float m_nu;						//!< 動粘性係数

public:
	//! コンストラクタ
	WaveSWE();
	
	//! デストラクタ
	~WaveSWE();

	//! グリッドインデックスの計算
	inline int IDX(int i, int j){ return i+j*(m_nx); }
	inline int IDX(int i, int j, int n){ return i+j*n; }

public:
	//! 初期化
	void Init(int n, float scale, float (*ground)(float,float));
	
	//! 更新
	int Update(int step, float dt);

	//! 描画
	void Draw(int draw = 2+4);

	//! アクセスメソッド
	rxPolygons& GetMesh(void){ return m_mesh; }
	void GetDim(glm::vec3 &cen, glm::vec3 &dim){ cen = 0.5f*(m_min+m_max); dim = (m_max-m_min); }

	//! 高さ値の直接変更
	void AddHeight(int idx, float h);
	void AddHeight(int i, int j, float h);
	void AddHeightArea(int idx, float h);

	// 頂点選択(レイと頂点(球)の交差判定)
	int IntersectRay(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, float &t, float rad = 0.05);

protected:
	//! n×nの頂点を持つメッシュ生成(x-z平面)
	void generateMesh(glm::vec3 c1, glm::vec3 c2, rxPolygons &poly);

	//! ハイトフィールドメッシュの更新
	void updateMesh(const vector<float> &h);

	//! SWEによるハイトフィールドの更新
	void advection(float *d_new, float *d, float *u_new, float *v_new, float *u, float *v, float dt);
	void pressure(float *d_new, float *d, float *u_new, float *v_new, float *u, float *v, float dt);
	void viscosity(float *u_new, float *v_new, float *u, float *v, float dt);

	//! 水深と地形の高さから水面の高さを更新
	void updateHeight(float *d);

	//! 境界条件
	void bnd(float *d);
	void bnd(float *u, float *v);

	//! 打ち寄せる波の生成
	void makeSurf(float t, float *h, float wave_height);

	//! 平均高さの算出
	double calAvarageHeight(void);

};




#endif // #ifdef _RX_WAVE_H_


