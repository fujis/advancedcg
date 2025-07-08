/*!
  @file deformation.h

  @brief 2Dメッシュ変形

  @author Makoto Fujisawa
  @date 2021-03
*/

#ifndef _RX_DEFORM_H_
#define _RX_DEFORM_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"

//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
using namespace std;


//-----------------------------------------------------------------------------
// rxMeshDeform2D
//-----------------------------------------------------------------------------
class rxMeshDeform2D
{
	vector<glm::vec2> m_vX;			//!< ポリゴン頂点座標
	vector<glm::vec2> m_vP;			//!< 初期座標
	vector<glm::vec2> m_vTC;		//!< テクスチャ座標

	vector<unsigned int> m_vTri;	//!< 三角形ポリゴン頂点インデックス
	int m_iNt, m_iNv;				//!< ポリゴン数，頂点数

	vector<int> m_vCP;				//!< 制御点リスト
	int m_iNcp;					//!< 制御点数

	GLuint m_vao_mesh;				//!< メッシュデータのVAO
	GLuint m_vao_fix;				//!< 固定頂点のためのVAO

public:
	double m_alpha;					//!< 重み計算のための係数(w=1/(v-p)^(2*alpha)
	int m_deform_type;				//!< 変形のタイプ(0:affine, 1:similarity, 2:rigid)

public:
	//! コンストラクタ
	rxMeshDeform2D();

	//! デストラクタ
	~rxMeshDeform2D();

public:
	//! 初期化
	void Init(int random_mesh = 0);
	//! 更新
	int Update(double dt);
	//! 描画
	void DrawMesh(void);
	void DrawPoints(void);
	void DrawFixPoints(void);
	void InitVAO(void);

	//! 近傍頂点,近傍制御点探索
	int Search(glm::vec2 pos, double h = 0.05);
	int SearchCP(glm::vec2 pos, double h = 0.05);

	//! 固定点設定
	void SetCP(int idx, glm::vec2 pos, bool move = false);
	//! 固定点解除
	void UnsetCP(int idx);

protected:
	//! n×nの頂点を持つメッシュ生成(x-z平面)
	void generateMesh(glm::vec2 c1, glm::vec2 c2, int nx, int ny);

	//! ポアソンディスクサンプリングによるメッシュ生成
	void generateRandomMesh(glm::vec2 c1, glm::vec2 c2, double min_dist, int n);

	//! 固定点座標VAOの更新
	void updateCPVAO(void);

	//! グリッドインデックスの計算
	inline int IDX(int i, int j, int n){ return i+j*n; }

	// MLS mesh deformation
	glm::vec2 affineDeformation(const glm::vec2 &v, const glm::vec2 &pc, const glm::vec2 &qc, const double alpha);		//!< Affine Deformation
	glm::vec2 similarityDeformation(const glm::vec2 &v, const glm::vec2 &pc, const glm::vec2 &qc, const double alpha);	//!< Similarity Deformation
	glm::vec2 rigidDeformation(const glm::vec2 &v, const glm::vec2 &pc, const glm::vec2 &qc, const double alpha);		//!< Rigid Deformation

};




#endif // #ifdef _RX_DEFORM_H_


