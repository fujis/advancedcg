/*!
  @file pbd.h
	
  @brief PBDによる弾性体シミュレーション
 
  @author Makoto Fujisawa
  @date 2021
*/

#ifndef _RX_PBD_STRAND_H_
#define _RX_PBD_STRAND_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"

#include "rx_mesh.h"

using namespace std;

//! 衝突判定用関数
typedef void (*CollisionFunc)(glm::vec3&, glm::vec3&, glm::vec3&, int);

//-----------------------------------------------------------------------------
// ElasticPBDクラス
//  - PBDによる弾性体シミュレーション
//-----------------------------------------------------------------------------
class ElasticPBD
{
protected:
	// 形状データ
	int m_iNumVertices;					//!< 頂点数
	vector<glm::vec3> m_vCurPos;		//!< 現在の頂点位置
	vector<glm::vec3> m_vNewPos;		//!< 次のステップの頂点位置

	vector<glm::vec3> m_vVel;			//!< 頂点速度
	vector<bool> m_vFix;				//!< 頂点固定フラグ
	vector<float> m_vMass;				//!< 頂点質量(変形時の重み)

	rxPolygonsE m_poly;					//!< エッジ情報を含む3角形ポリゴン(元の頂点位置含む)
	int m_iNumEdge;						//!< エッジ数
	int m_iNumTris;						//!< ポリゴン数

	vector< vector<int> > m_vTets;		//!< 四面体情報
	int m_iNumTets;						//!< 四面体数

	vector<float> m_vLengths;			//!< 頂点間の初期長さ(stretch constraint用)
	vector<float> m_vBends;				//!< エッジを挟むメッシュ間の初期角度(bending constraint用)
	vector<float> m_vVolumes;			//!< 四面体の初期体積(volume constraint用)

	vector<int> m_vInEdge;				//!< 内部エッジならtrue

	// シミュレーションパラメータ
	glm::vec3 m_v3Min, m_v3Max;			//!< シミュレーション空間の大きさ
	glm::vec3 m_v3Gravity;				//!< 重力加速度ベクトル

	// 衝突処理関数の関数ポインタ
	CollisionFunc m_fpCollision;

public:
	int m_iNmax;						//!< 最大反復回数
	float m_fK;							//!< stiffnessパラメータ[0,1] 
	float m_fWind;						//!< 風の強さ

	bool m_bUseInEdge;					//!< 物体内部のエッジを制約計算に使うかどうかのフラグ
	int m_iObjectNo;					//!< オブジェクト番号

public:
	//! コンストラクタとデストラクタ
	ElasticPBD(int obj);
	~ElasticPBD();

	void Clear();

	void AddVertex(const glm::vec3 &pos, float mass);
	void AddEdge(int v0, int v1);
	void AddTriangle(int v0, int v1, int v2);
	void AddTetra(int v0, int v1, int v2, int v3);

	void MakeEdge(void);
	void CalBends(void);

	// ストランド・メッシュ・四面体生成
	void GenerateStrand(glm::vec3 c1, glm::vec3 c2, int n);
	void GenerateMesh(glm::vec2 c1, glm::vec2 c2, int nx, int ny);
	void GenerateTetrahedralMeshFromFile(const string &filename);

	void Update(float dt);

	// OpenGL描画
	void Draw(int drw);

	// アクセスメソッド
	void SetSimulationSpace(glm::vec3 minp, glm::vec3 maxp){ m_v3Min = minp; m_v3Max = maxp; }
	void SetStiffness(float k){ m_fK = k; }
	void SetCollisionFunc(CollisionFunc func){ m_fpCollision = func; }

	glm::vec3 GetVertexPos(int i){ return ((i >= 0 && i < m_iNumVertices) ? m_vCurPos[i] : glm::vec3(0.0)); }

	// 固定点設定
	void FixVertex(int i, const glm::vec3 &pos);
	void FixVertex(int i);
	void UnFixVertex(int i);
	void UnFixAllVertex(void);
	bool IsFixed(int i) { return m_vFix[i]; }

	// 頂点選択(レイと頂点(球)の交差判定)
	int IntersectRay(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, float &t, float rad = 0.05);

protected:
	// PBD計算
	void calExternalForces(float dt);
	void genCollConstraints(float dt);
	void projectStretchingConstraint(float ks);
	void projectBendingConstraint(float ks);
	void projectVolumeConstraint(float ks);
	void integrate(float dt);


	void clamp(glm::vec3 &pos) const
	{
		if(pos[0] < m_v3Min[0]) pos[0] = m_v3Min[0];
		if(pos[0] > m_v3Max[0]) pos[0] = m_v3Max[0];
		if(pos[1] < m_v3Min[1]) pos[1] = m_v3Min[1];
		if(pos[1] > m_v3Max[1]) pos[1] = m_v3Max[1];
		if(pos[2] < m_v3Min[2]) pos[2] = m_v3Min[2];
		if(pos[2] > m_v3Max[2]) pos[2] = m_v3Max[2];
	}

};


//! 四面体の体積計算
static inline float calVolume(const glm::vec3 &p0, const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3)
{
	return abs((1.0/6.0)*glm::dot(glm::cross(p1-p0, p2-p0), p3-p0));
}




#endif // _RX_PBD_STRAND_H_
