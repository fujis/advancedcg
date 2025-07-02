/*!
  @file characteranimation.h
	
  @brief BVHファイルに基づくキャラクターアニメーション

  @author Makoto Fujisawa
  @date   2021-02
*/

#ifndef _RX_BVH_H_
#define _RX_BVH_H_


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------
#include "utils.h"
#include "dualquaternion.h"
#include "rx_mesh.h"

using namespace std;



// 間接自由度(BVHではCHANNEL)
// これに対応する数のモーションが記述されているのでJointとは別で定義する
class Channel
{
public:
	//! 間接自由度
	enum
	{
		X_ROT, Y_ROT, Z_ROT,
		X_POS, Y_POS, Z_POS,
	};

	int index;	//!< チャンネル番号
	int type;	//!< 間接自由度のタイプ
	int joint;	//!< 対応する間接ノード番号
};



//! 間接ノード格納用
class Joint
{
public:
	string name;			//!< 間接名
	int index;				//!< 間接番号
	int parent;				//!< 親ノードインデックス(-1で親なし)
	vector<int> children;	//!< 子ノードインデックス
	glm::vec3 offset;		//!< 間接位置(親ノードからのオフセット)
	glm::mat4 B;			//!< 間接位置のグローバル座標(ルート間接が原点)への変換行列(rest pose)
	glm::mat4 W;			//!< 間接の移動を含む変換行列

	bool is_site;			//!< 末端ノードかどうかのフラグ
	glm::vec3 site_offset;		//!< 末端位置(間接位置からのオフセット)

	vector<int> channels;	//!< 間接自由度(これ毎にモーションが定義される)

	Joint()
	{
		index = parent = -1;
		offset = site_offset = glm::vec3(0.0);
		is_site = false;
	}
};


//-----------------------------------------------------------------------------
// CharacterAnimationクラス
//  - スケルトンによるアニメーション
//  - BVH形式からのスケルトンデータの読み込み
//  - メッシュデータのスキニング
//    (メッシュデータ自体は頂点列を引数として取るだけでこのクラス内では管理しない)
//-----------------------------------------------------------------------------
class CharacterAnimation
{
	// VAOと頂点数,ポリゴン数をまとめて管理するための構造体:スケルトン描画用
	struct Primitive
	{
		int nvrts, ntris;
		GLuint vao;
	};

	vector<Joint> m_joints;		//!< 間接情報
	vector<Channel> m_channels;	//!< 各間接自由度情報
	vector<float> m_motions;	//!< 各間接自由度毎の動き(全フレーム分)

	int m_frames;				//!< アニメーションフレーム数
	float m_dt;					//!< アニメーションタイムステップ幅

	// スケルトン描画用メッシュ
	rxPolygons m_sphere;		//!< 関節部分の球体
	rxPolygons m_cylinder;		//!< ボーンを描画するための円筒

public:
	int m_skinning;				//!< スキニング方法
	bool m_rest_pose;			//!< 初期ポーズ描画フラグ

public:
	//! コンストラクタ
	CharacterAnimation();
	//! デストラクタ
	~CharacterAnimation();

	// BVHファイル読み込み
	bool Read(string file_name);
	void CheckData(void);

	// 描画＆情報取得
	void Draw(int step, float scale = 1.0);
	void AABB(glm::vec3 &minp, glm::vec3 &maxp, bool rest = false);// { calAABB(0, glm::vec3(0.0), minp, maxp); }
	float Dt(void) const { return m_dt; }

	// スキニング用重みの計算
	int Weight(const vector<glm::vec3> &vrts, vector< map<int, double> > &weights);

	// スキニング
	int Skinning(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights);

private:
	// 描画
	void drawJoint(const int joint_idx, float *motion, float scale = 1.0);
	void drawCapsule(glm::vec3 pos0, glm::vec3 pos1);
	void drawPolygon(rxPolygons &poly);


	// 各間接姿勢への変換行列の計算
	int calTransMatrices(int step, float scale);
	void calRestGlobalTrans(const int joint_idx, glm::mat4 mat, float scale);
	void calAnimatedGlobalTrans(const int joint_idx, const glm::mat4 mat, float *motion, float scale);
	
	// メッシュ頂点の重み計算
	float calWeight(const int joint_idx, const glm::vec3 &pos, const vector<glm::vec3> &posj);
	void calGlobalPos(const int joint_idx, glm::vec3 pos, vector<glm::vec3> &trans);

	// スキニング
	int skinningLBS(vector<glm::vec3> &vrts, const vector< map<int, double> > &weights);
	int skinningDQS(vector<glm::vec3> &vrts, const vector< map<int, double> > &weights);
};


/*!
 * 文字列中の改行記号削除
 * @param[inout] 文字列
 */
static inline void DeleteEOL(std::string &str)
{
	const char CR = '\r';
	const char LF = '\n';
	std::string str1;
	for(std::string::const_iterator itr = str.begin(); itr != str.end(); ++itr){
		if(*itr != CR && *itr != LF) str1 += *itr;
	}
	str = str1;
}


#endif // _RX_BVH_H_
