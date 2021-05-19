/*!
  @file rx_bvh.h
	
  @brief BVH File Input

  @author Makoto Fujisawa
  @date   2021-02
*/

#ifndef _RX_BVH_H_
#define _RX_BVH_H_


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------
#include "utils.h"
//#include "dualquaternion.h"

using namespace std;



// 間接自由度(BVHではCHANNEL)
// これに対応する数のモーションが記述されているのでJointとは別で定義する
class rxChannel
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
class rxJoint
{
public:
	string name;			//!< 間接名
	int index;				//!< 間接番号
	int parent;				//!< 親ノードインデックス(-1で親なし)
	vector<int> children;	//!< 子ノードインデックス
	glm::vec3 offset;			//!< 間接位置(親ノードからのオフセット)
	glm::mat4 global_trans;	//!< 間接位置のグローバル座標(ルート間接が原点)への変換行列(rest pose)
	glm::mat4 global_animated_trans;	//!< 間接の移動を含む変換行列

	//rxDualQuaternion dq_trans;
	//rxDualQuaternion dq_animated_trans;

	bool is_site;			//!< 末端ノードかどうかのフラグ
	glm::vec3 site_offset;		//!< 末端位置(間接位置からのオフセット)

	vector<int> channels;	//!< 間接自由度(これ毎にモーションが定義される)

	rxJoint()
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
	vector<rxJoint> m_joints;		//!< 間接情報
	vector<rxChannel> m_channels;	//!< 各間接自由度情報
	vector<float> m_motions;		//!< 各間接自由度毎の動き(全フレーム分)

	int m_frames;					//!< アニメーションフレーム数
	float m_dt;					//!< アニメーションタイムステップ幅

	map<int, glm::mat4> m_trans;		//!< 各間接毎のグローバル変換行列(4x4)
	vector<glm::mat4> m_vtrans;

public:
	// デバッグ用
	vector<glm::vec3> m_pos;

public:
	//! コンストラクタ
	CharacterAnimation();
	//! デストラクタ
	~CharacterAnimation();

	/*!
	 * BVHファイル読み込み
	 * @param[in] file_name ファイル名(フルパス)
	 */
	bool Read(string file_name);

	void Draw(int step, float scale = 1.0);
	void AABB(glm::vec3 &minp, glm::vec3 &maxp){ calAABB(0, glm::vec3(0.0), minp, maxp); }
	float Dt(void) const { return m_dt; }

	int Trans(int step, float scale);
	int Weight(const vector<glm::vec3> &vrts, vector< map<int, double> > &weights);

	int SkinningLBS(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights);
	int SkinningDQS(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights);

private:
	void drawJoint(const int joint_idx, float *motion, float scale = 1.0);
	void drawCapsule(glm::vec3 pos0, glm::vec3 pos1);

	void calAABB(const int joint_idx, glm::vec3 pos, glm::vec3 &minp, glm::vec3 &maxp);
	void calGlobalPos(const int joint_idx, glm::vec3 pos, vector<glm::vec3> &trans);
	void calGlobalAnimatedTrans(const int joint_idx, const glm::mat4 mat, float *motion, float scale);
	void calRestGlobalTrans(const int joint_idx, glm::mat4 mat);
	
	float calWeight(const int joint_idx, const glm::vec3 &pos, const vector<glm::vec3> &posj);
};



#endif // _RX_BVH_H_
