/*!
  @file rx_bvh.cpp

  @brief BVH File Input

  @author Makoto Fujisawa
  @date   2021-02
*/


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------
#include "characteranimation.h"

// STLのstack : ノード情報格納に使用
#include <stack>

#include <GL/glu.h>

/*!
 * 先頭の空白(スペース，タブ)を削除
 * @param[inout] buf 処理文字列
 */
inline void DeleteHeadSpace(string& buf)
{
	size_t pos;
	while ((pos = buf.find_first_of(" 　\t")) == 0) {
		buf.erase(buf.begin());
		if (buf.empty()) break;
	}
}

/*!
 * stringからpos以降で最初の区切り文字までを抽出
 * @param[in] src 元の文字列
 * @param[out] sub 抽出文字列
 * @param[in] pos 探索開始位置
 * @param[in] sep 区切り文字
 * @return 次の抽出開始位置(","の後にスペースがあればそのスペースの後)
 */
inline int GetNextString(const string& src, string& sub, size_t pos, string sep = ",")
{
	size_t i = src.find(sep, pos);
	if (i == string::npos) {
		sub = src.substr(pos, string::npos);
		return (int)string::npos;
	}
	else {
		int cnt = 1;
		while (src[i + cnt] == ' ') {    // 区切り文字の後のスペースを消す
			cnt++;
		}
		sub = src.substr(pos, i - pos);
		return (int)(i + cnt >= src.size() ? (int)string::npos : i + cnt);
	}
}


//-----------------------------------------------------------------------------
// CharacterAnimationクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 */
CharacterAnimation::CharacterAnimation(void)
{
	m_frames = 1;
	m_dt = 0.033;
}

/*!
 * デストラクタ
 */
CharacterAnimation::~CharacterAnimation()
{
}

/*!
 * BVHファイル読み込み
 *  参考: http://www.oshita-lab.org/software/bvh/
 *        http://www.cg.ces.kyutech.ac.jp/lecture/cg2/
 * @param[in] file_name ファイル名(フルパス)
 */
bool CharacterAnimation::Read(string file_name)
{
	ifstream file;

	file.open(file_name.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		std::cout << "CharacterAnimation::Read : Invalid file specified" << std::endl;
		return false;
	}

	stack<int> node;
	int cur_idx = -1;

	size_t pos;
	string buf;
	while(!file.eof()){
		getline(file, buf);

		//// '#'以降はコメントとして無視(BVHはコメントは入れられないので必要なし)
		//if( (comment_start = buf.find('#')) != string::size_type(-1) )
		//	buf = buf.substr(0, comment_start);

		// 行頭のスペース，タブを削除
		while ((pos = buf.find_first_of(" 　\t")) == 0) {
			buf.erase(buf.begin());
			if (buf.empty()) break;
		}

		// 空行は無視
		if(buf.empty())
			continue;

		// "{"が見つかったら間接ノード番号をスタックに格納
		if((pos = buf.find_first_of("{")) != string::npos) {
			node.push(cur_idx);
		}

		// "}"が見つかったら間接ノード番号スタックをポップ
		if((pos = buf.find_first_of("}")) != string::npos){
			node.pop();
			cur_idx = node.empty() ? -1 : node.top();
		}

		// "ROOT"or"JOINT"で始まる行が見つかったら間接ノード情報の始まり
		if(buf.find("ROOT") == 0 || buf.find("JOINT") == 0){
			// 間接ノードの作成
			rxJoint joint;
			joint.parent = cur_idx;
			cur_idx = int(m_joints.size());
			joint.index = cur_idx;
			if(joint.parent != -1) {
				// 親ノードへ子ノードとして追加
				m_joints[joint.parent].children.push_back(cur_idx);
			}

			// ノード名
			if((pos = buf.find_first_of(" ")) != string::npos) {
				joint.name = buf.substr(pos+1);
			}

			m_joints.push_back(joint);
		}

		// "End"から始まる行は末端位置情報の始まりを表す
		if(buf.find("End") == 0){
			m_joints[cur_idx].is_site = true;
		}

		// "OFFSET"から始まる行は間接位置情報
		if(buf.find("OFFSET") == 0){
			string sub;		// 部分文字列
			glm::vec3 p(0.0f);	// 位置情報の一時的な格納用
			pos = GetNextString(buf, sub, 0, " "); // 最初のOFFSET部分
			pos = GetNextString(buf, sub, pos, " "); p[0] = atof(sub.c_str()); // 1番目の数値
			pos = GetNextString(buf, sub, pos, " "); p[1] = atof(sub.c_str()); // 2番目の数値
			pos = GetNextString(buf, sub, pos, " "); p[2] = atof(sub.c_str()); // 3番目の数値

			if (m_joints[cur_idx].is_site) {
				m_joints[cur_idx].site_offset = p;
			}
			else {
				m_joints[cur_idx].offset = p;
			}
		}

		// "CHANNELS"から始まる行は間接自由度情報
		if(buf.find("CHANNELS") == 0){
			string sub;		// 部分文字列
			int nchannels = 0;
			pos = GetNextString(buf, sub, 0, " "); // 最初のCHANNELS部分
			pos = GetNextString(buf, sub, pos, " "); nchannels = atoi(sub.c_str()); // 間接自由度の数
			m_joints[cur_idx].channels.resize(nchannels);
			for (int i = 0; i < nchannels; ++i) {
				rxChannel channel;
				channel.joint = cur_idx;
				channel.index = int(m_channels.size());
				pos = GetNextString(buf, sub, pos, " ");
				if (sub == "Xposition") channel.type = rxChannel::X_POS;
				else if (sub == "Yposition") channel.type = rxChannel::Y_POS;
				else if (sub == "Zposition") channel.type = rxChannel::Z_POS;
				else if (sub == "Xrotation") channel.type = rxChannel::X_ROT;
				else if (sub == "Yrotation") channel.type = rxChannel::Y_ROT;
				else if (sub == "Zrotation") channel.type = rxChannel::Z_ROT;

				m_channels.push_back(channel);
				m_joints[cur_idx].channels[i] = channel.index;
			}

		}

		// "MOTION"で始まる行が見つかったらそれ以降は動きの情報なのでループを抜ける
		if(buf.find("MOTION") == 0){
			break;
		}
	}

	string sub;	// 部分文字列

	// "MOTION"1行目:総フレーム数(モーションデータ行数)
	getline(file, buf);
	pos = GetNextString(buf, sub, 0, ":"); // 最初の"Frames:"
	pos = GetNextString(buf, sub, pos, ":"); m_frames = atoi(sub.c_str()); // 1番目の数値

	// "MOTION"2行目:時間ステップ幅(fpsの逆数)
	getline(file, buf);
	pos = GetNextString(buf, sub, 0, ":"); // 最初の"Frame Time:"
	pos = GetNextString(buf, sub, pos, ":"); m_dt = atof(sub.c_str()); // 1番目の数値

	// モーションデータ配列の初期化
	size_t nchannels = int(m_channels.size());
	m_motions.resize(nchannels*m_frames);

	// "MOTION"3行目以降:モーションデータ(各行に全チャネルの回転角度)
	for(size_t i = 0; i < m_frames; ++i){
		getline(file, buf);
		pos = 0;
		for(int j = 0; j < nchannels; ++j){
			pos = GetNextString(buf, sub, pos, " ");
			m_motions[i*nchannels+j] = atof(sub.c_str());
		}
	}

	file.close();

	return true;
}

/*!
* ノードを再帰的にたどっていって各ノード座標からAABBを計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calAABB(const int joint_idx, glm::vec3 pos, glm::vec3 &minp, glm::vec3 &maxp)
{
	if(m_joints.empty()) return;

	// 最小・最大座標の初期化(最初の1回のみ)
	if(joint_idx == 0){
		minp = glm::vec3(1000.0); maxp = glm::vec3(-1000.0);
	}

	// ジョイント座標の計算
	const rxJoint& joint = m_joints[joint_idx];
	pos += joint.offset;
	if(joint.is_site){
		pos += joint.site_offset;
	}

	// 最小・最大座標の更新
	for(int i = 0; i < 3; ++i){
		if(pos[i] < minp[i]) minp[i] = pos[i];
		if(pos[i] > maxp[i]) maxp[i] = pos[i];
	}

	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		calAABB(joint.children[i], pos, minp, maxp);
	}

	// AABBで長さ0の辺が出来ないように調整(最後の1回のみ)
	if(joint_idx == 0){
		glm::vec3 dim = maxp-minp;
		float l = glm::max(dim[0], dim[1], dim[2]);
		if(l < 1.0e-6) l = 1.0;
		for(int i = 0; i < 3; ++i){
			if(maxp[i]-minp[i] < 1.0e-6){ 
				minp[i] -= 0.05*l; maxp[i] += 0.05*l;
			}
		}
	}
}


/*!
* ノードを再帰的に辿っていって各ボーンのrest poseでのグローバル座標を計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calGlobalPos(const int joint_idx, glm::vec3 pos, vector<glm::vec3> &trans)
{
	rxJoint& joint = m_joints[joint_idx];

	if(joint.parent == -1){
		//joint_pos[0] = motion[0];
		//joint_pos[1] = motion[1];
		//joint_pos[2] = motion[2];
	}
	else{
		pos[0] += joint.offset[0];
		pos[1] += joint.offset[1];
		pos[2] += joint.offset[2];
	}

	trans[joint_idx] = pos;

	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		calGlobalPos(joint.children[i], pos, trans);
	}

}



/*!
* ノードを再帰的に辿っていって各ボーンのrest poseでのグローバル座標変換行列の計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calRestGlobalTrans(const int joint_idx, glm::mat4 mat)
{
	if(m_joints.empty()) return;

	rxJoint& joint = m_joints[joint_idx];

	glm::mat4 trs = glm::mat4(1.0f); // 単位行列で初期化
	trs[0][3] = joint.offset[0];
	trs[1][3] = joint.offset[1];
	trs[2][3] = joint.offset[2];
	mat *= trs;

	joint.global_trans = mat;

	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		calRestGlobalTrans(joint.children[i], mat);
	}

}



/*!
* ノードを再帰的に辿っていって各ボーンの変換行列を計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calGlobalAnimatedTrans(const int joint_idx, const glm::mat4 parent_mat, float *motion, float scale)
{
	if(m_joints.empty()) return;

	rxJoint& joint = m_joints[joint_idx];
	glm::mat4 trs, rot;

	// 単位行列で初期化
	trs = glm::mat4(1.0f);
	rot = glm::mat4(1.0f);

	// ルート間接は単純に平行移動，子間接の場合は親ノードからオフセットを設定
	if(joint.parent == -1){
		trs[0][3] = motion[0]*scale;
		trs[1][3] = motion[1]*scale;
		trs[2][3] = motion[2]*scale;
	}
	else{
		trs[0][3] = joint.offset[0]*scale;
		trs[1][3] = joint.offset[1]*scale;
		trs[2][3] = joint.offset[2]*scale;
	}

	// 親間接からの回転
	for(int i = 0; i < joint.channels.size(); ++i){
		const rxChannel &channel = m_channels[joint.channels[i]];
		float ang = motion[channel.index]*RX_DEGREES_TO_RADIANS;
		float c = cos(ang), s = sin(ang);
		glm::mat4 ri = glm::mat4(1.0f);
		if(channel.type == rxChannel::X_ROT){
			ri[1][1] =  c; ri[1][2] = -s;
			ri[2][1] =  s; ri[2][2] =  c;
		}
		else if(channel.type == rxChannel::Y_ROT){
			ri[0][0] =  c; ri[0][2] =  s;
			ri[2][0] = -s; ri[2][2] =  c;
		}
		else if(channel.type == rxChannel::Z_ROT){
			ri[0][0] =  c; ri[0][1] = -s;
			ri[1][0] =  s; ri[1][1] =  c;
		}
		rot = rot*ri;
	}
	glm::mat4 global_trans = parent_mat*trs*rot;

	// 変換行列の格納
	joint.global_animated_trans = global_trans;


	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		calGlobalAnimatedTrans(joint.children[i], global_trans, motion, scale);
	}

}


/*!
* 各間接のrest poseでの変換行列および指定したステップでの動きを含む変換行列を計算ｎ
* @param[in] step 現在のステップ数
* @param[in] scale 描画スケール
*/
int CharacterAnimation::Trans(int step, float scale)
{
	if(m_joints.empty()) return 0;

	float *motion = &m_motions[0];
	size_t nchannels = m_channels.size();

	glm::mat4 mat(1.0f);
	calRestGlobalTrans(0, mat);

	mat = glm::mat4(1.0f);
	calGlobalAnimatedTrans(0, mat, motion+nchannels*(step%m_frames), scale);

	// DQS用に行列をDualQuaternionに変換して格納しておく
	int n_joints = static_cast<int>(m_joints.size());
	for(int i = 0; i < n_joints; ++i){
		glm::mat4 Bj = m_joints[i].global_trans;
		glm::mat4 Wj = m_joints[i].global_animated_trans;
		glm::mat4 Tj = Wj*(glm::inverse(Bj));
		//m_joints[i].dq_animated_trans.setMat(Tj);	// HACK:
	}


	return 1;
}

/*!
* 入力された頂点座標群からボーンの重みを距離に応じて計算
* @param[in] vrts 頂点座標が格納されたベクトル
* @param[out] weights vrtsと同じ大きさの配列で各頂点の重みを格納
* @return 
*/
float CharacterAnimation::calWeight(const int joint_idx, const glm::vec3 &p, const vector<glm::vec3> &posj)
{
	int parent_idx = m_joints[joint_idx].parent;
	if(parent_idx == -1){	// 親ジョイントなし
		return 1.0;
	}
	else{	// 親間接あり
		float wmax = 0.25*glm::length(m_pos[joint_idx]+posj[parent_idx]);
		if(fabs(wmax) <= 1.0e-6) return 1.0;
		float dist = glm::length(posj[parent_idx]-p);
		return (dist >= wmax ? 1.0f : 1.0f-dist/wmax);
	}

}




/*!
* 入力された頂点座標群からボーンの重みを距離に応じて計算
* @param[in] vrts 頂点座標が格納されたベクトル
* @param[out] weights vrtsと同じ大きさの配列で各頂点の重みを格納
* @return 
*/
int CharacterAnimation::Weight(const vector<glm::vec3> &vrts, vector< map<int, double> > &weights)
{
	if(m_joints.empty()) return 0;

	// 全ての間接のグローバル座標での位置を計算
	vector<glm::vec3> posj;
	posj.resize(m_joints.size(), glm::vec3(0.0));
	calGlobalPos(0, glm::vec3(0.0), posj);
	m_pos = posj; // デバッグ用

	// 各頂点毎に最近傍の関節を求める
	int n = (int)(posj.size()), vrt_idx = 0;
	for(const glm::vec3 &v: vrts){
		// 各関節との距離の計算
		float min_dist = RX_FEQ_INF;
		int min_joint = -1, joint_idx = 0;
		for(const glm::vec3 &p: posj){
			float dist = glm::length(v-p);
			if(dist <= min_dist){
				min_dist = dist; min_joint = joint_idx;
			}
			joint_idx++;
		}

		// 各間接の中点(ボーンの中心)との距離の計算
		float min_dist_b = RX_FEQ_INF;
		int min_bone = -1, bone_idx = 0;
		for(const glm::vec3 &p: posj){
			int parent_idx = m_joints[bone_idx].parent;

			// 親ジョイントを持つ場合，そことの中点をボーンの中心座標とする
			if(parent_idx != -1){
				// この場合のボーンの動きは親ジョイントのもの
				glm::vec3 cp = 0.5f*(p+posj[parent_idx]);
				float dist = glm::length(v-cp);
				if(dist <= min_dist_b){
					min_dist_b = dist; min_bone = parent_idx;
				}
			}

			// 末端ジョイントの場合は末端との中点も調べる
			if(m_joints[bone_idx].is_site){
				glm::vec3 cp = p+0.5f*m_joints[bone_idx].site_offset;
				float dist = glm::length(v-cp);
				if(dist <= min_dist_b){
					min_dist_b = dist; min_bone = bone_idx;
				}
			}

			bone_idx++;
		}
		min_dist_b *= 1.5;	// スムージングの範囲を大きくするために中点までの距離に係数をかける

		if(min_joint != -1 && min_dist < min_dist_b){
			// 一番近いのが間接の場合は，その間接を含むボーン(親と子すべて)を対応するボーンとして追加
			int parent_idx = m_joints[min_joint].parent;
			if(parent_idx == -1){	// 親ジョイントなし(自分自身のみ)
				weights[vrt_idx].insert(pair<int, float>(min_joint, 1.0));
			}
			else{	// 親間接あり
				int nc = static_cast<int>(m_joints[parent_idx].children.size()); // 子間接の数
				float w = 1.0/(nc+1.0);	// 重みは単純に全ボーン同じになるようにする(合計すると1)

				// 親間接
				weights[vrt_idx].insert(pair<int, float>(parent_idx, w));

				// 親間接に連なる子間接(自分自身(min_joint)を含む)
				for(int i = 0; i < m_joints[parent_idx].children.size(); ++i){
					int idx = m_joints[parent_idx].children[i];
					weights[vrt_idx].insert(pair<int, float>(idx, w));
				}

			}
			
		}
		else if(min_bone != -1 && min_dist_b < min_dist){
			// 一番近いのがボーン(間接間の中点)の場合は，そのボーン(親間接)のみを対応するボーンとして追加
			weights[vrt_idx].insert(pair<int, float>(min_bone, 1.0));
		}


		vrt_idx++;
	}
	


	return 1;
}

/*!
* スケルトンの動きに合わせてメッシュ頂点を移動
* @param[in] step 現在のステップ数
* @param[inout] vrts 頂点座標が格納されたベクトル
* @param[in] weights vrtsと同じ大きさの配列で各頂点の重みを格納
*/
int CharacterAnimation::SkinningLBS(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
	if(m_joints.empty()) return 0;

	// グローバル変換行列を計算
	Trans(step, 1.0);

	// 頂点毎に変換行列を重みをかけながら適用
	int nv = (int)(vrts.size());
	for(int i = 0; i < nv; ++i){
		glm::vec4 v(vrts[i][0], vrts[i][1], vrts[i][2], 1.0);

		glm::mat4 blend_mat(0.0f);
		map<int, double>::const_iterator itr = weights[i].begin();
		for(; itr != weights[i].end(); ++itr){
			int bone = itr->first;
			float wij = static_cast<float>(itr->second);
			glm::mat4 Bj = m_joints[bone].global_trans;
			glm::mat4 Wj = m_joints[bone].global_animated_trans;
			glm::mat4 Tj = Wj*(glm::inverse(Bj));

			blend_mat += wij*Tj;
		}

		glm::vec4 new_pos = blend_mat*v;
		vrts[i] = glm::vec3(new_pos[0], new_pos[1], new_pos[2]);
	}

	return 1;
}


/*!
* スケルトンの動きに合わせてメッシュ頂点を移動
* @param[in] step 現在のステップ数
* @param[inout] vrts 頂点座標が格納されたベクトル
* @param[in] weights vrtsと同じ大きさの配列で各頂点の重みを格納
*/
int CharacterAnimation::SkinningDQS(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
	if(m_joints.empty()) return 0;

	// グローバル変換行列を計算
	Trans(step, 1.0);

	//// 頂点毎に変換DQを重みをかけながら適用
	//int nv = (int)(vrts.size());
	//for(int i = 0; i < nv; ++i){
	//	const int n_joints = weights[i].size();
	//	rxDualQuaternion dq_blend;

	//	if(n_joints == 0){
	//		dq_blend = rxDualQuaternion(rxQuaternion(0, 0, 0, 1), glm::vec3(0, 0, 0));
	//	}
	//	else{
	//		// 1つ目の間接の動きの処理
	//		map<int, double>::const_iterator itr = weights[i].begin();
	//		int bone = itr->first;
	//		float wij = itr->second;
	//		dq_blend = wij*m_joints[bone].dq_animated_trans;

	//		rxQuaternion q0;	// 1つ目の間接の回転を表す四元数
	//		q0 = m_joints[bone].dq_animated_trans.getRotation();

	//		// 2つ目以降の間接の動きの処理(1つ目と回転方向が逆にならないようにする)
	//		itr++;
	//		for(; itr != weights[i].end(); ++itr){
	//			bone = itr->first;
	//			wij = itr->second;

	//			rxDualQuaternion dq = (bone == -1) ? rxDualQuaternion(rxQuaternion(0, 0, 0, 1), glm::vec3(0, 0, 0)) : m_joints[bone].dq_animated_trans;
	//			if(dq.getRotation().Dot(q0) < 0.0){
	//				wij *= -1.0;
	//			}

	//			dq_blend += wij*dq;
	//		}

	//		glm::vec3 new_pos = dq_blend.transform(vrts[i]);
	//		vrts[i] = new_pos;
	//	}
	//	
	//}

	return 1;
}


/*!
* 動きを含めたスケルトンの描画
* @param[in] step 現在のステップ数
* @param[in] scale 描画スケール
*/
void CharacterAnimation::Draw(int step, float scale)
{
	if(m_joints.empty()) return;

	size_t nchannels = m_channels.size();
	float *motion = &m_motions[0];

	//Trans(step, scale, true);

	drawJoint(0, motion+nchannels*(step%m_frames), scale);
}

/*!
 * 間接ノードの描画(子ノードも含めて再帰的に描画)
 * @param[in] joint_idx 間接ノードインデックス
 * @param[in] motion 間接自由度毎の動き(現フレームデータの先頭ポインタ)
 * @param[in] scale 描画スケール
 */
void CharacterAnimation::drawJoint(const int joint_idx, float *motion, float scale)
{

	glPushMatrix();

	const rxJoint& joint = m_joints[joint_idx];

	// ルート間接は単純に動きのみ，子間接の場合は親ノードからオフセットを設定
	if(joint.parent == -1) glTranslatef(motion[0]*scale, motion[1]*scale, motion[2]*scale);
	else glTranslatef(joint.offset[0]*scale, joint.offset[1]*scale, joint.offset[2]*scale);

	// 親間接からの回転
	for(int i = 0; i < joint.channels.size(); ++i){
		const rxChannel &channel = m_channels[joint.channels[i]];
		float ang = motion[channel.index];
		switch(channel.type){
		case rxChannel::X_ROT: glRotatef(ang, 1.0f, 0.0f, 0.0f); break;
		case rxChannel::Y_ROT: glRotatef(ang, 0.0f, 1.0f, 0.0f); break;
		case rxChannel::Z_ROT: glRotatef(ang, 0.0f, 0.0f, 1.0f); break;
		}
	}

	glPushMatrix();

	// 間接間のボーンの描画
	if(joint.children.size() == 0){	// 子間接なし=末端(site)あり
		drawCapsule(glm::vec3(0.0), joint.site_offset*scale);
	}
	else{	// 子間接が1個以上
		for(int i = 0; i < joint.children.size(); ++i){
			const rxJoint& child_joint = m_joints[joint.children[i]];
			drawCapsule(glm::vec3(0.0), child_joint.offset*scale);
		}
	}
	glPopMatrix();

	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		drawJoint(joint.children[i], motion, scale);
	}

	glPopMatrix();
}

/*!
 * カプセル描画(円筒の両端に半球をつけた形)
 * @param[in] pos0,pos1 両端の位置
 */
void CharacterAnimation::drawCapsule(glm::vec3 pos0, glm::vec3 pos1)
{
	GLUquadricObj* qobj;
	qobj = gluNewQuadric();

	glm::vec3 dir = pos1 - pos0;
	float len = glm::length(dir);
	if(len < 0.0001){
		dir = glm::vec3(0.0, 0.0, 1.0); len = 1.0;
	}
	else{
		dir /= len;
	}

	// 上方向の設定(y軸方向)
	glm::vec3 up(0.0, 1.0, 0.0);

	// ボーンの方向をz軸としてx軸方向を算出
	glm::vec3 side = glm::cross(up, dir);
	float side_len = glm::length(side);
	glm::normalize(side);
	if(side_len < 0.0001){	// dir==upだった場合の対策
		up = glm::vec3(1.0, 0.0, 0.0);
		side = glm::normalize(glm::cross(up, dir));
	}

	// x,z軸からy軸方向を算出
	up = glm::normalize(glm::cross(dir, side));

	// 回転行列の設定
	GLfloat m[16] = { side[0], side[1], side[2], 0.0,
					   up[0],   up[1],  up[2],    0.0,
					   dir[0],  dir[1], dir[2],   0.0,
					   0.0,     0.0,    0.0,      1.0 };
	
	glPushMatrix();

	glTranslatef(pos0[0], pos0[1], pos0[2]);
	glMultMatrixf(m);

	GLfloat rad = 0.02;
	GLuint slices = 8;
	GLuint stacks = 3;

	gluQuadricDrawStyle(qobj, GLU_FILL);
	gluQuadricNormals(qobj, GLU_SMOOTH);
	gluCylinder(qobj, rad, rad, len, slices, stacks);

	glPushMatrix();
	GLUquadricObj* qobj_s1 = gluNewQuadric();
	gluSphere(qobj_s1, 1.3*rad, slices, slices);
	//glutSolidSphere(1.3*rad, slices, slices);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0, 0.0, len);
	GLUquadricObj* qobj_s2 = gluNewQuadric();
	gluSphere(qobj_s2, 1.3*rad, slices, slices);
	//glutSolidSphere(1.3*rad, slices, slices);
	glPopMatrix();

	glPopMatrix();

}
