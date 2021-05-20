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

// 四元数
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>


/*!
 * 先頭の空白(スペース，タブ)を削除
 * @param[inout] buf 処理文字列
 */
inline void DeleteHeadSpace(string& buf)
{
	size_t pos;
	while((pos = buf.find_first_of(" 　\t")) == 0) {
		buf.erase(buf.begin());
		if(buf.empty()) break;
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
	m_skinning = 0;
	m_rest_pose = false;

	m_sphere.vao = MakeSphereVAO(m_sphere.nvrts, m_sphere.ntris, 0.5, 32, 16);
	m_cylinder.vao = MakeCylinderVAO(m_cylinder.nvrts, m_cylinder.ntris, 0.5, 1.0, 32);
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
			Joint joint;
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
				Channel channel;
				channel.joint = cur_idx;
				channel.index = int(m_channels.size());
				pos = GetNextString(buf, sub, pos, " ");
				if (sub == "Xposition") channel.type = Channel::X_POS;
				else if (sub == "Yposition") channel.type = Channel::Y_POS;
				else if (sub == "Zposition") channel.type = Channel::Z_POS;
				else if (sub == "Xrotation") channel.type = Channel::X_ROT;
				else if (sub == "Yrotation") channel.type = Channel::Y_ROT;
				else if (sub == "Zrotation") channel.type = Channel::Z_ROT;

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
* 動きを含めたスケルトンの描画
* @param[in] step 現在のステップ数
* @param[in] scale 描画スケール
*/
void CharacterAnimation::Draw(int step, float scale)
{
	if(m_joints.empty()) return;

	size_t nchannels = m_channels.size();
	float *motion = &m_motions[0];

	// ルート関節から順番に全ての間接のグローバル変換行列(回転+平行移動)を計算
	calTransMatrices(step, scale);

	// 各関節を再帰的に呼び出しながら描画
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

	const Joint& joint = m_joints[joint_idx];

	//// 間接のオフセット,回転(motion)情報を直接使う場合
	//// ルート間接は単純に動きのみ，子間接の場合は親ノードからオフセットを設定
	//if(joint.parent == -1) glTranslatef(motion[0]*scale, motion[1]*scale, motion[2]*scale);
	//else glTranslatef(joint.offset[0]*scale, joint.offset[1]*scale, joint.offset[2]*scale);


	//// 親間接からの回転
	//for(int i = 0; i < joint.channels.size(); ++i){
	//	const Channel &channel = m_channels[joint.channels[i]];
	//	float ang = motion[channel.index];
	//	switch(channel.type){
	//	case Channel::X_ROT: glRotatef(ang, 1.0f, 0.0f, 0.0f); break;
	//	case Channel::Y_ROT: glRotatef(ang, 0.0f, 1.0f, 0.0f); break;
	//	case Channel::Z_ROT: glRotatef(ang, 0.0f, 0.0f, 1.0f); break;
	//	}
	//}

	glPushMatrix();

	// 各間接の情報から算出された変換行列を使う場合(CalTransMatrices関数を先に呼び出す必要あり)
	if(m_rest_pose) glMultMatrixf(glm::value_ptr(joint.global_trans));
	else glMultMatrixf(glm::value_ptr(joint.global_animated_trans));

	// 間接間のボーンの描画
	if(joint.children.size() == 0){	// 子間接なし=末端(site)あり
		drawCapsule(glm::vec3(0.0), joint.site_offset*scale);
	}
	else{	// 子間接が1個以上
		for(int i = 0; i < joint.children.size(); ++i){
			const Joint& child_joint = m_joints[joint.children[i]];
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
 * 2つのベクトル間の回転を表す四元数
 *  - http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-17-quaternions/
 * @param[in] src,dst 2つのベクトル
 * @return 四元数
 */
glm::quat calRotationBetweenVectors(glm::vec3 src, glm::vec3 dst)
{
	src = glm::normalize(src);
	dst = glm::normalize(dst);

	float c = glm::dot(src, dst);
	glm::vec3 axis;

	if(c < -1.0f + 1e-6){
		// 2つのベクトルが反対方向を向いている場合
		axis = glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), src);
		if(glm::length2(axis) < 0.01) // 平行の場合
			axis = glm::cross(glm::vec3(1.0f, 0.0f, 0.0f), src);

		axis = glm::normalize(axis);
		return glm::quat(glm::radians(180.0f), axis);
	}

	axis = glm::cross(src, dst);

	float s = sqrt((1+c)*2);
	float inv_s = 1 / s;

	return glm::quat(s*0.5f, axis.x*inv_s, axis.y*inv_s, axis.z*inv_s);
}

/*!
 * カプセル描画(円筒の両端に半球をつけた形)
 * @param[in] pos0,pos1 両端の位置
 */
void CharacterAnimation::drawCapsule(glm::vec3 pos0, glm::vec3 pos1)
{
	if(!m_cylinder.vao || !m_sphere.vao) return;

	// 2端点間のベクトル
	glm::vec3 dir = pos1 - pos0;
	float len = glm::length(dir);
	if(len < 0.0001){
		return;
	}
	else{
		dir /= len;
	}

	// 円筒のデフォルトの向き(z軸方向)と2端点間ベクトルの間の回転を洗わす四元数を計算
	glm::quat q = calRotationBetweenVectors(glm::vec3(0.0f, 0.0f, 1.0f), dir);

	glPushMatrix();

	// 平行移動＆回転
	glTranslatef(pos0[0], pos0[1], pos0[2]);
	glMultMatrixf(glm::value_ptr(glm::mat4_cast(q)));

	GLfloat rad = 0.025;

	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glColor3d(0.1, 0.5, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// 円筒描画
	glPushMatrix();
	glTranslatef(0.0f, 0.0f, 0.5f*len);
	glScalef(2*rad, 2*rad, len);
	glBindVertexArray(m_cylinder.vao);
	glDrawElements(GL_TRIANGLES, m_cylinder.ntris*3, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	glPopMatrix();

	// 関節部分に球体を描画
	glPushMatrix();
	glScalef(2.6*rad, 2.6*rad, 2.6*rad);
	glBindVertexArray(m_sphere.vao);
	glDrawElements(GL_TRIANGLES, m_sphere.ntris*3, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0, 0.0, len);
	glScalef(2.6*rad, 2.6*rad, 2.6*rad);
	glBindVertexArray(m_sphere.vao);
	glDrawElements(GL_TRIANGLES, m_sphere.ntris*3, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	glPopMatrix();

	glPopMatrix();

}


/*!
* スケルトンのAABBを計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::AABB(glm::vec3 &minp, glm::vec3 &maxp, bool rest)
{
	// 各関節の変換行列を計算
	calTransMatrices(0, 1.0f);

	minp = glm::vec3(1000.0); maxp = glm::vec3(-1000.0);

	// 全間接についてグローバル座標での位置を求めてAABBを算出
	vector<Joint>::iterator itr = m_joints.begin();
	for(; itr != m_joints.end(); ++itr){
		const Joint& joint = *itr;
		glm::vec4 pos(0.0f, 0.0f, 0.0f, 1.0f);
		if(rest){
			pos = joint.global_trans*pos;
		}
		else{
			pos = joint.global_animated_trans*pos;
		}

		// 最小・最大座標の更新
		for(int i = 0; i < 3; ++i){
			if(pos[i] < minp[i]) minp[i] = pos[i];
			if(pos[i] > maxp[i]) maxp[i] = pos[i];
		}

		// 端点
		if(joint.is_site){
			// 最小・最大座標の更新
			for(int i = 0; i < 3; ++i){
				pos[i] += joint.site_offset[i];
				if(pos[i] < minp[i]) minp[i] = pos[i];
				if(pos[i] > maxp[i]) maxp[i] = pos[i];
			}
		}

	}

	// AABBで長さ0の辺が出来ないように調整(最後の1回のみ)
	glm::vec3 dim = maxp-minp;
	float l = glm::max(dim[0], dim[1], dim[2]);
	if(l < 1.0e-6) l = 1.0;
	for(int i = 0; i < 3; ++i){
		if(maxp[i]-minp[i] < 1.0e-6){
			minp[i] -= 0.05*l; maxp[i] += 0.05*l;
		}
	}
}


/*!
* ノードを再帰的に辿っていって各ボーンのrest poseでのグローバル座標を計算
* - 重み計算時に使う
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calGlobalPos(const int joint_idx, glm::vec3 pos, vector<glm::vec3> &trans)
{
	Joint& joint = m_joints[joint_idx];

	if(joint.parent != -1){
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
		float wmax = 0.25*glm::length(posj[joint_idx]+posj[parent_idx]);
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
* ノードを再帰的に辿っていって各ボーンのrest poseでのグローバル座標変換行列の計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calRestGlobalTrans(const int joint_idx, glm::mat4 mat, float scale)
{
	if(m_joints.empty()) return;

	Joint& joint = m_joints[joint_idx];

	glm::mat4 trs = glm::translate(glm::mat4(), joint.offset*scale); // 平行移動行列
	mat *= trs;

	joint.global_trans = mat;

	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		calRestGlobalTrans(joint.children[i], mat, scale);
	}

}

/*!
* ノードを再帰的に辿っていって各ボーンの変換行列を計算
* @param[in] joint_idx 間接ノードインデックス
*/
void CharacterAnimation::calAnimatedGlobalTrans(const int joint_idx, const glm::mat4 parent_mat, float *motion, float scale)
{
	if(m_joints.empty()) return;

	Joint& joint = m_joints[joint_idx];

	// ルート間接は単純に平行移動，子間接の場合は親ノードからオフセットを設定
	glm::vec3 offset(0.0f);
	if(joint.parent == -1){
		offset = glm::vec3(motion[0], motion[1], motion[2]);
	}
	else{
		offset = joint.offset;
	}
	glm::mat4 trs = glm::translate(glm::mat4(), offset*scale); // 平行移動行列

	// 親間接からの回転
	glm::mat4 rot(1.0f);	// 回転行列(単位行列で初期化)
	for(int i = 0; i < joint.channels.size(); ++i){
		const Channel &channel = m_channels[joint.channels[i]];
		float ang = motion[channel.index];// 角度は[deg](glm::rotateに渡す角度も[deg]なので変換なし)
		glm::mat4 ri = glm::mat4(1.0f);
		if(channel.type == Channel::X_ROT){
			ri = glm::rotate(ang, glm::vec3(1.0f, 0.0f, 0.0f));
		}
		else if(channel.type == Channel::Y_ROT){
			ri = glm::rotate(ang, glm::vec3(0.0f, 1.0f, 0.0f));
		}
		else if(channel.type == Channel::Z_ROT){
			ri = glm::rotate(ang, glm::vec3(0.0f, 0.0f, 1.0f));
		}
		rot = rot*ri;
	}
	glm::mat4 global_trans = parent_mat*trs*rot;

	// 変換行列の格納
	joint.global_animated_trans = global_trans;

	// 子間接を再帰呼び出し
	for(int i = 0; i < joint.children.size(); ++i){
		calAnimatedGlobalTrans(joint.children[i], global_trans, motion, scale);
	}

}

/*!
* 各間接のrest poseでの変換行列および指定したステップでの動きを含む変換行列を計算
* @param[in] step 現在のステップ数
* @param[in] scale 描画スケール
*/
int CharacterAnimation::calTransMatrices(int step, float scale)
{
	if(m_joints.empty()) return 0;

	float *motion = &m_motions[0];
	size_t nchannels = m_channels.size();

	glm::mat4 mat(1.0f);
	calRestGlobalTrans(0, mat, scale);

	mat = glm::mat4(1.0f);
	calAnimatedGlobalTrans(0, mat, motion+nchannels*(step%m_frames), scale);

	// DQS用に行列をDualQuaternionに変換して格納しておく
	int n_joints = static_cast<int>(m_joints.size());
	for(int i = 0; i < n_joints; ++i){
		glm::mat4 Bj = m_joints[i].global_trans;
		glm::mat4 Wj = m_joints[i].global_animated_trans;
		glm::mat4 Tj = Wj*(glm::inverse(Bj));
		m_joints[i].dq_animated_trans.setMat(Tj);
	}


	return 1;
}


/*!
* スケルトンの動きに合わせてメッシュ頂点を移動
* - LBSによるスキニング
* @param[in] step 現在のステップ数
* @param[inout] vrts 頂点座標が格納されたベクトル
* @param[in] weights vrtsと同じ大きさの配列で各頂点の重みを格納
*/
int CharacterAnimation::skinningLBS(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
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
* - DQSによるスキニング
* @param[in] step 現在のステップ数
* @param[inout] vrts 頂点座標が格納されたベクトル
* @param[in] weights vrtsと同じ大きさの配列で各頂点の重みを格納
*/
int CharacterAnimation::skinningDQS(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
	// 頂点毎に変換DQを重みをかけながら適用
	int nv = (int)(vrts.size());
	for(int i = 0; i < nv; ++i){
		const int n_joints = static_cast<int>(weights[i].size());

		//DualQuaternion dq_blend;
		//if(n_joints == 0){
		//	dq_blend = DualQuaternion(glm::quat(1, 0, 0, 0), glm::vec3(0, 0, 0));
		//}
		//else{
		//	// 1つ目の間接の動きの処理
		//	map<int, double>::const_iterator itr = weights[i].begin();
		//	int bone = itr->first;
		//	float wij = itr->second;
		//	dq_blend = wij*m_joints[bone].dq_animated_trans;

		//	glm::quat q0;	// 1つ目の間接の回転を表す四元数
		//	q0 = m_joints[bone].dq_animated_trans.getRotation();

		//	// 2つ目以降の間接の動きの処理(1つ目と回転方向が逆にならないようにする)
		//	itr++;
		//	for(; itr != weights[i].end(); ++itr){
		//		bone = itr->first;
		//		wij = itr->second;

		//		DualQuaternion dq = (bone == -1) ? DualQuaternion(glm::quat(1, 0, 0, 0), glm::vec3(0, 0, 0)) : m_joints[bone].dq_animated_trans;
		//		if(glm::dot(dq.getRotation(), q0) < 0.0){
		//			wij *= -1.0;
		//		}

		//		dq_blend += wij*dq;
		//	}

		//	glm::vec3 new_pos = dq_blend.transform(vrts[i]);
		//	vrts[i] = new_pos;
		//}

		// 行列をDual Quaternionにその都度変換するだけのシンプルな例
		DualQuaternion dq_blend(glm::quat(0, 0, 0, 0), glm::quat(0, 0, 0, 0));
		map<int, double>::const_iterator itr = weights[i].begin();
		for(; itr != weights[i].end(); ++itr){
			int bone = itr->first;
			float wij = static_cast<float>(itr->second);
			glm::mat4 Bj = m_joints[bone].global_trans;
			glm::mat4 Wj = m_joints[bone].global_animated_trans;
			glm::mat4 Tj = Wj*(glm::inverse(Bj));

			DualQuaternion dq;
			dq.setMat(Tj);
			dq_blend += wij*dq;
		}

		glm::vec3 new_pos = dq_blend.transform(vrts[i]);
		vrts[i] = new_pos;
	}

	return 1;
}

/*!
* スケルトンの動きに合わせてメッシュ頂点を移動
* @param[in] step 現在のステップ数
* @param[inout] vrts 頂点座標が格納されたベクトル
* @param[in] weights vrtsと同じ大きさの配列で各頂点の重みを格納
*/
int CharacterAnimation::Skinning(int step, vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
	if(m_joints.empty()) return 0;

	// グローバル変換行列を計算
	calTransMatrices(step, 1.0);

	// スキニング
	switch(m_skinning){
	default:
	case 0:
		return skinningLBS(step, vrts, weights);
	case 1:
		return skinningDQS(step, vrts, weights);
	}
}
