/*!
  @file characteranimation.cpp

  @brief BVHファイルに基づくキャラクターアニメーション

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


//-----------------------------------------------------------------------------
// 課題用関数
//-----------------------------------------------------------------------------
/*!
* スケルトンの動きに合わせてメッシュ頂点を移動
* - LBSによるスキニング
* @param[in] step 現在のステップ数
* @param[inout] vrts 頂点座標が格納されたベクトル
* @param[in] weights vrtsと同じ大きさの配列で各頂点の重みを格納
*/
int CharacterAnimation::skinningLBS(vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
	// 頂点毎に変換行列を重みをかけながら適用
	int nv = (int)(vrts.size());
	for(int i = 0; i < nv; ++i){
		const int n_joints = static_cast<int>(weights[i].size());	// 頂点vに対応するジョイント数
		glm::vec4 v(vrts[i][0], vrts[i][1], vrts[i][2], 1.0);	// 更新前(オリジナル)のスキンメッシュ頂点位置
		glm::vec4 v_new = v;	// 更新後のスキンメッシュ頂点位置

		// TODO:この部分にLBSによる頂点位置の計算を書く
		// ・変数v_newにスキンメッシュ頂点vのLBSによる更新後の位置を格納する
		// ・元の頂点座標は3次元ベクトル(glm::vec3)として格納されているが，
		//   4x4行列との演算のために4次元ベクトル(glm::vec4)に変換していることに注意
		// ・頂点vに対応するジョイント番号とその重みは以下のようにして取得できる
		//	map<int, double>::const_iterator itr = weights[i].begin();
		//	for(; itr != weights[i].end(); ++itr){
		//		int j = itr->first;	// ジョイント番号
		//		float wij = static_cast<float>(itr->second);	// 重み
		//		// ここにジョイントjに関する処理を書く
		//
		//	}
		// ・ジョイントjのrest(bind) poseでのワールド変換行列(スライドp28のBj)は
		//    m_joints[j].B
		//   に格納されており，回転を含むワールド変換行列(スライドp28のWj)は
		//    m_joints[j].W
		//   に格納されている(どちらもglm::mat4型)
		// [glmでのベクトル・行列演算について]
		// ・glm::mat4 M(0.0f) と定義時に引数に0を指定すると0で初期化，
		//   glm::mat4 M(1.0f) と指定すると単位行列で初期化される(LBSでは行列WB^-1を足していくので単位行列ではなく...)
		// ・glmでの行列Mとベクトルvの掛け算は単純に M*v でよい(結果のベクトルが返ってくる)
		// ・glmでの逆行列計算は glm::inverse() を使うと良い

		// ----課題ここから----



		// ----課題ここまで----

		vrts[i] = glm::vec3(v_new[0], v_new[1], v_new[2]);
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
int CharacterAnimation::skinningDQS(vector<glm::vec3> &vrts, const vector< map<int, double> > &weights)
{
	// 頂点毎に変換DQを重みをかけながら適用
	int nv = (int)(vrts.size());
	for(int i = 0; i < nv; ++i){
		const int n_joints = static_cast<int>(weights[i].size());	// 頂点vに対応するジョイント数
		glm::vec3 v = vrts[i];	// 更新前(オリジナル)のスキンメッシュ頂点位置
		glm::vec3 v_new = v;	// 更新後のスキンメッシュ頂点位置

		// TODO:この部分にDQSによる頂点位置の計算を書く
		// ・変数v_newにスキンメッシュ頂点vのDQSによる更新後の位置を格納する
		// ・LBSと違って4x4行列との掛け算はない(全て四元数との演算)ため，
		//   こちらではglm::vec3として頂点座標を取り出していることに注意
		// ・頂点vに対応するジョイント番号とその重みは以下のようにして取得できる
		//	map<int, double>::const_iterator itr = weights[i].begin();
		//	for(; itr != weights[i].end(); ++itr){
		//		int j = itr->first;	// ジョイント番号
		//		float wij = static_cast<float>(itr->second);	// 重み
		//		// ここにジョイントjに関する処理を書く
		//
		//	}
		// ・ジョイントjのrest(bind) poseでのワールド変換行列(スライドp28のBj)は
		//    m_joints[j].B
		//   に格納されており，回転を含むワールド変換行列(スライドp28のWj)は
		//    m_joints[j].W
		//   に格納されている(どちらもglm::mat4型)
		// ・四元数や行列については授業スライドの最後の方の説明を参照すること
		// ・Dual Quaternionを扱うための，DualQuaternion型を定義してある(dualquaternion.h)ので自由に使ってください．
		//   基本的な演算は定義してありますが，頂点vのDual Quaternionによる変換は定義されていないので自分でここに実装してください．
		// ・注意としてDual Quaternionを重みを掛けながら足し合わせていくという演算をする場合，Dual Quaternionの初期化が必要になります．
 		//   こちらで用意したDualQuaternionクラスはデフォルトではreal partを単位四元数(1,0,0,0)，dual partを(0,0,0,0)で初期化しますが，
		//   Dual Quaterionを足し合わせて行く場合は回転部分を表すreal partも(0,0,0,0)で初期化するか，
		//   1つめの重み&位置姿勢を初期値として代入して，2つめ以降の重み&位置姿勢を足し合わせていくといった方法を取る必要があります．

		// ----課題ここから----



		// ----課題ここまで----

		vrts[i] = glm::vec3(v_new[0], v_new[1], v_new[2]);
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
		return skinningLBS(vrts, weights);
	case 1:
		return skinningDQS(vrts, weights);
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

	// 本当はVAOで描画をしていたがmac上でimgui_opengl2+vaoが上手く動かないため，メッシュ生成に切り替え
	int nvrts, ntris;
	vector<int> tris;

	// 球形状を表すポリゴンを作ってVBO登録(関節描画用)
	MakeSphere(nvrts, m_sphere.vertices, m_sphere.normals, ntris, tris, 0.5, 32, 16);
	rxFace f;
	f.resize(3);
	for(int i = 0; i < ntris; ++i){
		f[0] = tris[3*i]; f[1] = tris[3*i+1]; f[2] = tris[3*i+2];
		m_sphere.faces.push_back(f);
	}
	tris.clear();

	// 円筒形状を表すポリゴンを作ってVBO登録(ボーン描画用)
	MakeCylinder(nvrts, m_cylinder.vertices, m_cylinder.normals, ntris, tris, 0.5, 1.0, 32, 16);
	for(int i = 0; i < ntris; ++i){
		f[0] = tris[3*i]; f[1] = tris[3*i+1]; f[2] = tris[3*i+2];
		m_cylinder.faces.push_back(f);
	}

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

		// 改行記号(EOL)削除(macだと改行記号が違うのでこれをしていないと最後の文字の処理時に問題が発生することがある)
		DeleteEOL(buf);

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
 * 読み込んだBVHデータの画面出力(デバッグ用)
 */
void CharacterAnimation::CheckData(void)
{
	for(const Joint &joint : m_joints){
		cout << "[" << joint.name << "]" << endl;
		cout << "  index : " << joint.index << endl;
		cout << "  offset : " << glm::to_string(joint.offset) << endl;
		if(joint.is_site){
			cout << "  site_offset : " << glm::to_string(joint.site_offset) << endl;
		}

		cout << "  channels : ";
		for(int c : joint.channels){
			const Channel &channel = m_channels[c];
			if(channel.type == Channel::X_POS) cout << "Xposition ";
			else if(channel.type == Channel::Y_POS) cout << "Yposition ";
			else if(channel.type == Channel::Z_POS) cout << "Zposition ";
			else if(channel.type == Channel::X_ROT) cout << "Xrotation ";
			else if(channel.type == Channel::Y_ROT) cout << "Yrotation ";
			else if(channel.type == Channel::Z_ROT) cout << "Zrotation ";
		}
		cout << endl;

		cout << "  parent : " << joint.parent << endl;
		if(!joint.children.empty()){
			cout << "  childrens : ";
			for(int child : joint.children){
				cout << child << " ";
			}
			cout << endl;
		}

	}

	cout << "[motion]" << endl;
	int nchannels = static_cast<int>(m_channels.size());
	int nmotions = static_cast<int>(m_motions.size());
	for(int i = 0; i < nmotions; i += nchannels){
		cout << "  ";
		for(int j = 0; j < nchannels; ++j){
			cout << m_motions[i+j] << " ";
		}
		cout << endl;
	}
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

	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

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
	if(m_rest_pose) glMultMatrixf(glm::value_ptr(joint.B));
	else glMultMatrixf(glm::value_ptr(joint.W));

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
	if(m_cylinder.vertices.empty() || m_sphere.vertices.empty()) return;

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

	glColor3d(0.1, 0.5, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// 円筒描画
	glPushMatrix();
	glTranslatef(0.0f, 0.0f, 0.5f*len);
	glScalef(2*rad, 2*rad, len);
	drawPolygon(m_cylinder);
	// VAOを使う場合
	//glBindVertexArray(m_cylinder.vao);
	//glDrawElements(GL_TRIANGLES, m_cylinder.ntris*3, GL_UNSIGNED_INT, 0);
	//glBindVertexArray(0);
	glPopMatrix();

	// 関節部分に球体を描画
	glPushMatrix();
	glScalef(2.6*rad, 2.6*rad, 2.6*rad);
	drawPolygon(m_sphere);
	//glBindVertexArray(m_sphere.vao);
	//glDrawElements(GL_TRIANGLES, m_sphere.ntris*3, GL_UNSIGNED_INT, 0);
	//glBindVertexArray(0);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0, 0.0, len);
	glScalef(2.6*rad, 2.6*rad, 2.6*rad);
	drawPolygon(m_sphere);
	//glBindVertexArray(m_sphere.vao);
	//glDrawElements(GL_TRIANGLES, m_sphere.ntris*3, GL_UNSIGNED_INT, 0);
	//glBindVertexArray(0);
	glPopMatrix();

	glPopMatrix();

}


/*!
 * ポリゴンの描画
 * @param[in] polys ポリゴンデータ
 * @param[in] draw 描画フラグ(下位ビットから頂点,エッジ,面,法線 - 1,2,4,8)
 */
void CharacterAnimation::drawPolygon(rxPolygons &poly)
{
	// 頂点数とポリゴン数
	int vn = (int)poly.vertices.size();
	int pn = (int)poly.faces.size();

	// すべてのポリゴンを描画
	for(int i = 0; i < pn; ++i){
		const rxFace *face = &(poly.faces[i]);
		int n = (int)face->vert_idx.size();

		// ポリゴン描画
		glBegin(GL_POLYGON);
		for(int j = 0; j < n; ++j){
			int idx = face->vert_idx[j];
			glNormal3fv(glm::value_ptr(poly.normals[idx]));
			glVertex3fv(glm::value_ptr(poly.vertices[idx]));
		}
		glEnd();
	}

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
			pos = joint.B*pos;
		}
		else{
			pos = joint.W*pos;
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

	joint.B = mat;

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
	joint.W = global_trans;

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

	// rest(bind) poseでのワールド変換行列の計算
	glm::mat4 mat(1.0f);
	calRestGlobalTrans(0, mat, scale);

	// アニメーション時(stepステップ)でのワールド変換行列の計算
	mat = glm::mat4(1.0f);
	calAnimatedGlobalTrans(0, mat, motion+nchannels*(step%m_frames), scale);

	return 1;
}


