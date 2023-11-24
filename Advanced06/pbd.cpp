/*! 
  @file pbd.cpp
	
  @brief PBDによる弾性体シミュレーション
 
  @author Makoto Fujisawa
  @date 2021
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "pbd.h"
#include "msh.h"


//-----------------------------------------------------------------------------
// 課題用関数
//-----------------------------------------------------------------------------
/*!
 * PBDのstretching constraintの適用
 *  - エッジ間の長さを保つようにm_vNewPosを更新
 * @param[in] ks 合成パラメータ(スライド中のk,現在の反復回数で変化するので引数として与える)
 */
void ElasticPBD::projectStretchingConstraint(float ks)
{
	if(m_iNumEdge <= 1) return;

	for(int i = 0; i < m_iNumEdge; ++i){
		// 四面体を使うときの内部エッジかどうかの判定＆内部エッジを使うかどうかのフラグチェック
		if(m_vInEdge[i] && !m_bUseInEdge) continue;

		// エッジ情報の取得とエッジ両端の頂点番号および質量の取得(固定点の質量は大きくする)
		const rxEdge &e = m_poly.edges[i];
		int v1 = e.v[0];	// エッジ
		int v2 = e.v[1];
		float m1 = m_vFix[v1] ? 30.0f*m_vMass[v1] : m_vMass[v1];
		float m2 = m_vFix[v2] ? 30.0f*m_vMass[v2] : m_vMass[v2];
		if(m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>()) continue;

		// 2頂点の位置ベクトル
		glm::vec3 p1 = m_vNewPos[v1];
		glm::vec3 p2 = m_vNewPos[v2];

		// 計算点間の元の長さ(制約条件)
		float d = m_vLengths[i];

		// TODO:重力等を考慮した後の2頂点座標(スライド中のp')はp1=m_vNewPos[v1],p2=m_vNewPos[v2]で得られるので，
		//      これらから制約を満たすような位置修正量dp1,dp2を求めて，m_vNewPos[v1],m_vNewPos[v2]に足し合わせる．
		//      ◎エッジの長さによってはゼロ割が発生することがある．エラーチェックを忘れずに！
		glm::vec3 dp1, dp2;

		// ----課題ここから----



		// ----課題ここまで----

		// 頂点位置を修正
		if(!m_vFix[v1]) m_vNewPos[v1] += ks*dp1;
		if(!m_vFix[v2]) m_vNewPos[v2] += ks*dp2;
	}
}

/*!
 * PBDのbending constraintの適用
 *  - 1つのエッジを共有する2つの三角形間の角度を保つようにm_vNewPosを更新
 * @param[in] ks 合成パラメータ(スライド中のk,現在の反復回数で変化するので引数として与える)
 */
void ElasticPBD::projectBendingConstraint(float ks)
{
	if(m_iNumTris <= 1 || m_iNumEdge <= 0 || m_vBends.empty()) return;

	for(int i = 0; i < m_iNumEdge; i++){
		// 2つのポリゴンに挟まれたエッジ情報の取得
		const rxEdge &e = m_poly.edges[i];
		if(e.f.size() < 2) continue;	// このエッジを含むポリゴン数が1なら処理をスキップ

		// 2つの三角形を構成する4頂点のインデックスを抽出
		set<int>::iterator itr = e.f.begin();
		const rxFace &f1 = m_poly.faces[*itr]; itr++;
		const rxFace &f2 = m_poly.faces[*itr];
		int v1 = e.v[0], v2 = e.v[1], v3, v4;
		for(int j = 0; j < 3; ++j){
			if(f2[j] != v1 && f2[j] != v2) v4 = f2[j];
			if(f1[j] != v1 && f1[j] != v2) v3 = f1[j];
		}
		float m1 = m_vFix[v1] ? 30.0f*m_vMass[v1] : m_vMass[v1];
		float m2 = m_vFix[v2] ? 30.0f*m_vMass[v2] : m_vMass[v2];
		float m3 = m_vFix[v3] ? 30.0f*m_vMass[v3] : m_vMass[v3];
		float m4 = m_vFix[v4] ? 30.0f*m_vMass[v4] : m_vMass[v4];
		if(m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>() || m3 < glm::epsilon<float>() || m4 < glm::epsilon<float>()) continue;

		// 4頂点の位置ベクトル(p2-p4はp1に対する相対位置ベクトル) -> スライドp36の^(ハット)付きのp2-p4の方
		glm::vec3 p1 = m_vNewPos[v1];
		glm::vec3 p2 = m_vNewPos[v2]-p1;
		glm::vec3 p3 = m_vNewPos[v3]-p1;
		glm::vec3 p4 = m_vNewPos[v4]-p1;

		// 2面間の初期角度
		float phi0 = m_vBends[i];

		// TODO:エッジを挟んだ4頂点座標p1,p2,p3,p4からbending constraintを満たす位置修正量dp1～dp4を求め，
		//      m_vNewPos[v1]～m_vNewPos[v4]に足し合わせる．
		//		・↑で定義しているp1～p4で，p2～p4はp1に対する相対座標(スライドp36の^(ハット)付きのp2-p4の方)にしてあるので注意
		//		・三角形ポリゴン間の角度の初期値φ0はm_vBends[i]で得られる(↑で変数phi0に代入済み)
		//      ◎ベクトルの大きさで割るという式が多いが，メッシュの変形によってはゼロ割が発生することがある．エラーチェックを忘れずに！
		//		◎スライドp36のdを計算するときに，-1～1の範囲にあるかをちゃんとチェックして，範囲外ならクランプするように！
		//      - 授業スライドに合わせるためにp1～p4など配列を使わずに書いている．
		//		  配列を使って書き換えても構わないが添え字の違い(配列は0から始まる)に注意．
		glm::vec3 dp1(0.0f), dp2(0.0f), dp3(0.0f), dp4(0.0f);

		// ----課題ここから----



		// ----課題ここまで----

		// 頂点位置を移動
		if(!m_vFix[v1]) m_vNewPos[v1] += ks*dp1;
		if(!m_vFix[v2]) m_vNewPos[v2] += ks*dp2;
		if(!m_vFix[v3]) m_vNewPos[v3] += ks*dp3;
		if(!m_vFix[v4]) m_vNewPos[v4] += ks*dp4;
	}
}

/*!
 * PBDのvolume constraintの適用
 *  - 四面体の体積を保つようにm_vNewPosを更新
 * @param[in] ks 合成パラメータ(スライド中のk,現在の反復回数で変化するので引数として与える)
 */
void ElasticPBD::projectVolumeConstraint(float ks)
{
	if(m_iNumTets <= 1) return;

	for(int i = 0; i < m_iNumTets; i++){
		// 四面体情報(四面体を構成する4頂点インデックス)の取得
		int v1 = m_vTets[i][0], v2 = m_vTets[i][1], v3 = m_vTets[i][2], v4 = m_vTets[i][3];

		// 四面体の4頂点座標と質量の取り出し
		glm::vec3 p1 = m_vNewPos[v1];
		glm::vec3 p2 = m_vNewPos[v2];
		glm::vec3 p3 = m_vNewPos[v3];
		glm::vec3 p4 = m_vNewPos[v4];
		float m1 = m_vFix[v1] ? 30.0f*m_vMass[v1] : m_vMass[v1];
		float m2 = m_vFix[v2] ? 30.0f*m_vMass[v2] : m_vMass[v2];
		float m3 = m_vFix[v3] ? 30.0f*m_vMass[v3] : m_vMass[v3];
		float m4 = m_vFix[v4] ? 30.0f*m_vMass[v4] : m_vMass[v4];
		if(m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>() || m3 < glm::epsilon<float>() || m4 < glm::epsilon<float>()) continue;

		// 四面体の元の体積
		float V0 = m_vVolumes[i];

		// TODO:四面体の4頂点座標p1,p2,p3,p4からvolume constraintを満たす位置修正量dp1～dp4を求め，
		//      m_vNewPos[v1]～m_vNewPos[v4]に足し合わせる．
		//      ◎ベクトルの大きさで割るという式が多いが，メッシュの変形によってはゼロ割が発生することがある．エラーチェックを忘れずに！
		//		- 四面体の体積はスライドに書いてある式を書くのでもよいし，
		//		  calVolume()という四面体の体積計算用関数も用意してあるのでこれを使っても良い
		//      - 授業スライドに合わせるためにp1～p4など配列を使わずに書いている．
		//		  配列を使って書き換えても構わないが添え字の違い(配列は0から始まる)に注意．
		glm::vec3 dp1(0.0f), dp2(0.0f), dp3(0.0f), dp4(0.0f);

		// ----課題ここから----



		// ----課題ここまで----

		// 頂点位置を移動
		if(!m_vFix[v1]) m_vNewPos[v1] += ks*dp1;
		if(!m_vFix[v2]) m_vNewPos[v2] += ks*dp2;
		if(!m_vFix[v3]) m_vNewPos[v3] += ks*dp3;
		if(!m_vFix[v4]) m_vNewPos[v4] += ks*dp4;
	}
}


/*!
 * シミュレーションステップを進める
 */
void ElasticPBD::Update(float dt)
{
	// 外力で速度を更新＆予測位置p'の計算
	calExternalForces(dt);
	// 衝突処理(p'を更新)
	genCollConstraints(dt);

	// 制約条件を満たすための反復処理
	for(int i = 0; i < m_iNmax; ++i){
		// 剛性係数(stiffness)の反復回数による修正
		float ks = 1.0f-pow(1.0f-m_fK, 1.0f/static_cast<float>((i+1)));

		// stretching,bending,volume constraintの適用
		projectStretchingConstraint(ks);
		projectBendingConstraint(ks);
		projectVolumeConstraint(ks);
	}

	// 速度と位置の更新
	integrate(dt);
}


//-----------------------------------------------------------------------------
// PBDクラスの実装
//-----------------------------------------------------------------------------
/*!
 * デフォルトコンストラクタ
 */
ElasticPBD::ElasticPBD(int obj)
{
	m_iObjectNo = obj;
	m_bUseInEdge = true;

	//m_v3Gravity = glm::vec3(0.0, -9.81, 0.0);
	m_v3Gravity = glm::vec3(0.0, -0.5, 0.0);
	m_fWind = 0.0f;

	m_v3Min = glm::vec3(-1.0);
	m_v3Max = glm::vec3(1.0);

	m_fK = 0.9;
	m_iNmax = 3;

	m_fpCollision = 0;

	Clear();
}

/*!
 * デストラクタ
 */
ElasticPBD::~ElasticPBD()
{
}

/*!
 * 全頂点を消去
 */
void ElasticPBD::Clear()
{
	m_iNumVertices = 0;
	m_iNumEdge = 0;
	m_iNumTris = 0;
	m_iNumTets = 0;
	m_vCurPos.clear();
	m_vNewPos.clear();
	m_vMass.clear();
	m_vVel.clear();
	m_vFix.clear();
	m_poly.Clear();
	m_vTets.clear();
	m_vLengths.clear();
	m_vBends.clear();
	m_vVolumes.clear();
	m_vInEdge.clear();
}

/*!
 * 描画関数
 * @param[in] drw 描画フラグ
 */
void ElasticPBD::Draw(int drw)
{
	if(drw & 0x02){
		// エッジ描画における"stitching"をなくすためのオフセットの設定
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}

	// 頂点描画
	if(drw & 0x0001){
		glDisable(GL_LIGHTING);
		glPointSize(5.0);
		glColor3d(1.0, 0.3, 0.3);
		for(int i = 0; i < m_iNumVertices; ++i){
			glm::vec3 v = m_vCurPos[i];
			glPushName(i);
			glBegin(GL_POINTS);
			glVertex3d(v[0], v[1], v[2]);
			glEnd();
			glPopName();
		}
	}

	// エッジ描画
	if(drw & 0x0002){
		glEnable(GL_LIGHTING);
		glColor3d(0.5, 0.9, 0.9);
		glLineWidth(1.0);
		glBegin(GL_LINES);
		for(int i = 0; i < m_iNumEdge; ++i){
			const rxEdge &e = m_poly.edges[i];
			glm::vec3 v0 = m_vCurPos[e.v[0]];
			glm::vec3 v1 = m_vCurPos[e.v[1]];
			glVertex3d(v0[0], v0[1], v0[2]);
			glVertex3d(v1[0], v1[1], v1[2]);
		}
		glEnd();
	}

	// 面描画
	if(drw & 0x0004){
		glEnable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);
		glColor3d(0.1, 0.5, 1.0);
		for(int i = 0; i < m_iNumTris; ++i){
			const rxFace &face = m_poly.faces[i];
			int n = (int)face.vert_idx.size();

			// ポリゴン描画
			glBegin(GL_POLYGON);
			for(int j = 0; j < n; ++j){
				int idx = face.vert_idx[j];
				glVertex3fv(glm::value_ptr(m_vCurPos[idx]));
			}
			glEnd();
		}
	}

	// 四面体描画(ワイヤフレーム)
	if(drw & 0x0008){
		glDisable(GL_LIGHTING);
		glColor3d(0.4, 0.7, 0.9);
		glLineWidth(3.0);
		int edges[6][2] = { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}, {3, 1} };
		for(int i = 0; i < m_iNumTets; ++i){
			vector<int> &tet = m_vTets[i];

			if(drw & 0x0010){
				glm::vec3 mc(0.0);
				// 重心を求める
				for(int j = 0; j < 4; ++j){
					mc += m_vCurPos[tet[j]];
				}
				mc /= 4.0f;
				if(mc[0] < 0.0) continue;
			}

			// 四面体をワイヤフレームで描画
			glBegin(GL_LINES);
			for(int j = 0; j < 6; ++j){
				glVertex3fv(glm::value_ptr(m_vCurPos[tet[edges[j][0]]]));
				glVertex3fv(glm::value_ptr(m_vCurPos[tet[edges[j][1]]]));
			}
			glEnd();
		}
	}

	// 四面体描画(断面)
	if(drw & 0x0010){
		glEnable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);
		glColor3d(0.4, 0.7, 0.4);
		int poly[4][3] = { {0, 2, 1}, {0, 3, 2}, {0, 1, 3}, {1, 2, 3} };
		for(int i = 0; i < m_iNumTets; ++i){
			vector<int> &tet = m_vTets[i];
			vector<glm::vec3> v(4);
			glm::vec3 mc(0.0), n(0.0);

			// 頂点座標を抽出しつつ，重心を求める
			for(int j = 0; j < 4; ++j){
				v[j] = m_vCurPos[tet[j]];
				mc += v[j];
			}
			mc /= 4.0f;

			if(mc[0] >= 0.0){
				// 四面体をポリゴンで描画
				for(int j = 0; j < 4; ++j){
					// 法線
					n = glm::normalize(glm::cross(v[poly[j][1]]-v[poly[j][0]], v[poly[j][2]]-v[poly[j][0]]));

					// 法線と頂点から重心へのベクトルの内積で裏表を判定
					float nv = glm::dot(n, mc-v[poly[j][0]]);
					if(nv > 0.0) n *= -1.0f;

					glBegin(GL_POLYGON);
					glNormal3fv(glm::value_ptr(n));
					for(int k = 0; k < 3; ++k){
						glVertex3fv(glm::value_ptr(v[poly[j][k]]));
					}
					glEnd();
				}
			}

		}
	}

}


/*!
 * 頂点を追加
 * @param[in] pos 頂点位置
 * @param[out] mass 頂点質量
 */
void ElasticPBD::AddVertex(const glm::vec3 &pos, float mass)
{
	m_poly.vertices.push_back(pos);
	m_vCurPos.push_back(pos);
	m_vNewPos.push_back(pos);
	m_vMass.push_back(mass);
	m_vVel.push_back(glm::vec3(0.0));
	m_vFix.push_back(false);
	m_iNumVertices++;
}

/*!
 * エッジを追加
 *  - ポリゴンとの関係などを計算する場合はMakeEdge関数を使うこと
 *  - ポリゴンがない1次元弾性用
 * @param[in] v0,v1 エッジを構成する頂点番号
 */
void ElasticPBD::AddEdge(int v0, int v1)
{
	rxEdge e;
	e.v[0] = v0; e.v[1] = v1;
	m_poly.edges.push_back(e);
	m_vLengths.push_back(glm::length(m_poly.vertices[e.v[0]] - m_poly.vertices[e.v[1]]));
	m_vInEdge.push_back(0);
	m_iNumEdge++;
}

//! 三角体の面積計算(このバージョンでは使わないが面積を保存するような拘束条件もあるので念のため残しておく)
static inline float calArea(const glm::vec3 &p0, const glm::vec3 &p1, const glm::vec3 &p2)
{
	return (1.0f/2.0f)*glm::length(glm::cross(p1-p0, p2-p0));
}

/*!
 * 三角形ポリゴンの追加
 * @param[in] v0,v1,v2 三角形を構成する4頂点のインデックス
 */
void ElasticPBD::AddTriangle(int v0, int v1, int v2)
{
	rxFace f;
	f.vert_idx.resize(3, -1);
	f.vert_idx[0] = v0; f.vert_idx[1] = v1; f.vert_idx[2] = v2;

	m_poly.faces.push_back(f);
	m_iNumTris++;
}

/*!
 * 四面体の追加
 *  - 初期体積の計算を同時に行う必要があるので追加する場合はなるべくこの関数を使う
 * @param[in] v0,v1,v2,v3 四面体を構成する4頂点のインデックス
 */
void ElasticPBD::AddTetra(int v0, int v1, int v2, int v3)
{
	vector<int> tet(4, -1);
	tet[0] = v0; tet[1] = v1; tet[2] = v2; tet[3] = v3;
	m_vTets.push_back(tet);
	m_vVolumes.push_back(calVolume(m_poly.vertices[v0], m_poly.vertices[v1], m_poly.vertices[v2], m_poly.vertices[v3]));
	m_iNumTets++;
}

/*!
 * エッジ重複チェック用のキー(hash値)の計算
 * @param[in] v[2] エッジの2頂点のインデックスを格納した配列
 * @param[in] n hash値計算に使う定数(これが頂点数以上になるようにすること)
 * @return キー(hash値)
 */
static inline long calKey(vector<int> v, int n = 10000)
{
	if(v[0] > v[1]) RX_SWAP(v[0], v[1]);
	return v[0]+v[1]*n;
}

/*!
 * エッジ情報作成
 */
void ElasticPBD::MakeEdge(void)
{
	SearchVertexFace(m_poly);
	SearchEdge(m_poly);
	m_iNumEdge = static_cast<int>(m_poly.edges.size());
	m_vLengths.resize(m_iNumEdge, 0.0);
	m_vInEdge.resize(m_iNumEdge, 0);
	for(int i = 0; i < m_iNumEdge; ++i){
		const rxEdge &e = m_poly.edges[i];
		m_vLengths[i] = glm::length(m_poly.vertices[e.v[0]] - m_poly.vertices[e.v[1]]);
	}

	// 四面体が定義されていたらそのエッジも加える
	if(m_iNumTets){
		map<long, vector<int> > edges;	// エッジリスト(重複チェック用)
		// 重複チェック用に表面エッジをリストに入れておく
		for(int i = 0; i < m_iNumEdge; ++i){
			const rxEdge &e = m_poly.edges[i];
			vector<int> v(2);
			v[0] = e.v[0]; v[1] = e.v[1];
			edges.insert(pair<long, vector<int> >(calKey(v), v));
		}
		// 四面体のエッジを格納
		int edge_idx[6][2] = { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}, {3, 1} };
		for(int i = 0; i < m_iNumTets; ++i){
			vector<int> &tet = m_vTets[i];
			for(int j = 0; j < 6; ++j){
				vector<int> v(2);
				v[0] = tet[edge_idx[j][0]]; v[1] = tet[edge_idx[j][1]];
				long key = calKey(v);
				if(edges.find(key) == edges.end()){
					rxEdge e;
					e.v[0] = v[0]; e.v[1] = v[1];
					m_poly.edges.push_back(e);
					m_vLengths.push_back(glm::length(m_poly.vertices[e.v[0]] - m_poly.vertices[e.v[1]]));
					m_vInEdge.push_back(1);
					m_iNumEdge++;
					edges.insert(pair<long, vector<int> >(key, v));
				}
			}
		}
	}

	CalBends();
}

/*!
 * 三角形ポリゴン間の角度情報作成
 */
void ElasticPBD::CalBends(void)
{
	if(m_iNumTris <= 1 || m_iNumEdge <= 0) return;
	m_vBends.resize(m_iNumEdge, glm::pi<float>());
	for(int i = 0; i < m_iNumEdge; i++){
		const rxEdge &e = m_poly.edges[i];
		if(e.f.size() < 2) continue;

		set<int>::iterator itr = e.f.begin();
		const rxFace &f1 = m_poly.faces[*itr++];
		const rxFace &f2 = m_poly.faces[*itr++];

		// 2つの三角形を構成する4頂点を抽出
		int v[4];
		v[0] = e.v[0];
		v[1] = e.v[1];
		for(int j = 0; j < 3; ++j){
			if(f1[j] != v[0] && f1[j] != v[1]) v[2] = f1[j];
			if(f2[j] != v[0] && f2[j] != v[1]) v[3] = f2[j];
		}

		// 4頂点の位置ベクトル(p2-p4はp1に対する相対位置ベクトル)
		glm::vec3 p1 = m_poly.vertices[v[0]];
		glm::vec3 p2 = m_poly.vertices[v[1]]-p1;
		glm::vec3 p3 = m_poly.vertices[v[2]]-p1;
		glm::vec3 p4 = m_poly.vertices[v[3]]-p1;


		// 2面の法線ベクトル
		glm::vec3 n1 = glm::normalize(glm::cross(p2, p3));
		glm::vec3 n2 = glm::normalize(glm::cross(p2, p4));
		float d = glm::dot(n1, n2);

		// 2面間の角度
		m_vBends[i] = glm::pi<float>();// acos(d);
	}
}


/*!
 * nの頂点を持つ1次元弾性体生成
 * @param[in] c1,c2 2端点座標
 * @param[in] n エッジ分割数(頂点数ではないので注意.頂点数はn+1になる)
 */
void ElasticPBD::GenerateStrand(glm::vec3 c1, glm::vec3 c2, int n)
{
	Clear();

	glm::vec3 dir = c2-c1;
	float l = glm::length(dir);
	float dx = l/(n-1.0f);

	// 頂点座標生成
	double mass = 0.1;				// 各頂点の質量
	for(int i = 0; i < n+1; ++i){
		AddVertex(c1+dir*(static_cast<float>(i)/static_cast<float>(n)), mass);
	}
	for(int i = 0; i < n; ++i){
		AddEdge(i, i+1);
	}
}

/*!
 * n×nの頂点を持つメッシュ生成(x-y平面)
 * @param[in] c1,c2 2端点座標
 * @param[in] nx,ny 生成するメッシュ解像度
 */
void ElasticPBD::GenerateMesh(glm::vec2 c1, glm::vec2 c2, int nx, int ny)
{
	Clear();

	// 頂点座標生成
	float dx = (c2[0]-c1[0])/(nx-1.0);
	float dz = (c2[1]-c1[1])/(ny-1.0);
	double mass = 0.1;				// 各頂点の質量
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i < nx; ++i){
			glm::vec3 pos;
			pos[0] = c1[0]+i*dx;
			pos[1] = c1[1]+j*dz;
			pos[2] = 0.0;
			AddVertex(pos, mass);
		}
	}

	// メッシュ作成
	for(int j = 0; j < ny-1; ++j){
		for(int i = 0; i < nx-1; ++i){
			AddTriangle(i+j*nx, (i+1)+j*nx, (i+1)+(j+1)*nx);
			AddTriangle(i+j*nx, (i+1)+(j+1)*nx, i+(j+1)*nx);
		}
	}

	// エッジ生成
	MakeEdge();
}

/*!
 * MSHファイルから四面体構造を読み込み
 * @param[in] filename MSHファイル名(パス)
 */
void ElasticPBD::GenerateTetrahedralMeshFromFile(const string &filename)
{
	rxMSH msh;
	vector<glm::vec3> points;
	vector< vector<int> > tris, tets;
	if(!msh.Read(filename, points, tris, tets)){
		return;
	}
	if(points.empty() || tets.empty()) return;

	// 頂点生成
	double mass = 0.1;				// 各頂点の質量
	for(int i = 0; i < points.size(); ++i){
		AddVertex(points[i], mass);
	}

	// メッシュ作成
	for(int i = 0; i < tris.size(); ++i){
		AddTriangle(tris[i][0], tris[i][1], tris[i][2]);
	}

	// 四面体生成
	for(int i = 0; i < tets.size(); ++i){
		AddTetra(tets[i][0], tets[i][1], tets[i][2], tets[i][3]);
	}

	// エッジ生成
	MakeEdge();

}


/*!
 * 固定頂点を設定(位置変更も含む)
 * @param[in] i 頂点インデックス
 * @param[in] pos 固定位置
 */
void ElasticPBD::FixVertex(int i, const glm::vec3 &pos)
{
	m_vCurPos[i] = pos;
	m_vNewPos[i] = pos;
	m_vFix[i] = true;
}

/*!
 * 固定頂点を設定
 * @param[in] i 頂点インデックス
 */
void ElasticPBD::FixVertex(int i)
{
	m_vFix[i] = true;
}


/*!
 * 頂点の固定を解除
 * @param[in] i 頂点インデックス
 */
void ElasticPBD::UnFixVertex(int i)
{
	m_vFix[i] = false;
}

/*!
 * 全頂点の固定を解除
 */
void ElasticPBD::UnFixAllVertex(void)
{
	for(int i = 0; i < m_iNumVertices; ++i)	m_vFix[i] = false;
}

/*!
 * 頂点選択(レイと頂点(球)の交差判定)
 * @param[in]  ray_origin,ray_dir レイ(光線)の原点と方向ベクトル
 * @param[out] t 交差があったときの原点から交差点までの距離(媒介変数の値)
 * @param[in]  rad 球の半径(これを大きくすると頂点からマウスクリック位置が多少離れていても選択されるようになる)
 * @return 交差していればその頂点番号，交差していなければ-1を返す
 */
int ElasticPBD::IntersectRay(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, float &t, float rad)
{
	int v = -1;
	float min_t = 1.0e6;
	float rad2 = rad*rad;
	float a = glm::length2(ray_dir);
	if(a < 1.0e-6) return -1;

	for(int i = 0; i < m_iNumVertices; ++i){
		glm::vec3 cen = m_vCurPos[i];
		glm::vec3 s = ray_origin-cen;
		float b = 2.0f*glm::dot(s, ray_dir);
		float c = glm::length2(s)-rad2;

		float D = b*b-4.0f*a*c;
		if(D < 0.0f) continue;	// 交差なし

		float t0 = (-b-sqrt(D))/(2.0*a);
		float t1 = (-b+sqrt(D))/(2.0*a);
		if(t0 > 0.0 && t1 > 0.0 && t0 < min_t){	// 2交点がある場合
			v = i; min_t = t0;
		}
		else if(t0 < 0.0 && t1 > 0.0 && t1 < min_t){	// 1交点のみの場合(光線の始点が球内部にある)
			v = i; min_t = t1;
		}
	}
	if(v != -1) t = glm::length(ray_origin-m_vCurPos[v]);
	return v;
}



/*!
 * 外力
 *  - 重力と境界壁からの力の影響
 */
void ElasticPBD::calExternalForces(float dt)
{
	// 重力の影響を付加
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_vFix[i]) continue;
		m_vVel[i] += m_v3Gravity*dt+glm::vec3(m_fWind, 0.0f, 0.0f);
		m_vNewPos[i] = m_vCurPos[i]+m_vVel[i]*dt;
	}
}

/*!
* 衝突判定
* @param[in] dt タイムステップ幅
*/
void ElasticPBD::genCollConstraints(float dt)
{
	// 境界壁の影響
	float res = 0.9;	// 反発係数
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_vFix[i]) continue;
		glm::vec3 &p = m_vCurPos[i];
		glm::vec3 &np = m_vNewPos[i];
		glm::vec3 &v = m_vVel[i];
		if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
			np[0] = p[0]-v[0]*dt*res;
			np[1] = p[1];
			np[2] = p[2];
		}
		if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
			np[1] = p[1]-v[1]*dt*res;
			np[0] = p[0];
			np[2] = p[2];
		}
		if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
			np[2] = p[2]-v[2]*dt*res;
			np[0] = p[0];
			np[1] = p[1];
		}
		clamp(m_vNewPos[i]);
	}

	// 他のオブジェクトとの衝突
	if(m_fpCollision != 0){
		for(int i = 0; i < m_iNumVertices; ++i){
			if(m_vFix[i]) continue;
			glm::vec3 &p = m_vCurPos[i];
			glm::vec3 &np = m_vNewPos[i];
			glm::vec3 &v = m_vVel[i];
			m_fpCollision(p, np, v, m_iObjectNo);
		}
	}
}


/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
void ElasticPBD::integrate(float dt)
{
	float dt1 = 1.0f/dt;
	for(int i = 0; i < m_iNumVertices; ++i){
		m_vVel[i] = (m_vNewPos[i]-m_vCurPos[i])*dt1;
		m_vCurPos[i] = m_vNewPos[i];
	}
}
