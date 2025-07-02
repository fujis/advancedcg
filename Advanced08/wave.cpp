/*! 
  @file wave.cpp
	
  @brief SWEを使った波のシミュレーション
 
  @author Makoto Fujisawa
  @date 2021-05
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "wave.h"


//-----------------------------------------------------------------------------
// 課題用関数
//-----------------------------------------------------------------------------
/*!
 * SWEによるハイトフィールドの更新
 *  - 移流項のセミラグランジュ法による計算
 *  - *_newの方が更新後の値で出力値となる
 * @param[out] d_new 更新後の水深(デプス)値を格納する配列
 * @param[in] d 更新前の水深(デプス)値を格納した配列
 * @param[out] u_new,v_new 更新後の速度場を格納する配列
 * @param[in] u,v 更新前の速度場を格納した配列
 * @param[in] dt 時間ステップ幅
 */
void WaveSWE::advection(float *d_new, float *d, float *u_new, float *v_new, float *u, float *v, float dt)
{
	float dx = m_dx;
	float dy = m_dy;

	// TODO:この部分にSWE計算での移流項の計算コードを書く(セミラグランジュ法で計算すること)
	// ・配列d_new,u_new,v_newにそれぞれ移流項を適用した後の水深,x方向速度,y方向速度を格納する
	// ・移流項適用前の水深,速度は配列d,u,vに格納されている
	//   ⇒ 1次元配列を使っているので，取り出すときはd[IDX(i,j)]のようにIDX関数を使うと便利
	// ・各グリッドセルの幅は変数dx,dy(もしくはm_dx,m_dy)に入っている
	// ・グリッドセル数はm_nx,m_nyで，境界を除いた部分を処理した後，境界処理関数bndを呼び出すのが基本的な流れ
	//  for(int j = 1; j < m_ny-1; ++j){
	//  	for(int i = 1; i < m_nx-1; ++i){
	//  		// d,u,vの更新処理
	//			// ⇒d_new,u_new,v_newに結果を格納
	//  	}
	//  }
	//  bnd(d_new);
	//  bnd(u_new, v_new);
	//  ・移流計算でのd,u,vは同じ反復ループ内で計算してもOK．
	//    ただし，バックトレース位置を求める速度には更新したu_new,v_newではなく移流前の値u,vを使うこと
	//  ・バックトレースした位置がグリッド範囲外に成らないようにチェックが必要(glm::clamp<int>(x, min, max)を上手く使おう)
	//  ・高さや速度が定義されているのはグリッドセル中心ですが，最小インデックスをもつグリッドセルの中心を原点とした座標系を仮定しているので，
	//    グリッドセル中心座標は(i*dx, j*dx)でOK((i+0.5)*dxとかにしなくてもよい)

	// ----課題ここから----



	// ----課題ここまで----

}


/*!
 * SWEによるハイトフィールドの更新
 *  - 圧力項の計算
 *  - *_newの方が更新後の値で出力値となる
 * @param[out] d_new 更新後の水深(デプス)値を格納する配列
 * @param[in] d 更新前の水深(デプス)値を格納した配列
 * @param[out] u_new,v_new 更新後の速度場を格納する配列
 * @param[in] u,v 更新前の速度場を格納した配列
 * @param[in] dt 時間ステップ幅
 */
void WaveSWE::pressure(float *d_new, float *d, float *u_new, float *v_new, float *u, float *v, float dt)
{
	float dx = m_dx;
	float dy = m_dy;
	float g = m_gravity;

	// TODO:この部分にSWE計算での圧力項の計算コードを書く(中心差分で計算すること)
	// ・配列d_new,u_new,v_newにそれぞれ移流項を適用した後の水深,x方向速度,y方向速度を格納する
	// ・移流項適用前の水深,速度は配列d,u,vに格納されている
	//   ⇒ 1次元配列を使っているので，取り出すときはd[IDX(i,j)]のようにIDX関数を使うと便利
	// ・水面高さhはメンバ変数m_hに格納されていて，Update関数内で呼ばれているupdateHeight関数でその値が更新されている．
	// ・各グリッドセルの幅は変数dx,dy(もしくはm_dx,m_dy)に入っている
	// ・グリッドセル数はm_nx,m_nyで，境界を除いた部分を処理した後，境界処理関数bndを呼び出すのが基本的な流れ
	//  for(int j = 1; j < m_ny-1; ++j){
	//  	for(int i = 1; i < m_nx-1; ++i){
	//  		// u,vの更新処理
	//			// ⇒u_new,v_newに結果を格納
	//  	}
	//  }
	//  bnd(u_new, v_new);
	//  for(int j = 1; j < m_ny-1; ++j){
	//  	for(int i = 1; i < m_nx-1; ++i){
	//  		// dの更新処理(u_new,v_newを使う)
	//			// ⇒d_newに結果を格納
	//  	}
	//  }
	//  bnd(d_new);
	//  ・こちらは移流計算と違ってu,vを更新するループを回し終わった後に，別のループでdを更新することになる．
	//    dを更新するときは更新済みのu_new,v_newを使うこと

	// ----課題ここから----



	// ----課題ここまで----

}



/*!
 * ハイトフィールドの更新
 * @param[in] step 現在の計算ステップ数
 * @param[in] dt 時間ステップ幅
 */
int WaveSWE::Update(int step, float dt)
{
	// 強制的な波の生成
	if(m_surf) makeSurf(step*dt, &m_dprev[0], 0.1);

	// SWEによるハイトフィールドの更新

	// 移流項(*_prev ⇒ *)
	advection(&m_d[0], &m_dprev[0], &m_u[0], &m_v[0], &m_uprev[0], &m_vprev[0], dt);

	// 粘性項(速度u,vは* ⇒ *_prev)
	viscosity(&m_uprev[0], &m_vprev[0], &m_u[0], &m_v[0], dt);

	// 水面高さ場hの更新(h=b+d)
	updateHeight(&m_d[0]);

	// 圧力項(水深dは* ⇒ *_prev，速度u,vは*_prev ⇒ *)
	pressure(&m_dprev[0], &m_d[0], &m_u[0], &m_v[0], &m_uprev[0], &m_vprev[0], dt);

	// 水面高さ場hの再更新と描画用メッシュの更新
	updateHeight(&m_dprev[0]);
	updateMesh(m_h);

	// 次のステップのためにu_prev,v_prevを更新しておく
	for(int i = 0; i < m_nx*m_ny; ++i){ m_uprev[i] = m_u[i]; m_vprev[i] = m_v[i]; }

	return 1;
}





//-----------------------------------------------------------------------------
// WaveSWEクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 */
WaveSWE::WaveSWE()
{
	// 重力
	m_gravity = 9.81;
	m_surf = false;
	m_gs = false;
	m_nu = 1.0e-3;
}

/*!
 * デストラクタ
 */
WaveSWE::~WaveSWE()
{
}

/*!
 * 波の初期化
 *  - 水面がy=0になるように設定する
 * @param[in] n グリッド数
 * @param[in] scale 全体のスケール
 * @param[in] ground 水底の高さを与える関数ポインタ
 */
void WaveSWE::Init(int n, float scale, float(*ground)(float, float))
{
	// ハイトフィールドの解像度と全体の大きさ
	m_nx = n;
	m_ny = n;
	m_scale = scale;

	// 水底の高さを与える関数
	m_ground = ground;

	// 配列の初期化
	m_h.resize(m_nx*m_ny, 0.0);
	m_d.resize(m_nx*m_ny, 0.0);
	m_dprev.resize(m_nx*m_ny, 0.0);
	m_u.resize(m_nx*m_ny, 0.0);
	m_v.resize(m_nx*m_ny, 0.0);
	m_uprev.resize(m_nx*m_ny, 0.0);
	m_vprev.resize(m_nx*m_ny, 0.0);

	// メッシュ作成(グリッド幅m_dx,m_dyの計算もgenerateMeshでやっている)
	glm::vec3 c1(-m_scale/2.0f, 0.0f, -m_scale/2.0f);
	glm::vec3 c2(m_scale/2.0f, 0.0f, m_scale/2.0f);
	generateMesh(c1, c2, m_mesh);
	generateMesh(c1, c2, m_mesh_ground);

	//! 波の初期化
	for(int j = 0; j < m_ny; ++j){
		for(int i = 0; i < m_nx; ++i){
			int idx = IDX(i, j);
			float b = m_ground(i*m_dx, j*m_dy);
			m_d[idx] = m_dprev[idx] = -b;	// 水面の高さ0なので水深は水底の高さに-1を掛けたものになる
			m_u[idx] = m_uprev[idx] = 0.0;
			m_v[idx] = m_vprev[idx] = 0.0;
		}
	}

	// 水面高さの更新
	updateHeight(&m_d[0]);

	// 水面ハイトフィールドメッシュの更新
	updateMesh(m_h);

	// 水底メッシュの更新(水底地形は変わらないと仮定しているので最初だけでよい)
	for(int j = 0; j < m_ny; ++j){
		for(int i = 0; i < m_nx; ++i){
			int idx = IDX(i, j);
			m_mesh_ground.vertices[idx][1] = m_ground(i*m_dx, j*m_dy);
		}
	}

	// 法線再計算
	CalVertexNormals(m_mesh);
	CalVertexNormals(m_mesh_ground);
}



/*!
 * OpenGLによるハイトフィールドメッシュの描画
 * @param[in] drw 描画フラグ
 */
void WaveSWE::Draw(int draw)
{
	if(draw & 0x02){
		// エッジ描画における"stitching"をなくすためのオフセットの設定
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}

	if(draw & 0x04){
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_LIGHTING);

		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnable(GL_COLOR_MATERIAL);

		// 水面をポリゴンで描画
		glColor3d(0.2, 0.2, 0.7);
		for(const rxFace &face : m_mesh.faces){
			glBegin(GL_POLYGON);
			for(int idx : face.vert_idx){
				glNormal3fv(glm::value_ptr(m_mesh.normals[idx]));
				glVertex3fv(glm::value_ptr(m_mesh.vertices[idx]));
			}
			glEnd();
		}
		// 水底をポリゴンで描画
		glColor3d(0.7, 0.4, 0.2);
		for(const rxFace &face : m_mesh_ground.faces){
			glBegin(GL_POLYGON);
			for(int idx : face.vert_idx){
				glNormal3fv(glm::value_ptr(m_mesh_ground.normals[idx]));
				glVertex3fv(glm::value_ptr(m_mesh_ground.vertices[idx]));
			}
			glEnd();
		}
	}


	// 頂点描画
	if(draw & 0x01){
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glDisable(GL_LIGHTING);
		glPointSize(5.0);
		glColor3d(1.0, 0.3, 0.3);
		glBegin(GL_POINTS);
		for(const glm::vec3 &v : m_mesh.vertices){
			glVertex3fv(glm::value_ptr(v));
		}
		glEnd();
	}

	// エッジ描画
	if(draw & 0x02){
		glDisable(GL_LIGHTING);
		glColor3d(0.5, 0.5, 0.5);
		glLineWidth(1.0);
		for(const rxFace &face : m_mesh.faces){
			glBegin(GL_LINE_LOOP);
			for(int idx : face.vert_idx){
				glVertex3fv(glm::value_ptr(m_mesh.vertices[idx]));
			}
			glEnd();
		}
	}
}


/*!
 * n×nの頂点を持つメッシュ生成(x-z平面)
 * @param[in] c1,c2 2端点座標

 */
void WaveSWE::generateMesh(glm::vec3 c1, glm::vec3 c2, rxPolygons &poly)
{
	if(!poly.vertices.empty()){
		poly.vertices.clear();
		poly.faces.clear();
	}

	// 頂点座標生成
	float dx = (c2[0]-c1[0])/static_cast<float>(m_nx-1);
	float dz = (c2[2]-c1[2])/static_cast<float>(m_ny-1);
	poly.vertices.resize(m_nx*m_ny);
	for(int k = 0; k < m_ny; ++k){
		for(int i = 0; i < m_nx; ++i){
			glm::vec3 pos;
			pos[0] = c1[0]+i*dx;
			pos[1] = c1[1];
			pos[2] = c1[2]+k*dz;
			poly.vertices[IDX(i, k)] = pos;
		}
	}

	m_min = c1+glm::vec3(0.0f, -0.5f*dx*m_nx, 0.0f);
	m_max = c1+glm::vec3((m_nx-1)*dx, 0.5f*dx*m_nx, (m_ny-1)*dz);

	m_dx = dx;
	m_dy = dz;

	// メッシュ作成
	for(int k = 0; k < m_ny-1; ++k){
		for(int i = 0; i < m_nx-1; ++i){
			rxFace face;
			face.resize(3);

			face[0] = IDX(i, k);
			face[1] = IDX(i+1, k+1);
			face[2] = IDX(i+1, k);
			poly.faces.push_back(face);

			face[0] = IDX(i, k);
			face[1] = IDX(i, k+1);
			face[2] = IDX(i+1, k+1);
			poly.faces.push_back(face);
		}
	}

	// 頂点法線の更新
	CalVertexNormals(poly.vertices, poly.vertices.size(), poly.faces, poly.faces.size(), poly.normals);

}


/*!
 * ハイトフィールドに従って描画用メッシュ頂点のy座標値を更新
 * @param[in] h ハイトフィールド(m_nx*m_ny)
 */
void WaveSWE::updateMesh(const vector<float> &h)
{
	for(int k = 0; k < m_ny; ++k){
		for(int i = 0; i < m_nx; ++i){
			int idx = IDX(i, k);
			m_mesh.vertices[idx][1] = h[idx];
		}
	}

	// 頂点法線の更新
	CalVertexNormals(m_mesh.vertices, m_mesh.vertices.size(), m_mesh.faces, m_mesh.faces.size(), m_mesh.normals);
}



/*!
 * 頂点選択(レイと頂点(球)の交差判定)
 * @param[in]  ray_origin,ray_dir レイ(光線)の原点と方向ベクトル
 * @param[out] t 交差があったときの原点から交差点までの距離(媒介変数の値)
 * @param[in]  rad 球の半径(これを大きくすると頂点からマウスクリック位置が多少離れていても選択されるようになる)
 * @return 交差していればその頂点番号，交差していなければ-1を返す
 */
int WaveSWE::IntersectRay(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, float &t, float rad)
{
	if(m_mesh.vertices.empty()) return -1;
	int v = -1;
	float min_t = 1.0e6;
	float rad2 = rad*rad;
	float a = glm::length2(ray_dir);
	if(a < 1.0e-6) return -1;

	glm::vec3 origin(m_mesh.vertices[0][0], 0.0, m_mesh.vertices[0][2]);
	for(int j = 1; j < m_ny-1; ++j){
		for(int i = 1; i < m_nx-1; ++i){
			int idx = IDX(i, j);
			float x = i*m_dx, z = j*m_dy; // グリッド位置

			glm::vec3 cen = glm::vec3(x, m_h[IDX(i, j)], z)+origin;
			glm::vec3 s = ray_origin-cen;
			float b = 2.0f*glm::dot(s, ray_dir);
			float c = glm::length2(s)-rad2;

			float D = b*b-4.0f*a*c;
			if(D < 0.0f) continue;	// 交差なし

			float t0 = (-b-sqrt(D))/(2.0*a);
			float t1 = (-b+sqrt(D))/(2.0*a);
			if(t0 > 0.0 && t1 > 0.0 && t0 < min_t){	// 2交点がある場合
				v = idx; min_t = t0;
			}
			else if(t0 < 0.0 && t1 > 0.0 && t1 < min_t){	// 1交点のみの場合(光線の始点が球内部にある)
				v = idx; min_t = t1;
			}
		}
	}
	return v;
}

/*!
 * 高さ値の直接変更
 * @param[in] idx 頂点インデックス
 * @param[in] h 追加する高さ
 */
void WaveSWE::AddHeight(int idx, float dh)
{
	if(idx < 0 || idx >= m_nx*m_ny) return;
	m_d[idx] += dh;
	m_dprev[idx] = m_d[idx];
	float b = m_ground((idx%m_nx)*m_dx, (idx/m_nx)*m_dy);	// 水底の地形の高さ
	m_h[idx] = m_d[idx]+b;
	m_mesh.vertices[idx][1] = m_h[idx];
}
void WaveSWE::AddHeight(int i, int j, float h)
{
	AddHeight(IDX(i, j), h);
}
void WaveSWE::AddHeightArea(int idx, float dh)
{
	if(idx < 0 || idx >= m_nx*m_ny) return;
	int s = 3;
	int i0 = idx%m_nx;
	int j0 = idx/m_nx;
	for(int j = 0; j < m_ny; ++j){
		for(int i = 0; i < m_nx; ++i){
			if(i >= i0-s && i <= i0+s && j >= j0-s && j <= j0+s){
				int idx = IDX(i, j);
				float b = m_ground(i*m_dx, j*m_dy);	// 水底の地形の高さ
				m_d[idx] = m_dprev[idx] = m_avg_h-b + dh;
			}
		}
	}
	updateHeight(&m_d[0]);
	updateMesh(m_h);
}

/*!
 * 平均高さの算出
 * @return 平均高さ
 */
double WaveSWE::calAvarageHeight(void)
{
	int n = m_nx*m_ny;
	if(n <= 0) return 0.0f;

	float hsum = 0.0f;
	for(const float &h : m_h){
		hsum += h;
	}
	return hsum/static_cast<float>(n);
}



/*!
 * ハイトフィールドの更新
 * @param[in] t 現在のシミュレーション時間(step*dt)
 * @param[out] h ハイトフィールド
 * @param[in] wave_height 設定する波の高さ
 */
void WaveSWE::makeSurf(float t, float *h, float wave_height)
{
	float ht = wave_height*sin(4*glm::pi<float>()*t);
	for(int i = 0; i < m_nx; ++i){
		float b = m_ground(i*m_dx, 1*m_dy);	// 水底の地形の高さ
		h[IDX(i, 1)] = -b + ht;
	}
}




/*!
 * 周囲境界条件の処理(水深場)
 * @param[in] d 水深の値が格納された配列
 */
void WaveSWE::bnd(float *d)
{
	for(int i = 0; i < m_nx; ++i){
		d[IDX(i, 0)] = d[IDX(i, 1)]; d[IDX(i, m_ny-1)] = d[IDX(i, m_ny-2)];
	}
	for(int j = 0; j < m_ny; ++j){
		d[IDX(0, j)] = d[IDX(1, j)]; d[IDX(m_nx-1, j)] = d[IDX(m_nx-2, j)];
	}
}
/*!
 * 周囲境界条件の処理(速度場)
 * @param[in] u,v 速度が格納された配列
 */
void WaveSWE::bnd(float *u, float *v)
{
	for(int i = 0; i < m_nx; ++i){
		u[IDX(i, 0)] = u[IDX(i, 1)]; u[IDX(i, m_ny-1)] = u[IDX(i, m_ny-2)];	// 壁に平行な成分は単純にコピー
		//v[IDX(i, 0)] = v[IDX(i, m_ny-1)] = 0.0f;
		v[IDX(i, 0)] = v[IDX(i, 1)]; v[IDX(i, m_ny-1)] = v[IDX(i, m_ny-2)];	// 手前と奥の面は波が抜けるようにする
	}
	for(int j = 0; j < m_ny; ++j){
		u[IDX(0, j)] = u[IDX(m_nx-1, j)] = 0.0f;	// 左右の面は壁面
		v[IDX(0, j)] = v[IDX(1, j)]; v[IDX(m_nx-1, j)] = v[IDX(m_nx-2, j)];	// 壁に平行な成分は単純にコピー
	}
}


/*!
 * 水深と地形の高さから水面の高さを更新
 * @param[in] d 水深の値が格納された配列
 */
void WaveSWE::updateHeight(float *d)
{
	// デプス値dから地形を含めた高さhを計算
	for(int j = 0; j < m_ny; ++j){
		for(int i = 0; i < m_nx; ++i){
			int idx = IDX(i, j);
			float b = m_ground(i*m_dx, j*m_dy);	// 水底の地形の高さ
			m_h[idx] = d[idx]+b;
		}
	}

	// マウス入出力のために平均高さを求めておく
	//  - 境界条件によっては高さが変わってしまう
	m_avg_h = calAvarageHeight();
}

/*!
 * SWEによるハイトフィールドの更新
 *  - 速度の粘性拡散項の計算
 *  - *_newの方が更新後の値で出力値となる
 *  - この項は速度のみに適用する
 * @param[in] u,v 更新前の速度場を格納した配列
 * @param[out] u_new,v_new 更新後の速度場を格納する配列
 * @param[in] dt 時間ステップ幅
 */
void WaveSWE::viscosity(float *u_new, float *v_new, float *u, float *v, float dt)
{
	if(m_nu <= 0.0){ m_nu = 0.0; return; }
	float dx = m_dx;
	float dy = m_dy;

	// u_new,v_newに値をコピーしておく
	for(int i = 0; i < m_nx*m_ny; ++i){ u_new[i] = u[i]; v_new[i] = v[i]; }

	// ガウス・ザイデル法で粘性拡散項を解く
	for(int k = 0; k < 10; ++k){	// ガウス・ザイデル反復
		for(int j = 1; j < m_ny-1; ++j){
			for(int i = 1; i < m_nx-1; ++i){
				int idx = IDX(i, j);

				// 正方グリッドを仮定する(dx=dy)ならもう少し式をまとめられるが，分かりやすくするために中心差分式そのまま計算
				u_new[idx] = u_new[idx] + dt*m_nu*((u_new[IDX(i+1, j)]-2*u_new[idx]+u_new[IDX(i-1, j)])/(dx*dx) + (u_new[IDX(i, j+1)]-2*u_new[idx]+u_new[IDX(i, j-1)])/(dy*dy));
				v_new[idx] = v_new[idx] + dt*m_nu*((v_new[IDX(i+1, j)]-2*v_new[idx]+v_new[IDX(i-1, j)])/(dx*dx) + (v_new[IDX(i, j+1)]-2*v_new[idx]+v_new[IDX(i, j-1)])/(dy*dy));
			}
		}
		bnd(u_new, v_new);
	}
}

