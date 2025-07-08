/*! 
  @file rx_deform.cpp
	
	@brief 2Dメッシュ変形

  @author Makoto Fujisawa
  @date 2021-03
  */


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "deformation.h"

#include "rx_sampler.h"
#include "rx_delaunay.h"


//-----------------------------------------------------------------------------
// 課題用関数
//-----------------------------------------------------------------------------
/*!
* メッシュ変形 by MLS
*  - Affine Deformation
* @param[in] v 変形する頂点座標
* @param[in] pc 制御点重心
* @param[in] qc 変形後の制御点重心
* @param[in] alpha 重みwを計算するための係数(w=1/(|p-v|^(2*alpha)))
* @return 変形後の座標(f(v))
*/
glm::vec2 rxMeshDeform2D::affineDeformation(const glm::vec2 &v, const glm::vec2 &pc, const glm::vec2 &qc, const double alpha)
{
	// TODO:この部分で 入力頂点座標v と 制御点座標p,q から入力頂点の変形後の座標を計算するコードを書く
	// - Ajを計算してから変形後の座標を計算するのでも，スライドp43の式を直接計算するのでもどちらでもOK
	//   (Ajの前計算は今回は行わないでもよい)
	// - glmでの行列×縦ベクトルの計算は行列とベクトルをそれぞれglm::mat2,glm::vec2で定義しているならば単純な掛け算として計算出来る
	//    例) glm::mat2 m(1.0f);
	//        glm::vec2 v(0.0f);
	//        glm::vec2 mv = m*v;
	//   ただし，横ベクトル×行列の計算部分はglmのオペレータ*は使わない方がよい(glmは縦ベクトルを前提としている)．
	//        glm::vec2 vm;
	// 	      vm[0] = v[0]*m[0][0]+v[1]*m[1][0];
	//        vm[1] = v[0]*m[0][1]+v[1]*m[1][1];
	// - 逆行列計算にはglm::inverse関数が使えるが，行列の正則性チェックはしていないようなので，
	//   glm::determinant関数などで行列式の値を求めて0でないかをチェックするエラー処理を忘れないようにしよう．

	//
	// - 制御点でループして相対座標を計算するまでのコード例: 
	//for(int k = 0; k < m_iNcp; ++k){	// 制御点数(m_iNcp)でループを回す
	//	int j = m_vCP[k];	// 制御点の頂点インデックス
	//
	//	// 重心を中心とした制御点の相対座標
	//	// 各頂点の座標は配列m_vPとm_vXに格納されている(それぞれ初期形状と変形後の形状)
	//	glm::vec2 p = m_vP[j]-pc;	// 初期形状での座標
	//	glm::vec2 q = m_vX[j]-qc;	// 変形後の座標
	//
	//	// ここに色々計算するコードを書く
	//
	//}

	// 変形後の頂点vの座標
	glm::vec2 fa = v;	// ここも書き換えること

	// ----課題ここから----



	// ----課題ここまで----

	return fa;
}

/*!
* メッシュ変形 by MLS
*  - Similarity Deformation
* @param[in] v 変形する頂点座標
* @param[in] pc 制御点重心
* @param[in] qc 変形後の制御点重心
* @param[in] alpha 重みwを計算するための係数(w=1/(|p-v|^(2*alpha)))
* @return 変形後の座標(f(v))
*/
glm::vec2 rxMeshDeform2D::similarityDeformation(const glm::vec2 &v, const glm::vec2 &pc, const glm::vec2 &qc, const double alpha)
{
	// TODO:この部分で 入力頂点座標v と 制御点座標p,q から入力頂点の変形後の座標を計算するコードを書く
	// - μsを先に計算してから変形後の頂点位置を計算
	// - 行列A_iの前計算は今回は行わないでもよい
	// - 横ベクトル×行列の計算部分はglmのオペレータ*は使わない方がよい(glmは縦ベクトルを前提としている)

	// 変形後の頂点vの座標
	glm::vec2 fsv = v;	// ここも書き換えること

	// ----課題ここから----



	// ----課題ここまで----

	return fsv;
}

/*!
* メッシュ変形 by MLS
*  - Rigid Deformation
* @param[in] v 変形する頂点座標
* @param[in] pc 制御点重心
* @param[in] qc 変形後の制御点重心
* @param[in] alpha 重みwを計算するための係数(w=1/(|p-v|^(2*alpha)))
* @return 変形後の座標(f(v))
*/
glm::vec2 rxMeshDeform2D::rigidDeformation(const glm::vec2 &v, const glm::vec2 &pc, const glm::vec2 &qc, const double alpha)
{
	// TODO:この部分で 入力頂点座標v と 制御点座標p,q から入力頂点の変形後の座標を計算するコードを書く
	// - μfを先に計算してから変形後の頂点位置を計算
	// - 行列A_iの前計算は今回は行わないでもよい
	// - 横ベクトル×行列の計算部分はglmのオペレータ*は使わない方がよい(glmは縦ベクトルを前提としている)

	// 変形後の頂点vの座標
	glm::vec2 frv = v;	// ここも書き換えること

	// ----課題ここから----



	// ----課題ここまで----

	return frv;
}



/*!
* メッシュ更新
* @param[in] dt 時間ステップ幅(このメッシュ変形法では使わない)
*/
int rxMeshDeform2D::Update(double dt)
{
	if(m_iNcp <= 1) return 0;

	// 各頂点を変形
	for(int i = 0; i < m_iNv; ++i){
		// 制御点はユーザー入力位置で固定なので処理をスキップ
		if(std::find(m_vCP.begin(), m_vCP.end(), i) != m_vCP.end()) continue;

		// 頂点の初期座標
		const glm::vec2 &v = m_vP[i];

		// 固定点の移動前，移動後の重み付き中心p*,q*の計算
		glm::vec2 pc(0.0), qc(0.0);
		double wsum = 0.0;
		for(int k = 0; k < m_iNcp; ++k){
			int j = m_vCP[k];
			const glm::vec2 &p = m_vP[j];
			const glm::vec2 &q = m_vX[j];

			// 固定点と計算点の間の距離に基づく重み
			double dist = glm::length2(p-v);
			float w = (dist > 1.0e-6) ? 1.0f/pow(dist, m_alpha) : 0.0f;

			pc += w*p;
			qc += w*q;
			wsum += w;
		}
		pc /= wsum;
		qc /= wsum;

		// MLS Deformations
		switch(m_deform_type){
		case 0: m_vX[i] = affineDeformation(v, pc, qc, m_alpha); break;
		case 1: m_vX[i] = similarityDeformation(v, pc, qc, m_alpha); break;
		case 2: m_vX[i] = rigidDeformation(v, pc, qc, m_alpha); break;
		}
	}

	return 1;
}



//-----------------------------------------------------------------------------
// rxMeshDeform2Dクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 * @param[in] n グリッド数
 */
rxMeshDeform2D::rxMeshDeform2D()
{
	m_iNv = m_iNt = m_iNcp = 0;

	m_vao_mesh = 0;
	m_vao_fix = 0;

	m_alpha = 1.0;
	m_deform_type = 0;

	Init(0);
}

/*!
 * デストラクタ
 */
rxMeshDeform2D::~rxMeshDeform2D()
{
	if(m_vao_mesh) glDeleteVertexArrays(1, &m_vao_mesh);
	if(m_vao_fix) glDeleteVertexArrays(1, &m_vao_fix);
}

/*!
 * メッシュ初期化
 */
void rxMeshDeform2D::Init(int random_mesh)
{
	// メッシュ作成
	glm::vec2 c1(-1.0, -1.0);
	glm::vec2 c2(1.0, 1.0);
	if(random_mesh) generateRandomMesh(c1, c2, 0.07, 800);
	else generateMesh(c1, c2, 32, 32);
	
	// 頂点配列オブジェクトの作成
	if(m_vao_mesh != 0) glDeleteVertexArrays(1, &m_vao_mesh);
	m_vao_mesh = CreateVAO((GLfloat*)&m_vX[0], m_iNv, 2, &m_vTri[0], m_iNt, 0, 0, 0, 0, (GLfloat*)&m_vTC[0], m_iNv);

	// 固定点の設定
	m_vCP.clear(); m_iNcp = 0;
	updateCPVAO();
}


/*!
 * 近傍頂点探索
 *  - マウスピック用
 *  - 探索半径h内で最も近い点を返す
 * @param[in] pos 探索位置
 * @param[in] h 探索半径
 * @return 近傍頂点インデックス
 */
int rxMeshDeform2D::Search(glm::vec2 pos, double h)
{
	int idx = -1;
	double min_d2 = RX_FEQ_INF;
	double h2 = h*h;
	for(int i = 0; i < m_iNv; ++i){
		double d2 = glm::length2(m_vX[i]-pos);
		if(d2 < h2 && d2 < min_d2){
			min_d2 = d2;
			idx = i;
		}
	}
	return idx;
}
/*!
* 近傍固定頂点探索
*  - マウスピック用
*  - 探索半径h内で最も近い点を返す
* @param[in] pos 探索位置
* @param[in] h 探索半径
* @return 近傍頂点インデックス
*/
int rxMeshDeform2D::SearchCP(glm::vec2 pos, double h)
{
	int idx = -1;
	double min_d2 = RX_FEQ_INF;
	double h2 = h*h;
	for(int k = 0; k < m_iNcp; ++k){
		int i = m_vCP[k];
		double d2 = glm::length2(m_vX[i]-pos);
		if(d2 < h2 && d2 < min_d2){
			min_d2 = d2;
			idx = i;
		}
	}
	return idx;
}

/*!
* 制御点座標VAOの更新
*  - 制御点位置を変更した場合，その描画に用いているVAOの中身も更新する必要がある
*/
void rxMeshDeform2D::updateCPVAO(void)
{
	if(m_vCP.empty()){
		if(m_vao_fix != 0){
			glDeleteVertexArrays(1, &m_vao_fix);
			m_vao_fix = 0;
		}
	}
	else{
		vector<glm::vec2> fixpos;
		for(int i : m_vCP) fixpos.push_back(m_vX[i]);
		if(m_vao_fix != 0) glDeleteVertexArrays(1, &m_vao_fix);
		m_vao_fix = CreateVAO((GLfloat*)&fixpos[0], m_iNcp, 2);
	}
}


//! 固定点設定
void rxMeshDeform2D::SetCP(int idx, glm::vec2 pos, bool move)
{
	if(std::find(m_vCP.begin(), m_vCP.end(), idx) == m_vCP.end()){
		m_vCP.push_back(idx);
		m_iNcp++;
		updateCPVAO();
	}
	else if(move){
		m_vX[idx] = pos;
		updateCPVAO();
	}
}
//! 固定点解除
void rxMeshDeform2D::UnsetCP(int idx)
{
	if(std::find(m_vCP.begin(), m_vCP.end(), idx) != m_vCP.end()){
		std::remove(m_vCP.begin(), m_vCP.end(), idx);
		m_iNcp--;
		updateCPVAO();
	}
}

/*!
 * OpenGLによるメッシュ,頂点,固定頂点描画
 */
void rxMeshDeform2D::DrawMesh(void)
{
	UpdateDataVAO(m_vao_mesh, 0, (GLfloat*)&m_vX[0], 2*m_iNv);
	glBindVertexArray(m_vao_mesh);
	glDrawElements(GL_TRIANGLES, m_vTri.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}
void rxMeshDeform2D::DrawPoints(void)
{
	UpdateDataVAO(m_vao_mesh, 0, (GLfloat*)&m_vX[0], 2*m_iNv);
	glBindVertexArray(m_vao_mesh);
	glDrawArrays(GL_POINTS, 0, m_iNv);
	glBindVertexArray(0);
}
void rxMeshDeform2D::DrawFixPoints(void)
{
	glBindVertexArray(m_vao_fix);
	glDrawArrays(GL_POINTS, 0, m_iNcp);
	glBindVertexArray(0);
}
void rxMeshDeform2D::InitVAO(void)
{
	// 頂点配列オブジェクトの作成
	if(m_vao_mesh != 0) glDeleteVertexArrays(1, &m_vao_mesh);
	m_vao_mesh = CreateVAO((GLfloat*)&m_vX[0], m_iNv, 2, &m_vTri[0], m_iNt, 0, 0, 0, 0, (GLfloat*)&m_vTC[0], m_iNv);
}



/*!
 * n×nの頂点を持つメッシュ生成(x-z平面)
 * @param[in] c1,c2 2端点座標
 * @param[in] nx,ny 格子分割数
 */
void rxMeshDeform2D::generateMesh(glm::vec2 c1, glm::vec2 c2, int nx, int ny)
{
	if(!m_vX.empty()){
		m_vX.clear();
		m_vTC.clear();
		m_vP.clear();
		m_vTri.clear();
	}

	// グリッド状に頂点座標生成
	// dx,dz:格子分割幅
	double dx = (c2[0]-c1[0])/(nx-1.0);
	double dz = (c2[1]-c1[1])/(ny-1.0);
	// 頂点情報格納配列の容量確保
	m_iNv = nx*ny;
	m_vX.resize(m_iNv);		// 頂点位置座標
	m_vP.resize(m_iNv);		// 頂点位置座標(初期座標として保存しておくための配列)
	m_vTC.resize(m_iNv);	// テクスチャ座標
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i < nx; ++i){
			// 格子頂点位置
			glm::vec2 pos;
			pos[0] = c1[0]+i*dx;
			pos[1] = c1[1]+j*dz;

			// 位置とテクスチャ座標を格納
			int idx = IDX(i, j, nx);
			m_vX[idx] = m_vP[idx] = pos;
			m_vTC[idx] = glm::vec2((pos[0]+c1[0])/(c2[0]-c1[0]), (pos[1]+c1[1])/(c2[1]-c1[1]));
		}
	}

	// メッシュ作成
	for(int j = 0; j < ny-1; ++j){
		for(int i = 0; i < nx-1; ++i){
			// 各グリッドセル(四角形)を2つの三角形メッシュに分割
			m_vTri.push_back(IDX(i, j, nx));
			m_vTri.push_back(IDX(i+1, j+1, nx));
			m_vTri.push_back(IDX(i+1, j, nx));

			m_vTri.push_back(IDX(i, j, nx));
			m_vTri.push_back(IDX(i, j+1, nx));
			m_vTri.push_back(IDX(i+1, j+1, nx));
		}
	}

	// 三角形メッシュ数
	m_iNt = (int)m_vTri.size()/3;
}

/*!
 * 線分(を含む直線)と点の距離
 * @param[in] v0,v1 線分の両端点座標
 * @param[in] p 点の座標
 * @return 距離
 */
inline double segment_point_dist(const glm::vec2 &v0, const glm::vec2 &v1, const glm::vec2 &p, glm::vec2 &ip)
{
	glm::vec2 v = glm::normalize(v1-v0);
	glm::vec2 vp = p-v0;
	glm::vec2 vh = glm::dot(vp, v)*v;
	ip = v0+vh;
	return glm::length(vp-vh);
}


/*!
 * n頂点のランダムメッシュ生成(x-z平面)
 * @param[in] c1,c2 2端点座標
 * @param[in] min_dist メッシュ頂点間の最低距離
 * @param[in] n メッシュ頂点数
 */
void rxMeshDeform2D::generateRandomMesh(glm::vec2 c1, glm::vec2 c2, double min_dist, int n)
{
	if(!m_vX.empty()){
		m_vX.clear();
		m_vTC.clear();
		m_vP.clear();
		m_vTri.clear();
	}

	glm::vec2 minp = c1, maxp = c2;

	// メッシュ生成範囲を配列に格納
	vector<glm::vec2> c(4);
	c[0] = glm::vec2(minp[0], minp[1]);
	c[1] = glm::vec2(maxp[0], minp[1]);
	c[2] = glm::vec2(maxp[0], maxp[1]);
	c[3] = glm::vec2(minp[0], maxp[1]);

	// 4隅の点を追加(メッシュ全体形状制御用)
	vector<glm::vec2> points;
	for(int i = 0; i < 4; ++i) points.push_back(c[i]);

	// 境界エッジ上にmin_distを基準に点を追加(メッシュ全体形状制御用)
	// - 全体の形状が四角なので4つの境界辺にそれぞれ順番にランダムな感覚で点を設置していく
	double d = 0.0;
	for(int j = 0; j < 4; ++j){
		glm::vec2 v0 = c[j];
		glm::vec2 v1 = c[(j == 3 ? 0 : j+1)];
		glm::vec2 edir = v1-v0;
		double len = glm::length(edir);
		while(d < len){
			// 境界辺上の距離をmin_dist*乱数分進める
			double t = min_dist*RX_RAND(1.0, 1.4);
			d += t;	// 4隅から境界辺上の距離
			if(len-d < 0.3*min_dist) break;	// 辺のもう一方の端点を超えたらループ終了
			points.push_back(v0+float(d)*glm::normalize(edir));	// 頂点座標を格納
		}
		d = 0;
	}

	// ポアソンディスクサンプリングで内部に点を生成
	rxUniformPoissonDiskSampler sampler(minp, maxp, n, 10, min_dist);
	sampler.Generate(points);

	// ドロネー三角形分割で三角形を生成
	vector< vector<int> > tris;
	CreateDelaunayTriangles(points, tris);

	// 頂点数と三角形数
	m_iNv = (int)points.size();
	m_iNt = (int)tris.size();

	// 頂点座標配列，テクスチャ座標配列の生成
	m_vX = m_vP = points;
	m_vTC.resize(m_iNv);
	for(int i = 0; i < m_iNv; ++i){
		m_vTC[i] = glm::vec2((points[i][0]+c1[0])/(c2[0]-c1[0]), (points[i][1]+c1[1])/(c2[1]-c1[1]));
	}
	// メッシュ情報配列の生成
	m_vTri.resize(m_iNt*3);
	for(int i = 0; i < m_iNt; ++i){
		m_vTri[3*i+0] = tris[i][0];
		m_vTri[3*i+1] = tris[i][1];
		m_vTri[3*i+2] = tris[i][2];
	}
}

