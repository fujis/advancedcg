/*!
  @file tetmesh.h

  @brief 表面三角形メッシュ(OBJなど)からの四面体メッシュ生成
		 - 外部ライブラリ(TetGenなど)を使わない簡易版
		 - 表面上の点(Poisson Disk Sampling)と内部の格子点(BCC格子)を3次元ドロネー分割(Bowyer-Watson法)し，
		   モデルの外側にある四面体や表面付近の扁平な四面体を除いて四面体メッシュを作る
		 - 四面体メッシュの表面は元の表面メッシュを頂点間隔程度の解像度で近似したものになる
		   (元の表面メッシュの頂点・三角形がそのまま残るわけではない)
		 - 依存: glm, STLのみ (OpenGLには非依存)

  @author Makoto Fujisawa
  @date 2026-10
*/

#ifndef _RX_TETMESH_H_
#define _RX_TETMESH_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <functional>

#include "glm/glm.hpp"


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
//! 四面体メッシュ生成パラメータ
struct rxTetMeshParams
{
	double spacing;			//!< 頂点間隔(正規化後のモデル座標での長さ)
	double min_quality;		//!< 表面にある四面体の品質(正四面体で1)がこれ未満なら除く
	bool use_surface_vertices;	//!< 元の表面メッシュの頂点を優先して表面点に使う(角を残したい場合など)

	rxTetMeshParams() : spacing(0.1), min_quality(0.05), use_surface_vertices(true) {}
};

/*!
 * 表面三角形メッシュから四面体メッシュを生成
 *  - 表面メッシュは閉じていることが望ましい(内外判定に使うため．多少の穴は3軸方向の判定の多数決で吸収する)
 * @param[in] vrts 表面メッシュの頂点座標
 * @param[in] tris 表面メッシュの三角形(3つで1ポリゴン)
 * @param[in] prm 生成パラメータ
 * @param[out] points 四面体メッシュの頂点座標
 * @param[out] surf 四面体メッシュの表面三角形(3つで1ポリゴン，外向きが表になるように反時計回り)
 * @param[out] tets 四面体(4つで1つの四面体)
 * @return 生成に成功したらtrue
 */
inline bool GenerateTetMesh(const std::vector<glm::vec3> &vrts, const std::vector<int> &tris, const rxTetMeshParams &prm,
							std::vector<glm::vec3> &points, std::vector<int> &surf, std::vector<int> &tets);

/*!
 * 四面体の品質(正四面体で1，つぶれるほど0に近づく)
 *  - 12*(3V)^(2/3)/Σl^2 (Vは体積，lは6辺の長さ)
 */
inline double TetQuality(const glm::dvec3 &p0, const glm::dvec3 &p1, const glm::dvec3 &p2, const glm::dvec3 &p3);



//=============================================================================
// 実装(ヘッダオンリー)
//=============================================================================
namespace rx_tetmesh
{

//! 向き(p0,p1,p2から見てp3が表側(反時計回りの法線方向)にあれば正，体積の6倍)
inline double orient(const glm::dvec3 &p0, const glm::dvec3 &p1, const glm::dvec3 &p2, const glm::dvec3 &p3)
{
	return glm::dot(glm::cross(p1-p0, p2-p0), p3-p0);
}


//-----------------------------------------------------------------------------
// 三角形メッシュによる内外判定
//  - 点からx,y,z軸の正方向にレイを飛ばして三角形との交差数の偶奇を調べ，3方向の多数決で判定する
//  - 三角形は各軸に垂直な平面上の2次元格子に登録しておき，レイが通るセルの三角形のみ調べる
//-----------------------------------------------------------------------------
class rxInsideTester
{
	struct AxisGrid
	{
		int n[2];							//!< 格子の分割数
		double minp[2], h[2];				//!< 格子の最小座標とセル幅
		std::vector< std::vector<int> > cells;	//!< 各セルに含まれる三角形
	};

	std::vector<glm::dvec3> m_v;		//!< 頂点
	std::vector<int> m_tri;				//!< 三角形
	AxisGrid m_grid[3];					//!< 各軸方向のレイ用の格子
	double m_eps;						//!< レイの位置をずらす量(三角形の辺や頂点をちょうど通るのを避ける)

public:
	void Set(const std::vector<glm::vec3> &vrts, const std::vector<int> &tris)
	{
		m_v.resize(vrts.size());
		for(size_t i = 0; i < vrts.size(); ++i) m_v[i] = glm::dvec3(vrts[i]);
		m_tri = tris;

		glm::dvec3 minp(1e30), maxp(-1e30);
		for(auto &v : m_v){ minp = glm::min(minp, v); maxp = glm::max(maxp, v); }
		glm::dvec3 ext = maxp-minp;
		m_eps = 1.0e-7*std::max(ext[0], std::max(ext[1], ext[2]));

		int nt = (int)m_tri.size()/3;
		int ndiv = std::max(1, std::min(256, (int)(2.0*sqrt((double)nt))));
		for(int a = 0; a < 3; ++a){
			AxisGrid &g = m_grid[a];
			for(int k = 0; k < 2; ++k){
				int b = (a+1+k)%3;
				g.n[k] = ndiv;
				g.minp[k] = minp[b];
				g.h[k] = std::max(ext[b], 1e-12)/ndiv;
			}
			g.cells.assign(g.n[0]*g.n[1], std::vector<int>());
			for(int t = 0; t < nt; ++t){
				int i0[2], i1[2];
				for(int k = 0; k < 2; ++k){
					int b = (a+1+k)%3;
					double lo = 1e30, hi = -1e30;
					for(int o = 0; o < 3; ++o){
						double x = m_v[m_tri[3*t+o]][b];
						lo = std::min(lo, x); hi = std::max(hi, x);
					}
					i0[k] = cellIndex(g, k, lo-m_eps);
					i1[k] = cellIndex(g, k, hi+m_eps);
				}
				for(int j = i0[1]; j <= i1[1]; ++j){
					for(int i = i0[0]; i <= i1[0]; ++i){
						g.cells[i+j*g.n[0]].push_back(t);
					}
				}
			}
		}
	}

	//! 点pがメッシュの内側にあるかどうか
	bool Inside(const glm::dvec3 &p) const
	{
		int votes = 0;
		for(int a = 0; a < 3; ++a){
			if(crossings(a, p)%2 == 1) votes++;
		}
		return votes >= 2;
	}

private:
	static int cellIndex(const AxisGrid &g, int k, double x)
	{
		int i = (int)floor((x-g.minp[k])/g.h[k]);
		return std::max(0, std::min(g.n[k]-1, i));
	}

	//! 点pから軸aの正方向に飛ばしたレイと三角形の交差数
	int crossings(int a, const glm::dvec3 &p) const
	{
		const AxisGrid &g = m_grid[a];
		int b = (a+1)%3, c = (a+2)%3;
		// 辺や頂点をちょうど通らないように少しずらす(ずらす量は軸ごとに変える)
		double px = p[b]+m_eps*(0.7548776662+0.13*a), py = p[c]+m_eps*(0.5698402910+0.17*a);
		if(px < g.minp[0] || py < g.minp[1] || px > g.minp[0]+g.n[0]*g.h[0] || py > g.minp[1]+g.n[1]*g.h[1]) return 0;
		const std::vector<int> &cell = g.cells[cellIndex(g, 0, px)+cellIndex(g, 1, py)*g.n[0]];
		int cnt = 0;
		for(int t : cell){
			const glm::dvec3 &v0 = m_v[m_tri[3*t]], &v1 = m_v[m_tri[3*t+1]], &v2 = m_v[m_tri[3*t+2]];
			double d = (v1[c]-v2[c])*(v0[b]-v2[b])+(v2[b]-v1[b])*(v0[c]-v2[c]);
			if(fabs(d) < 1e-300) continue;
			double l0 = ((v1[c]-v2[c])*(px-v2[b])+(v2[b]-v1[b])*(py-v2[c]))/d;
			double l1 = ((v2[c]-v0[c])*(px-v2[b])+(v0[b]-v2[b])*(py-v2[c]))/d;
			double l2 = 1.0-l0-l1;
			if(l0 < 0.0 || l1 < 0.0 || l2 < 0.0) continue;
			double z = l0*v0[a]+l1*v1[a]+l2*v2[a];
			if(z > p[a]) cnt++;
		}
		return cnt;
	}
};


//-----------------------------------------------------------------------------
// 点の近傍探索用のハッシュ格子
//-----------------------------------------------------------------------------
class rxPointHash
{
	double m_h;
	std::unordered_map<long long, std::vector<int> > m_cells;
	const std::vector<glm::dvec3> *m_pts;

	long long key(int i, int j, int k) const
	{
		return ((long long)(i+(1<<20))<<42)|((long long)(j+(1<<20))<<21)|(long long)(k+(1<<20));
	}
	void cell(const glm::dvec3 &p, int c[3]) const
	{
		for(int d = 0; d < 3; ++d) c[d] = (int)floor(p[d]/m_h);
	}

public:
	rxPointHash(double h, const std::vector<glm::dvec3> *pts) : m_h(h), m_pts(pts) {}

	void Add(int idx)
	{
		int c[3]; cell((*m_pts)[idx], c);
		m_cells[key(c[0], c[1], c[2])].push_back(idx);
	}

	//! pから距離r(<=格子幅)以内に点があるかどうか
	bool HasNeighbor(const glm::dvec3 &p, double r) const
	{
		int c[3]; cell(p, c);
		double r2 = r*r;
		for(int k = c[2]-1; k <= c[2]+1; ++k){
			for(int j = c[1]-1; j <= c[1]+1; ++j){
				for(int i = c[0]-1; i <= c[0]+1; ++i){
					auto it = m_cells.find(key(i, j, k));
					if(it == m_cells.end()) continue;
					for(int q : it->second){
						glm::dvec3 d = (*m_pts)[q]-p;
						if(glm::dot(d, d) < r2) return true;
					}
				}
			}
		}
		return false;
	}
};


//-----------------------------------------------------------------------------
// 3次元ドロネー分割(Bowyer-Watson法)
//-----------------------------------------------------------------------------
class rxDelaunay3D
{
	struct Tet
	{
		int v[4];			//!< 頂点(orient(v0,v1,v2,v3) > 0)
		int nb[4];			//!< 頂点v[i]の対面で隣接する四面体(なければ-1)
		glm::dvec3 cc;		//!< 外接球の中心
		double r2;			//!< 外接球の半径の2乗
		bool alive;
	};

	std::vector<glm::dvec3> m_p;	//!< 点(最後の4つはスーパー四面体の頂点)
	std::vector<Tet> m_t;			//!< 四面体
	std::vector<int> m_free;		//!< 再利用可能な四面体番号
	int m_last;						//!< 直前に作った四面体(点の探索の開始位置)

	// 作業用
	std::vector<int> m_mark;		//!< 空洞(cavity)に含まれるかどうかのマーク(点番号+1を入れる)

public:
	/*!
	 * ドロネー分割
	 * @param[in] pts 点(重複がないこと．格子状などの退化を避けるため少しずらしておくとよい)
	 * @param[out] tets 四面体(4つで1つ．orient > 0)
	 */
	void Triangulate(const std::vector<glm::dvec3> &pts, std::vector<int> &tets)
	{
		int n = (int)pts.size();
		tets.clear();
		if(n < 4) return;

		// スーパー四面体(全点を十分大きく囲む四面体)
		glm::dvec3 minp(1e30), maxp(-1e30);
		for(auto &p : pts){ minp = glm::min(minp, p); maxp = glm::max(maxp, p); }
		glm::dvec3 cen = 0.5*(minp+maxp);
		double R = 0.5*glm::length(maxp-minp)+1e-6;
		double s = 100.0*R;
		m_p = pts;
		m_p.push_back(cen+s*glm::dvec3(1, 1, 1));
		m_p.push_back(cen+s*glm::dvec3(-1, -1, 1));
		m_p.push_back(cen+s*glm::dvec3(-1, 1, -1));
		m_p.push_back(cen+s*glm::dvec3(1, -1, -1));

		m_t.clear(); m_free.clear();
		m_mark.clear();
		int t0 = newTet(n, n+1, n+2, n+3);
		for(int i = 0; i < 4; ++i) m_t[t0].nb[i] = -1;
		m_last = t0;

		// 点を1つずつ挿入
		std::vector<int> order(n);
		for(int i = 0; i < n; ++i) order[i] = i;
		std::mt19937 rng(12345);
		std::shuffle(order.begin(), order.end(), rng);
		for(int i : order) insert(i);

		// スーパー四面体の頂点を含まない四面体を出力
		for(auto &t : m_t){
			if(!t.alive) continue;
			if(t.v[0] >= n || t.v[1] >= n || t.v[2] >= n || t.v[3] >= n) continue;
			for(int k = 0; k < 4; ++k) tets.push_back(t.v[k]);
		}
	}

private:
	int newTet(int a, int b, int c, int d)
	{
		int idx;
		if(!m_free.empty()){ idx = m_free.back(); m_free.pop_back(); }
		else{ idx = (int)m_t.size(); m_t.push_back(Tet()); m_mark.push_back(0); }
		Tet &t = m_t[idx];
		t.v[0] = a; t.v[1] = b; t.v[2] = c; t.v[3] = d;
		t.nb[0] = t.nb[1] = t.nb[2] = t.nb[3] = -1;
		t.alive = true;
		m_mark[idx] = 0;
		calCircumsphere(t);
		return idx;
	}

	void calCircumsphere(Tet &t)
	{
		const glm::dvec3 &p0 = m_p[t.v[0]];
		glm::dvec3 a = m_p[t.v[1]]-p0, b = m_p[t.v[2]]-p0, c = m_p[t.v[3]]-p0;
		double d = 2.0*glm::dot(a, glm::cross(b, c));
		if(fabs(d) < 1e-300){	// 退化した四面体(通常は生じない)
			t.cc = 0.25*(p0+m_p[t.v[1]]+m_p[t.v[2]]+m_p[t.v[3]]);
			t.r2 = 1e300;
			return;
		}
		glm::dvec3 o = (glm::dot(a, a)*glm::cross(b, c)+glm::dot(b, b)*glm::cross(c, a)+glm::dot(c, c)*glm::cross(a, b))/d;
		t.cc = p0+o;
		t.r2 = glm::dot(o, o);
	}

	bool inSphere(int t, const glm::dvec3 &p) const
	{
		glm::dvec3 d = p-m_t[t].cc;
		return glm::dot(d, d) < m_t[t].r2;
	}

	//! 四面体tの頂点iの対面の3頂点(tの外側から見て反時計回り)
	void face(int t, int i, int f[3]) const
	{
		static const int fi[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
		for(int k = 0; k < 3; ++k) f[k] = m_t[t].v[fi[i][k]];
	}

	//! 点pを外接球に含む四面体を探す(四面体をたどって点を含む四面体を探す)
	int locate(const glm::dvec3 &p)
	{
		int t = m_last;
		if(t < 0 || !m_t[t].alive){
			for(t = (int)m_t.size()-1; t >= 0 && !m_t[t].alive; --t);
		}
		int steps = 0, maxsteps = (int)m_t.size()+100;
		while(steps++ < maxsteps){
			int next = -1;
			for(int i = 0; i < 4; ++i){
				int f[3]; face(t, i, f);
				// 面の外側(四面体と反対側)にpがあればその隣へ
				if(orient(m_p[f[0]], m_p[f[1]], m_p[f[2]], p) > 0.0 && m_t[t].nb[i] >= 0){
					next = m_t[t].nb[i];
					break;
				}
			}
			if(next < 0) break;
			t = next;
		}
		if(inSphere(t, p)) return t;

		// 見つからなかった場合は全探索
		for(int k = 0; k < (int)m_t.size(); ++k){
			if(m_t[k].alive && inSphere(k, p)) return k;
		}
		return -1;
	}

	void insert(int pi)
	{
		const glm::dvec3 &p = m_p[pi];
		int t0 = locate(p);
		if(t0 < 0) return;

		// 外接球にpを含む四面体の集合(空洞)を探索
		int stamp = pi+1;
		std::vector<int> cavity, stack;
		stack.push_back(t0);
		m_mark[t0] = stamp;
		while(!stack.empty()){
			int t = stack.back(); stack.pop_back();
			cavity.push_back(t);
			for(int i = 0; i < 4; ++i){
				int nb = m_t[t].nb[i];
				if(nb < 0 || m_mark[nb] == stamp) continue;
				if(inSphere(nb, p)){
					m_mark[nb] = stamp;
					stack.push_back(nb);
				}
			}
		}

		// 空洞の境界面を集める
		//  - 数値誤差で空洞が星型にならない(境界面からpが見えない)場合は外側の四面体も空洞に加えてやり直す
		struct BFace { int f[3]; int outside; int from; };
		std::vector<BFace> bfaces;
		bool ok = false;
		for(int iter = 0; iter < 100 && !ok; ++iter){
			ok = true;
			bfaces.clear();
			for(size_t c = 0; c < cavity.size() && ok; ++c){
				int t = cavity[c];
				for(int i = 0; i < 4; ++i){
					int nb = m_t[t].nb[i];
					if(nb >= 0 && m_mark[nb] == stamp) continue;
					BFace bf;
					face(t, i, bf.f);
					bf.outside = nb; bf.from = t;
					// 新しい四面体(f0,f2,f1,p)の向きが正でなければ空洞を広げる
					if(orient(m_p[bf.f[0]], m_p[bf.f[2]], m_p[bf.f[1]], p) <= 0.0 && nb >= 0){
						m_mark[nb] = stamp;
						cavity.push_back(nb);
						ok = false;
						break;
					}
					bfaces.push_back(bf);
				}
			}
		}

		// 空洞の四面体を削除
		//  - 番号の再利用(m_free)は新しい四面体を作り終えてから
		//    (外側の四面体の隣接情報を空洞の四面体番号で探して書き換えるため)
		for(int t : cavity) m_t[t].alive = false;

		// 境界面とpを結んで新しい四面体を作る
		std::unordered_map<long long, std::pair<int, int> > edges;	// 辺 -> (四面体, 局所番号)
		long long np = (long long)m_p.size();
		for(const BFace &bf : bfaces){
			int nt = newTet(bf.f[0], bf.f[2], bf.f[1], pi);	// pから見て境界面が反時計回り
			m_mark[nt] = 0;
			// pの対面は境界面 -> 外側の四面体と隣接
			m_t[nt].nb[3] = bf.outside;
			if(bf.outside >= 0){
				for(int i = 0; i < 4; ++i){
					if(m_t[bf.outside].nb[i] == bf.from) m_t[bf.outside].nb[i] = nt;
				}
			}
			// 他の3面(pを含む面)は辺を共有する新しい四面体と隣接
			for(int i = 0; i < 3; ++i){
				int a = m_t[nt].v[(i+1)%3], b = m_t[nt].v[(i+2)%3];
				long long key = std::min(a, b)*np+std::max(a, b);
				auto it = edges.find(key);
				if(it == edges.end()){
					edges[key] = std::make_pair(nt, i);
				}
				else{
					m_t[nt].nb[i] = it->second.first;
					m_t[it->second.first].nb[it->second.second] = nt;
					edges.erase(it);
				}
			}
			m_last = nt;
		}
		for(int t : cavity) m_free.push_back(t);
	}
};

} // namespace rx_tetmesh


/*!
 * 四面体の品質
 */
inline double TetQuality(const glm::dvec3 &p0, const glm::dvec3 &p1, const glm::dvec3 &p2, const glm::dvec3 &p3)
{
	double V = fabs(rx_tetmesh::orient(p0, p1, p2, p3))/6.0;
	glm::dvec3 e[6] = { p1-p0, p2-p0, p3-p0, p2-p1, p3-p1, p3-p2 };
	double l2 = 0.0;
	for(int i = 0; i < 6; ++i) l2 += glm::dot(e[i], e[i]);
	if(l2 < 1e-300) return 0.0;
	return 12.0*pow(3.0*V, 2.0/3.0)/l2;
}


/*!
 * 表面三角形メッシュから四面体メッシュを生成
 */
inline bool GenerateTetMesh(const std::vector<glm::vec3> &vrts, const std::vector<int> &tris, const rxTetMeshParams &prm,
							std::vector<glm::vec3> &points, std::vector<int> &surf, std::vector<int> &tets)
{
	using namespace rx_tetmesh;

	int nv = (int)vrts.size(), nt = (int)tris.size()/3;
	if(nv < 4 || nt < 4 || prm.spacing <= 0.0) return false;
	double h = prm.spacing;
	std::mt19937 rng(1234);
	std::uniform_real_distribution<double> rnd(0.0, 1.0);

	// 内外判定
	rxInsideTester inside;
	inside.Set(vrts, tris);

	glm::dvec3 minp(1e30), maxp(-1e30);
	for(auto &v : vrts){ minp = glm::min(minp, glm::dvec3(v)); maxp = glm::max(maxp, glm::dvec3(v)); }

	// 1. 表面上の点(Poisson Disk Sampling: 候補点のうち既存の点から距離r以上離れたものだけを採用)
	//    - 候補点は元メッシュの頂点(use_surface_vertices=trueの場合)と三角形上のランダムな点
	double r = 0.9*h;
	std::vector<glm::dvec3> pts;
	rxPointHash hash(r, &pts);
	auto trySurfacePoint = [&](const glm::dvec3 &p){
		if(hash.HasNeighbor(p, r)) return;
		pts.push_back(p);
		hash.Add((int)pts.size()-1);
	};
	if(prm.use_surface_vertices){
		for(auto &v : vrts) trySurfacePoint(glm::dvec3(v));
	}
	std::vector<glm::dvec3> cand;
	for(int t = 0; t < nt; ++t){
		glm::dvec3 v0(vrts[tris[3*t]]), v1(vrts[tris[3*t+1]]), v2(vrts[tris[3*t+2]]);
		double area = 0.5*glm::length(glm::cross(v1-v0, v2-v0));
		double m = 10.0*area/(r*r);	// 1点あたりの面積がr^2/10程度になるように候補点を置く
		int k = (int)m+(rnd(rng) < m-(int)m ? 1 : 0);
		for(int i = 0; i < k; ++i){
			double a = sqrt(rnd(rng)), b = rnd(rng);
			cand.push_back((1.0-a)*v0+a*(1.0-b)*v1+a*b*v2);
		}
	}
	std::shuffle(cand.begin(), cand.end(), rng);
	for(auto &p : cand) trySurfacePoint(p);
	int nsurf = (int)pts.size();

	// 2. 内部の点(BCC格子点のうち，内部にあって表面の点から離れたもの)
	//    - 格子定数aは最近接の格子点間の距離(a√3/2)が頂点間隔hになるように設定
	//    - 格子点はドロネー分割の退化を避けるために少しずらす
	double a = 2.0*h/sqrt(3.0);
	for(int sub = 0; sub < 2; ++sub){
		double o = 0.5*a*sub;
		for(double z = minp[2]+o; z <= maxp[2]; z += a){
			for(double y = minp[1]+o; y <= maxp[1]; y += a){
				for(double x = minp[0]+o; x <= maxp[0]; x += a){
					glm::dvec3 p(x, y, z);
					p += 1.0e-3*h*glm::dvec3(rnd(rng)-0.5, rnd(rng)-0.5, rnd(rng)-0.5);
					if(!inside.Inside(p)) continue;
					if(hash.HasNeighbor(p, 0.7*h)) continue;
					pts.push_back(p);
				}
			}
		}
	}
	if((int)pts.size() < 4) return false;

	// 3. ドロネー分割
	std::vector<int> dt;
	rxDelaunay3D delaunay;
	delaunay.Triangulate(pts, dt);
	int ntet = (int)dt.size()/4;

	// 4. モデルの外側にある四面体を除く
	//    - 重心と，各面・各辺の中点を重心側に少し寄せた点がすべて内側にあるものを残す
	//      (凹部をまたぐ四面体を除くため．面や辺の中点そのものは平らな表面上では内外判定が不安定なので寄せる)
	std::vector<char> keep(ntet, 0);
	for(int t = 0; t < ntet; ++t){
		const int *v = &dt[4*t];
		glm::dvec3 p[4] = { pts[v[0]], pts[v[1]], pts[v[2]], pts[v[3]] };
		glm::dvec3 c = 0.25*(p[0]+p[1]+p[2]+p[3]);
		if(!inside.Inside(c)) continue;
		bool ok = true;
		bool all_surf = (v[0] < nsurf && v[1] < nsurf && v[2] < nsurf && v[3] < nsurf);
		if(all_surf || v[0] < nsurf || v[1] < nsurf || v[2] < nsurf || v[3] < nsurf){
			for(int i = 0; i < 4 && ok; ++i){
				glm::dvec3 fc = (p[0]+p[1]+p[2]+p[3]-p[i])/3.0;
				if(!inside.Inside(0.8*fc+0.2*c)) ok = false;
			}
			for(int i = 0; i < 4 && ok; ++i){
				for(int j = i+1; j < 4 && ok; ++j){
					if(!(v[i] < nsurf && v[j] < nsurf)) continue;
					glm::dvec3 m = 0.5*(p[i]+p[j]);
					if(!inside.Inside(0.8*m+0.2*c)) ok = false;
				}
			}
		}
		keep[t] = ok ? 1 : 0;
	}

	// 面 -> 四面体の対応(隣接関係の計算用)
	auto faceKey = [&](int a, int b, int c){
		int s[3] = { a, b, c };
		std::sort(s, s+3);
		long long n = (long long)pts.size();
		return ((long long)s[0]*n+s[1])*n+s[2];
	};
	static const int fi[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };	// 頂点iの対面(外から見て反時計回り)
	std::vector<int> nbr(4*ntet, -1);
	{
		std::unordered_map<long long, int> fmap;
		for(int t = 0; t < ntet; ++t){
			for(int i = 0; i < 4; ++i){
				long long key = faceKey(dt[4*t+fi[i][0]], dt[4*t+fi[i][1]], dt[4*t+fi[i][2]]);
				auto it = fmap.find(key);
				if(it == fmap.end()) fmap[key] = 4*t+i;
				else{ nbr[4*t+i] = it->second/4; nbr[it->second] = t; fmap.erase(it); }
			}
		}
	}
	auto isBoundary = [&](int t, int i){ int n = nbr[4*t+i]; return n < 0 || !keep[n]; };

	// 5. 表面にある扁平な四面体(sliver)を除く(除くと新たに表面に出る四面体もあるので数回繰り返す)
	for(int it = 0; it < 5; ++it){
		int nrem = 0;
		for(int t = 0; t < ntet; ++t){
			if(!keep[t]) continue;
			int nb = 0;
			for(int i = 0; i < 4; ++i) if(isBoundary(t, i)) nb++;
			if(nb == 0) continue;
			const int *v = &dt[4*t];
			double q = TetQuality(pts[v[0]], pts[v[1]], pts[v[2]], pts[v[3]]);
			// 3面以上が表面に出ている四面体は品質が低めでも除く(表面のとげ)
			if(q < prm.min_quality || (nb >= 3 && q < 4.0*prm.min_quality)){
				keep[t] = 0; nrem++;
			}
		}
		if(nrem == 0) break;
	}

	// 6. 表面の非多様体の辺(3つ以上の表面三角形が共有する辺)を解消する
	//    - その辺を含む表面三角形を持つ四面体のうち品質が最も低いものを除く
	for(int it = 0; it < 20; ++it){
		std::unordered_map<long long, std::vector<int> > emap;	// 辺 -> その辺を含む表面三角形を持つ四面体
		long long n = (long long)pts.size();
		for(int t = 0; t < ntet; ++t){
			if(!keep[t]) continue;
			for(int i = 0; i < 4; ++i){
				if(!isBoundary(t, i)) continue;
				for(int k = 0; k < 3; ++k){
					int a = dt[4*t+fi[i][k]], b = dt[4*t+fi[i][(k+1)%3]];
					emap[std::min(a, b)*n+std::max(a, b)].push_back(t);
				}
			}
		}
		int nrem = 0;
		for(auto &e : emap){
			if(e.second.size() <= 2) continue;
			int worst = -1; double qmin = 1e30;
			for(int t : e.second){
				if(!keep[t]) continue;
				const int *v = &dt[4*t];
				double q = TetQuality(pts[v[0]], pts[v[1]], pts[v[2]], pts[v[3]]);
				if(q < qmin){ qmin = q; worst = t; }
			}
			if(worst >= 0){ keep[worst] = 0; nrem++; }
		}
		if(nrem == 0) break;
	}

	// 7. 最大の連結成分(面で接する四面体)のみ残す
	{
		std::vector<int> comp(ntet, -1);
		int best = -1, bestsize = 0;
		for(int s = 0; s < ntet; ++s){
			if(!keep[s] || comp[s] >= 0) continue;
			int size = 0;
			std::vector<int> stack(1, s);
			comp[s] = s;
			while(!stack.empty()){
				int t = stack.back(); stack.pop_back();
				size++;
				for(int i = 0; i < 4; ++i){
					int n = nbr[4*t+i];
					if(n >= 0 && keep[n] && comp[n] < 0){ comp[n] = s; stack.push_back(n); }
				}
			}
			if(size > bestsize){ bestsize = size; best = s; }
		}
		for(int t = 0; t < ntet; ++t) if(keep[t] && comp[t] != best) keep[t] = 0;
	}

	// 8. 出力(使われている頂点のみ詰めて出力，表面は四面体の境界面)
	std::vector<int> newidx(pts.size(), -1);
	points.clear(); surf.clear(); tets.clear();
	auto vid = [&](int v){
		if(newidx[v] < 0){ newidx[v] = (int)points.size(); points.push_back(glm::vec3(pts[v])); }
		return newidx[v];
	};
	for(int t = 0; t < ntet; ++t){
		if(!keep[t]) continue;
		for(int i = 0; i < 4; ++i) tets.push_back(vid(dt[4*t+i]));
		for(int i = 0; i < 4; ++i){
			if(!isBoundary(t, i)) continue;
			for(int k = 0; k < 3; ++k) surf.push_back(vid(dt[4*t+fi[i][k]]));
		}
	}

	return !tets.empty();
}


#endif // #ifdef _RX_TETMESH_H_
