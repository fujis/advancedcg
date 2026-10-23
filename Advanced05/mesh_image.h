/*!
  @file mesh_image.h

  @brief 画像中のオブジェクト形状に沿った2D三角形メッシュの生成
		 - 背景(白)以外の画素をオブジェクトとし，その輪郭に沿った点とPoisson Disk Samplingによる内部点を
		   ドロネー三角形分割してメッシュを作る
		 - 依存: glm, rx_sampler.h, rx_delaunay.h (OpenGLには非依存)

  @author Makoto Fujisawa
  @date 2026-10
*/

#ifndef _RX_MESH_IMAGE_H_
#define _RX_MESH_IMAGE_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <algorithm>
#include <list>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <functional>

#include "glm/glm.hpp"
#include "glm/gtc/random.hpp"	// rx_sampler.hで使うglm::linearRand
#include "rx_sampler.h"
#include "rx_delaunay.h"


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
//! メッシュ生成パラメータ
struct rxImageMeshParams
{
	double spacing;			//!< 頂点間隔(ワールド座標)
	int white_threshold;	//!< RGBすべてがこの値以上の画素を背景(白)とみなす
	int blur;				//!< 輪郭抽出前のマスクのぼかし回数(輪郭のギザギザを抑える)
	int smooth_iter;		//!< 内部頂点の平滑化(+再三角形分割)の反復回数

	rxImageMeshParams() : spacing(0.07), white_threshold(240), blur(2), smooth_iter(5) {}
};

//! 生成されたメッシュの品質(三角形の最小内角)
struct rxMeshQuality
{
	double min_angle;		//!< 全三角形での最小内角[deg]
	double avg_min_angle;	//!< 各三角形の最小内角の平均[deg] (正三角形なら60)
	int num_bad;			//!< 最小内角が20度未満の三角形数
};


/*!
 * 画像からオブジェクト形状に沿ったメッシュを生成
 * @param[in] img 画像データ(1画素c byte, 1行wstep byte, 1行目が画像の下端)
 * @param[in] w,h,c,wstep 画像の幅,高さ,チャンネル数,1行のバイト数
 * @param[in] c1,c2 メッシュを配置する領域(画像の縦横比を保ってこの中に収める)
 * @param[in] prm 生成パラメータ
 * @param[out] points 頂点座標
 * @param[out] texcoords 頂点のテクスチャ座標(画像全体が[0,1]x[0,1])
 * @param[out] tris 三角形の頂点インデックス(3つで1ポリゴン, 反時計回り)
 * @return 生成に成功したらtrue
 */
inline bool GenerateImageMesh(const unsigned char* img, int w, int h, int c, int wstep,
					   glm::vec2 c1, glm::vec2 c2, const rxImageMeshParams& prm,
					   std::vector<glm::vec2>& points, std::vector<glm::vec2>& texcoords, std::vector<unsigned int>& tris);

/*!
 * 三角形メッシュの品質(内角)の計算
 */
inline rxMeshQuality CalMeshQuality(const std::vector<glm::vec2>& points, const std::vector<unsigned int>& tris);



//=============================================================================
// 実装(ヘッダオンリー)
//=============================================================================
//-----------------------------------------------------------------------------
// 内部で使う関数(名前の衝突を避けるため名前空間に入れる)
//-----------------------------------------------------------------------------
namespace rx_imgmesh
{

//-----------------------------------------------------------------------------
// 画像から作るスカラー場(オブジェクト内部で1,背景で0)
//  - 画像の周囲に1画素分の余白(値0)を付けて，輪郭が必ず閉曲線になるようにする
//  - 格子点(I,J)は画素(I-1,J-1)の中心に対応
//-----------------------------------------------------------------------------
struct rxImageField
{
	int W, H;					//!< 格子点数(画像サイズ+2)
	std::vector<float> f;		//!< 格子点の値
	glm::vec2 org;				//!< 画素(0,0)の左下隅のワールド座標
	double s;					//!< 1画素の大きさ(ワールド座標)

	float& at(int I, int J){ return f[I+J*W]; }
	float at(int I, int J) const { return f[I+J*W]; }

	//! 格子点(I,J)のワールド座標
	glm::vec2 pos(double I, double J) const { return org+glm::vec2(float(s*(I-0.5)), float(s*(J-0.5))); }

	//! ワールド座標pでの値(双線形補間)
	float value(const glm::vec2& p) const
	{
		double x = (p[0]-org[0])/s+0.5, y = (p[1]-org[1])/s+0.5;
		int I = (int)floor(x), J = (int)floor(y);
		if(I < 0 || J < 0 || I >= W-1 || J >= H-1) return 0.0f;
		float a = float(x-I), b = float(y-J);
		return (1-a)*(1-b)*at(I, J)+a*(1-b)*at(I+1, J)+(1-a)*b*at(I, J+1)+a*b*at(I+1, J+1);
	}

	//! ワールド座標pがオブジェクト内部かどうか
	bool inside(const glm::vec2& p) const { return value(p) >= 0.5f; }
};


/*!
 * 2値マスクの連結成分ラベリング(4近傍)
 * @param[in] m マスク(W×H)
 * @param[in] val ラベリング対象の値
 * @param[out] label 各画素のラベル(対象外は-1)
 * @param[out] area 各ラベルの画素数
 */
inline void labeling(const std::vector<char>& m, int W, int H, char val, std::vector<int>& label, std::vector<int>& area)
{
	label.assign(W*H, -1);
	area.clear();
	std::vector<int> stack;
	for(int s = 0; s < W*H; ++s){
		if(m[s] != val || label[s] >= 0) continue;
		int l = (int)area.size();
		area.push_back(0);
		stack.push_back(s);
		label[s] = l;
		while(!stack.empty()){
			int k = stack.back(); stack.pop_back();
			area[l]++;
			int x = k%W, y = k/W;
			int nb[4][2] = { { x-1, y }, { x+1, y }, { x, y-1 }, { x, y+1 } };
			for(auto& n : nb){
				if(n[0] < 0 || n[1] < 0 || n[0] >= W || n[1] >= H) continue;
				int k2 = n[0]+n[1]*W;
				if(m[k2] == val && label[k2] < 0){
					label[k2] = l;
					stack.push_back(k2);
				}
			}
		}
	}
}

/*!
 * 画像からスカラー場を作成
 *  - 白以外の画素をオブジェクトとし，最大の連結成分のみ残す
 *  - オブジェクト内部の小さな穴(背景画素の塊)は埋める
 *  - 輪郭を滑らかにするためにぼかす
 */
inline bool makeField(const unsigned char* img, int w, int h, int c, int wstep, glm::vec2 c1, glm::vec2 c2,
					  const rxImageMeshParams& prm, rxImageField& fld)
{
	int W = w+2, H = h+2;
	std::vector<char> m(W*H, 0);
	for(int j = 0; j < h; ++j){
		for(int i = 0; i < w; ++i){
			const unsigned char* p = img+j*wstep+i*c;
			bool white = true;
			int nc = (c >= 3) ? 3 : 1;	// 色のチャンネル数(グレースケールなら1)
			for(int k = 0; k < nc; ++k) if(p[k] < prm.white_threshold) white = false;
			if((c == 4 || c == 2) && p[c-1] < 128) white = true;	// アルファ付きなら透明部分も背景
			m[(i+1)+(j+1)*W] = white ? 0 : 1;
		}
	}

	// オブジェクトの最大連結成分のみ残す
	std::vector<int> label, area;
	labeling(m, W, H, 1, label, area);
	if(area.empty()) return false;
	int lmax = (int)(std::max_element(area.begin(), area.end())-area.begin());
	for(int k = 0; k < W*H; ++k) m[k] = (label[k] == lmax) ? 1 : 0;

	// 画素サイズとワールド座標の対応(縦横比を保って領域c1-c2の中央に配置)
	glm::vec2 ext = c2-c1;
	fld.s = std::min(ext[0]/w, ext[1]/h);
	fld.org = 0.5f*(c1+c2)-0.5f*float(fld.s)*glm::vec2(w, h);

	// 小さな穴(画像の外周につながらない背景成分で面積が頂点間隔^2の数倍以下)を埋める
	double sp_px = prm.spacing/fld.s;
	labeling(m, W, H, 0, label, area);
	int lout = label[0];	// 外周(余白)の背景成分
	for(int k = 0; k < W*H; ++k){
		if(label[k] >= 0 && label[k] != lout && area[label[k]] < 4.0*sp_px*sp_px) m[k] = 1;
	}

	// ぼかし(3x3の平均)
	fld.W = W; fld.H = H;
	fld.f.resize(W*H);
	for(int k = 0; k < W*H; ++k) fld.f[k] = m[k];
	std::vector<float> tmp(W*H);
	for(int b = 0; b < prm.blur; ++b){
		for(int J = 0; J < H; ++J){
			for(int I = 0; I < W; ++I){
				float sum = 0.0f; int n = 0;
				for(int y = std::max(0, J-1); y <= std::min(H-1, J+1); ++y){
					for(int x = std::max(0, I-1); x <= std::min(W-1, I+1); ++x){
						sum += fld.at(x, y); n++;
					}
				}
				tmp[I+J*W] = sum/n;
			}
		}
		fld.f.swap(tmp);
	}
	// 外周は必ず0(輪郭を閉曲線にするため)
	for(int I = 0; I < W; ++I) fld.at(I, 0) = fld.at(I, H-1) = 0.0f;
	for(int J = 0; J < H; ++J) fld.at(0, J) = fld.at(W-1, J) = 0.0f;

	return true;
}


//-----------------------------------------------------------------------------
// 輪郭抽出(marching squares)
//-----------------------------------------------------------------------------
/*!
 * スカラー場の等値線(値0.5)を閉曲線のリストとして抽出
 * @param[out] loops 閉曲線(ワールド座標の点列)
 */
inline void extractContours(const rxImageField& fld, std::vector< std::vector<glm::vec2> >& loops)
{
	const float iso = 0.5f;
	int W = fld.W, H = fld.H;

	// 格子のエッジID : 水平エッジ(I,J)-(I+1,J) -> 2(I+JW), 垂直エッジ(I,J)-(I,J+1) -> 2(I+JW)+1
	auto hid = [&](int I, int J){ return 2LL*(I+(long long)J*W); };
	auto vid = [&](int I, int J){ return 2LL*(I+(long long)J*W)+1; };

	std::unordered_map<long long, glm::vec2> epos;					// エッジ上の交点座標
	std::unordered_map<long long, std::vector<long long> > adj;		// 交点のつながり

	auto crossing = [&](long long id) -> void {
		if(epos.count(id)) return;
		long long k = id/2;
		int I = (int)(k%W), J = (int)(k/W);
		int I2 = (id%2 == 0) ? I+1 : I, J2 = (id%2 == 0) ? J : J+1;
		float a = fld.at(I, J), b = fld.at(I2, J2);
		double t = (iso-a)/(b-a);
		epos[id] = fld.pos(I+t*(I2-I), J+t*(J2-J));
	};
	auto segment = [&](long long e0, long long e1){
		crossing(e0); crossing(e1);
		adj[e0].push_back(e1);
		adj[e1].push_back(e0);
	};

	for(int J = 0; J < H-1; ++J){
		for(int I = 0; I < W-1; ++I){
			float v0 = fld.at(I, J), v1 = fld.at(I+1, J), v2 = fld.at(I+1, J+1), v3 = fld.at(I, J+1);
			int cs = (v0 >= iso ? 1 : 0) | (v1 >= iso ? 2 : 0) | (v2 >= iso ? 4 : 0) | (v3 >= iso ? 8 : 0);
			if(cs == 0 || cs == 15) continue;
			long long e0 = hid(I, J), e1 = vid(I+1, J), e2 = hid(I, J+1), e3 = vid(I, J);	// 下,右,上,左
			bool cen = 0.25f*(v0+v1+v2+v3) >= iso;
			switch(cs){
			case 1: case 14: segment(e3, e0); break;
			case 2: case 13: segment(e0, e1); break;
			case 3: case 12: segment(e3, e1); break;
			case 4: case 11: segment(e1, e2); break;
			case 6: case 9:  segment(e0, e2); break;
			case 7: case 8:  segment(e3, e2); break;
			case 5:	// 鞍点(v0,v2が内部)
				if(cen){ segment(e0, e1); segment(e2, e3); }
				else{ segment(e3, e0); segment(e1, e2); }
				break;
			case 10:	// 鞍点(v1,v3が内部)
				if(cen){ segment(e3, e0); segment(e1, e2); }
				else{ segment(e0, e1); segment(e2, e3); }
				break;
			}
		}
	}

	// 交点をたどって閉曲線にする
	std::unordered_set<long long> visited;
	for(auto& a : adj){
		if(visited.count(a.first)) continue;
		std::vector<glm::vec2> loop;
		long long prev = -1, cur = a.first;
		while(cur >= 0 && !visited.count(cur)){
			visited.insert(cur);
			loop.push_back(epos[cur]);
			long long next = -1;
			for(long long n : adj[cur]){
				if(n != prev && !visited.count(n)){ next = n; break; }
			}
			prev = cur;
			cur = next;
		}
		if(loop.size() >= 3) loops.push_back(loop);
	}
}

/*!
 * 閉曲線を弧長に沿って等間隔に再サンプリング
 * @param[in] loop 閉曲線
 * @param[in] sp 間隔
 * @param[out] out 再サンプリングされた点列(短すぎる閉曲線では空)
 */
inline void resampleLoop(const std::vector<glm::vec2>& loop, double sp, std::vector<glm::vec2>& out)
{
	out.clear();
	int n = (int)loop.size();
	std::vector<double> len(n+1, 0.0);
	for(int i = 0; i < n; ++i) len[i+1] = len[i]+glm::length(loop[(i+1)%n]-loop[i]);
	double L = len[n];
	if(L < 3.0*sp) return;

	int m = std::max(3, (int)(L/sp+0.5));
	double dl = L/m;
	int seg = 0;
	for(int k = 0; k < m; ++k){
		double l = k*dl;
		while(seg < n-1 && len[seg+1] < l) seg++;
		double t = (len[seg+1] > len[seg]) ? (l-len[seg])/(len[seg+1]-len[seg]) : 0.0;
		out.push_back(loop[seg]+float(t)*(loop[(seg+1)%n]-loop[seg]));
	}
}


//-----------------------------------------------------------------------------
// 三角形分割
//-----------------------------------------------------------------------------
/*!
 * ドロネー三角形分割して背景にはみ出した三角形を除く
 * @param[in] pts 頂点座標
 * @param[in] bedges 輪郭上で隣り合う頂点の組(このエッジは中点の内外判定を行わない)
 * @param[out] tris 三角形(3つで1ポリゴン)
 */
inline void triangulate(std::vector<glm::vec2>& pts, const std::unordered_set<long long>& bedges,
						const rxImageField& fld, std::vector<unsigned int>& tris)
{
	long long np = (long long)pts.size();
	std::vector< std::vector<int> > dt;
	CreateDelaunayTriangles(pts, dt);

	tris.clear();
	for(auto& t : dt){
		glm::vec2 p0 = pts[t[0]], p1 = pts[t[1]], p2 = pts[t[2]];

		// 重心が背景なら除く
		if(!fld.inside((p0+p1+p2)/3.0f)) continue;

		// 輪郭上のエッジ以外で，中点が背景にあるエッジを持つ三角形を除く(凹部をまたぐ三角形)
		bool ok = true;
		for(int o = 0; o < 3 && ok; ++o){
			int a = t[o], b = t[(o+1)%3];
			if(bedges.count(std::min(a, b)*np+std::max(a, b))) continue;
			if(!fld.inside(0.5f*(pts[a]+pts[b]))) ok = false;
		}
		if(!ok) continue;

		// 反時計回りにそろえる
		double cr = (p1[0]-p0[0])*(p2[1]-p0[1])-(p1[1]-p0[1])*(p2[0]-p0[0]);
		if(fabs(cr) < 1e-12) continue;
		if(cr > 0){ tris.push_back(t[0]); tris.push_back(t[1]); tris.push_back(t[2]); }
		else      { tris.push_back(t[0]); tris.push_back(t[2]); tris.push_back(t[1]); }
	}
}

/*!
 * 最大の連結成分の三角形のみ残し，使われていない頂点を除いて番号を詰める
 */
inline void compactMesh(std::vector<glm::vec2>& pts, std::vector<unsigned int>& tris)
{
	int nv = (int)pts.size(), nt = (int)tris.size()/3;

	// 頂点を共有する三角形の連結成分(union-find)
	std::vector<int> parent(nv);
	for(int i = 0; i < nv; ++i) parent[i] = i;
	std::function<int(int)> find = [&](int x){ return parent[x] == x ? x : (parent[x] = find(parent[x])); };
	for(int t = 0; t < nt; ++t){
		int a = find(tris[3*t]), b = find(tris[3*t+1]), c = find(tris[3*t+2]);
		parent[b] = a; parent[find(c)] = find(a);
	}
	std::vector<int> cnt(nv, 0);
	for(int t = 0; t < nt; ++t) cnt[find(tris[3*t])]++;
	int root = (int)(std::max_element(cnt.begin(), cnt.end())-cnt.begin());

	std::vector<int> newidx(nv, -1);
	std::vector<glm::vec2> npts;
	std::vector<unsigned int> ntris;
	for(int t = 0; t < nt; ++t){
		if(find(tris[3*t]) != root) continue;
		for(int o = 0; o < 3; ++o){
			int v = tris[3*t+o];
			if(newidx[v] < 0){ newidx[v] = (int)npts.size(); npts.push_back(pts[v]); }
			ntris.push_back(newidx[v]);
		}
	}
	pts.swap(npts);
	tris.swap(ntris);
}

} // namespace rx_imgmesh


//-----------------------------------------------------------------------------
// メッシュ生成
//-----------------------------------------------------------------------------
/*!
 * 画像からオブジェクト形状に沿ったメッシュを生成
 */
inline bool GenerateImageMesh(const unsigned char* img, int w, int h, int c, int wstep,
					   glm::vec2 c1, glm::vec2 c2, const rxImageMeshParams& prm,
					   std::vector<glm::vec2>& points, std::vector<glm::vec2>& texcoords, std::vector<unsigned int>& tris)
{
	using namespace rx_imgmesh;

	if(!img || w <= 0 || h <= 0) return false;
	double sp = prm.spacing;

	// 画像 -> スカラー場
	rxImageField fld;
	if(!makeField(img, w, h, c, wstep, c1, c2, prm, fld)) return false;

	// 輪郭抽出と等間隔サンプリング -> 境界頂点
	std::vector< std::vector<glm::vec2> > loops;
	extractContours(fld, loops);

	std::vector<glm::vec2> pts;
	std::unordered_set<long long> bedges;	// 輪郭上で隣り合う頂点の組(番号は後で頂点数を掛けて符号化)
	std::vector< std::pair<int, int> > bpairs;
	for(auto& loop : loops){
		std::vector<glm::vec2> rs;
		resampleLoop(loop, sp, rs);
		int base = (int)pts.size(), m = (int)rs.size();
		for(int k = 0; k < m; ++k){
			pts.push_back(rs[k]);
			bpairs.push_back(std::make_pair(base+k, base+(k+1)%m));
		}
	}
	int nb = (int)pts.size();
	if(nb < 3) return false;

	// 境界頂点を種にしてPoisson Disk Samplingで内部頂点を追加
	glm::vec2 minp = pts[0], maxp = pts[0];
	for(auto& p : pts){ minp = glm::min(minp, p); maxp = glm::max(maxp, p); }
	minp -= glm::vec2(float(sp)); maxp += glm::vec2(float(sp));
	rxUniformPoissonDiskSampler sampler(minp, maxp, 1000000, 30, sp);
	sampler.Generate(pts);

	// 背景にある点を除く
	std::vector<glm::vec2> tmp(pts.begin(), pts.begin()+nb);
	for(int i = nb; i < (int)pts.size(); ++i){
		if(fld.inside(pts[i])) tmp.push_back(pts[i]);
	}
	pts.swap(tmp);

	long long np = (long long)pts.size();
	for(auto& e : bpairs) bedges.insert(std::min(e.first, e.second)*np+std::max(e.first, e.second));

	// 三角形分割
	triangulate(pts, bedges, fld, tris);

	// 内部頂点の平滑化(隣接頂点の平均へ移動)と再三角形分割で三角形の形状を改善
	for(int it = 0; it < prm.smooth_iter; ++it){
		std::vector<glm::vec2> sum(pts.size(), glm::vec2(0.0f));
		std::vector<int> num(pts.size(), 0);
		for(size_t t = 0; t < tris.size(); t += 3){
			for(int o = 0; o < 3; ++o){
				int a = tris[t+o], b = tris[t+(o+1)%3];
				sum[a] += pts[b]; num[a]++;
				sum[b] += pts[a]; num[b]++;
			}
		}
		for(int i = nb; i < (int)pts.size(); ++i){
			if(num[i] == 0) continue;
			glm::vec2 p = sum[i]/float(num[i]);
			if(fld.inside(p)) pts[i] = p;
		}
		triangulate(pts, bedges, fld, tris);
	}

	// 最大連結成分のみ残して頂点番号を詰める
	compactMesh(pts, tris);
	if(tris.empty()) return false;

	// テクスチャ座標(画像全体が[0,1]x[0,1])
	points = pts;
	texcoords.resize(points.size());
	for(size_t i = 0; i < points.size(); ++i){
		texcoords[i] = glm::vec2(float((points[i][0]-fld.org[0])/(fld.s*w)), float((points[i][1]-fld.org[1])/(fld.s*h)));
	}

	return true;
}


/*!
 * 三角形メッシュの品質(内角)の計算
 */
inline rxMeshQuality CalMeshQuality(const std::vector<glm::vec2>& points, const std::vector<unsigned int>& tris)
{
	rxMeshQuality q;
	q.min_angle = 180.0; q.avg_min_angle = 0.0; q.num_bad = 0;
	int nt = (int)tris.size()/3;
	if(nt == 0){ q.min_angle = 0.0; return q; }
	const double RAD2DEG = 180.0/3.14159265358979323846;
	for(int t = 0; t < nt; ++t){
		double amin = 180.0;
		for(int o = 0; o < 3; ++o){
			glm::vec2 a = points[tris[3*t+(o+1)%3]]-points[tris[3*t+o]];
			glm::vec2 b = points[tris[3*t+(o+2)%3]]-points[tris[3*t+o]];
			double ang = atan2(fabs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])*RAD2DEG;
			amin = std::min(amin, ang);
		}
		q.min_angle = std::min(q.min_angle, amin);
		q.avg_min_angle += amin;
		if(amin < 20.0) q.num_bad++;
	}
	q.avg_min_angle /= nt;
	return q;
}



#endif // #ifdef _RX_MESH_IMAGE_H_
