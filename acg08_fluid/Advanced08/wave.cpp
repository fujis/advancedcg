/*! 
  @file rx_wave.cpp
	
  @brief SWEを使った波のシミュレーション
 
  @author Makoto Fujisawa
  @date 2018-12
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "wave.h"


//-----------------------------------------------------------------------------
// rxWaveクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 * @param[in] n グリッド数
 * @param[in] wind_speed 風の速度
 * @param[in] wind_dir 風の方向
 */
rxWave::rxWave(int n, double scale, double depth)
{
	// ハイトフィールドの解像度と全体の大きさ
	m_iNx = n;
	m_iNy = n;
	m_fScale = scale;

	// 水深
	m_fWaterDepth = depth;

	// 重力
	m_fGravity = 9.81;

	// 配列の初期化
	m_vH.resize(m_iNx*m_iNy, 0.0);
	m_vHprev.resize(m_iNx*m_iNy, 0.0);
	m_vU.resize(m_iNx*m_iNy, 0.0);
	m_vV.resize(m_iNx*m_iNy, 0.0);
	m_vUprev.resize(m_iNx*m_iNy, 0.0);
	m_vVprev.resize(m_iNx*m_iNy, 0.0);

	Init();
}

/*!
 * デストラクタ
 */
rxWave::~rxWave()
{
}

/*!
 * 波の初期化
 */
void rxWave::Init(void)
{
	// ハイトフィールドの初期化
	InitWave();

	// メッシュ作成
	Vec3 c1(-m_fScale/2, 0.0, -m_fScale/2);
	Vec3 c2(m_fScale/2, 0.0, m_fScale/2);
	generateMesh(c1, c2);

	// ハイトフィールドメッシュの更新
	updateMesh(m_vH);

	// 法線再計算
	CalVertexNormals(m_Mesh);

	m_iStep = 0;
}


/*!
 * ハイトフィールドの更新
 * @param[in] dt 時間ステップ幅
 */
int rxWave::Update(double dt)
{
	// 前ステップの値のコピー
	std::copy(m_vH.begin(), m_vH.end(), m_vHprev.begin());
	std::copy(m_vU.begin(), m_vU.end(), m_vUprev.begin());
	std::copy(m_vV.begin(), m_vV.end(), m_vVprev.begin());

	// SWEによるハイトフィールドの更新
	updateHeightBySWE(&m_vH[0], &m_vHprev[0], &m_vU[0], &m_vV[0], &m_vUprev[0], &m_vVprev[0], dt);

	// ハイトフィールドメッシュの更新
	updateMesh(m_vH);
	//updateMesh(m_vHtmp);

	// 法線再計算
	CalVertexNormals(m_Mesh);

	m_iStep++;
	return m_iStep;
}

/*!
 * OpenGLによるハイトフィールドメッシュの描画
 * @param[in] drw 描画フラグ
 */
void rxWave::Draw(int drw)
{
	m_Mesh.Draw(drw, 0.02, false);
}

/*!
* OpenGLによるハイトフィールドメッシュの描画
*  - 頂点ピック用描画
* @param[in] drw 描画フラグ
*/
void rxWave::DrawForPick(void)
{
	double dx = 0.5*m_fDx, dz = 0.5*m_fDy;
	for(int k = 0; k < m_iNy; ++k){
		for(int i = 0; i < m_iNx; ++i){
			int idx = IDX(i, k);
			glLoadName(idx);
			glBegin(GL_QUADS);
			glVertex3dv(m_Mesh.vertices[idx]+Vec3(-dx, 0, -dz));
			glVertex3dv(m_Mesh.vertices[idx]+Vec3( dx, 0, -dz));
			glVertex3dv(m_Mesh.vertices[idx]+Vec3( dx, 0,  dz));
			glVertex3dv(m_Mesh.vertices[idx]+Vec3(-dx, 0,  dz));
			glEnd();
		}
	}
}


/*!
 * 初期ハイトフィールドh0の計算
 */
void rxWave::InitWave(void)
{
	double offset = m_iNx/4;
	for(int j = 0; j < m_iNy; ++j){
		for(int i = 0; i < m_iNx; ++i){
			int idx = IDX(i, j);
			m_vH[idx] = m_fWaterDepth;
			m_vHprev[idx] = m_fWaterDepth;
			m_vU[idx] = 0.0;
			m_vV[idx] = 0.0;
			m_vUprev[idx] = 0.0;
			m_vVprev[idx] = 0.0;

			if(i >= m_iNx/2+offset-3 && i <= m_iNx/2+offset+2 && j >= m_iNy/2+offset-3 && j <= m_iNy/2+offset+2){
				m_vH[idx] = 0.3+m_fWaterDepth;
			}
		}
	}
}

/*!
 * 高さ値の直接変更
 * @param[in] idx 頂点インデックス
 * @param[in] h 追加する高さ
 */
void rxWave::AddHeight(int idx, double h)
{
	if(idx < 0 || idx >= m_iNx*m_iNy) return;
	m_vH[idx] += h;
	m_Mesh.vertices[idx][1] = m_vH[idx]-m_fWaterDepth;
}
void rxWave::AddHeight(int i, int j, double h){ AddHeight(IDX(i, j), h); }


/*!
 * SWEによるハイトフィールドの更新
 * @param[in] dt 時間ステップ幅
 */
void rxWave::updateHeightBySWE(double *h_new, double *h, double *u_new, double *v_new, double *u, double *v, double dt)
{
	double dx = m_fDx;
	double dy = m_fDy;
	double g = m_fGravity;

	// 移流項をバックトレースで解く
	//  - [Layton2002]のFig.2前後の説明にあるsemi-Lagrangian法
	for(int j = 1; j < m_iNy-1; ++j){
		for(int i = 1; i < m_iNx-1; ++i){
			int idx = IDX(i, j);
			double x = i*dx, y = j*dy; // グリッド位置
			
			// バックトレースした位置
			x -= dt*u[idx];
			y -= dt*v[idx];

			// バックトレースした位置の周囲のグリッド情報
			int i0 = (int)(x/dx), j0 = (int)(y/dy);
			i0 = RX_CLAMP(i0, 0, m_iNx-2);
			j0 = RX_CLAMP(j0, 0, m_iNy-2);
			double s = (x-i0*dx)/dx, t = (y-j0*dy)/dy;

			// バックトレースした位置でのu,v,hを線形補間で求める
			u_new[idx] = u[IDX(i0, j0)]*(1-s)*(1-t) + u[IDX(i0+1, j0)]*s*(1-t) + u[IDX(i0, j0+1)]*(1-s)*t + u[IDX(i0+1, j0+1)]*s*t;
			v_new[idx] = v[IDX(i0, j0)]*(1-s)*(1-t) + v[IDX(i0+1, j0)]*s*(1-t) + v[IDX(i0, j0+1)]*(1-s)*t + v[IDX(i0+1, j0+1)]*s*t;
			h_new[idx] = h[IDX(i0, j0)]*(1-s)*(1-t) + h[IDX(i0+1, j0)]*s*(1-t) + h[IDX(i0, j0+1)]*(1-s)*t + h[IDX(i0+1, j0+1)]*s*t;
		}
	}

	// その他の項を解く
	//  - [Layton2002]の式(17)(18)で地形の高さbが変化しない(db/dx=db/dy=0)として計算している
	double d = m_fWaterDepth;
	for(int j = 1; j < m_iNy-1; ++j){
		for(int i = 1; i < m_iNx-1; ++i){
			int idx = IDX(i, j);
			// 式(17)によるu,vの更新
			u_new[idx] -= 0.5*dt*(g*(h_new[IDX(i+1, j)]-h_new[IDX(i-1, j)])/dx);
			v_new[idx] -= 0.5*dt*(g*(h_new[IDX(i, j+1)]-h_new[IDX(i, j-1)])/dy);

			// 式(18)によるhの更新(本当は式(25)のように陰解法で解いた方がよいがここでは簡易的に陽的に解いている)
			h_new[idx] -= 0.5*dt*(d+h_new[idx])*((u_new[IDX(i+1, j)]-u_new[IDX(i-1, j)])/dx+(v_new[IDX(i, j+1)]-v_new[IDX(i, j-1)])/dy);
		}
	}


	// 境界処理
	for(int i = 0; i < m_iNx; ++i){
		h_new[IDX(i, 0)] = h_new[IDX(i, 1)]; h_new[IDX(i, m_iNy-1)] = h_new[IDX(i, m_iNy-2)];
		u_new[IDX(i, 0)] = u_new[IDX(i, 1)]; u_new[IDX(i, m_iNy-1)] = u_new[IDX(i, m_iNy-2)];
		v_new[IDX(i, 0)] = v_new[IDX(i, 1)]; v_new[IDX(i, m_iNy-1)] = v_new[IDX(i, m_iNy-2)];
	}
	for(int j = 0; j < m_iNy; ++j){
		h_new[IDX(0, j)] = h_new[IDX(1, j)]; h_new[IDX(m_iNx-1, j)] = h_new[IDX(m_iNx-2, j)];
		u_new[IDX(0, j)] = u_new[IDX(1, j)]; u_new[IDX(m_iNx-1, j)] = u_new[IDX(m_iNx-2, j)];
		v_new[IDX(0, j)] = v_new[IDX(1, j)]; v_new[IDX(m_iNx-1, j)] = v_new[IDX(m_iNx-2, j)];
	}
}

/*!
* ガウスフィルタの係数行列の計算
 * @param[out] w 係数行列(size*size)
 * @param[in] size フィルタサイズ(3,5,7,...)
 * @param[in] sigma 分散
 */
inline void calGaussianWeight(vector<double> &w, int size, double sigma)
{
	const double PI = 3.14159265359;
	w.resize(size*size);
	double sum = 0.0;
	int c = size/2; // 中心画素
	double x, y;
	double s2 = sigma*sigma;
	double a = 0.5/(PI*s2);
	for(int j = 0; j < size; ++j){
		y = j-c;
		for(int i = 0; i < size; ++i){
			x = i-c;

			w[i+size*j] = a*exp(-0.5*(x*x+y*y)/s2);
			sum += w[i+size*j];
		}
	}

	// 正規化
	for(int i = 0; i < size*size; ++i){
		w[i] /= sum;
	}
}
void rxWave::gaussianHeight(const vector<double> &h, vector<double> &hs, int size, int sigma)
{
	// ガウスフィルタの準備
	vector<double> w(size*size, 0.0);
	calGaussianWeight(w, size, sigma);
	int n = (int)(size/2);

	hs.resize(m_iNx*m_iNy);
	for(int i = 0; i < m_iNx*m_iNy; ++i) hs[i] = h[i];

	// ガウスフィルタを掛けたハイトフィールドの作成
	for(int k = 0; k < 2; ++k){
		bool non_zero = true;
		for(int j = 0; j < m_iNy; ++j){
			for(int i = 0; i < m_iNx; ++i){
				int idx = IDX(i, j);
				double hws = 0.0;
				for(int t = 0; t < size; ++t){
					int ja = j+(t-n);
					if(ja < 0) ja = 0;
					if(ja >= m_iNy) ja = m_iNy-1;
					for(int s = 0; s < size; ++s){
						int ia = i+(s-n);
						if(ia < 0) ia = 0;
						if(ia >= m_iNx) ia = m_iNx-1;
						//if(ia >= 0 && ia < m_iNx && ja >= 0 && ja < m_iNy){
						hws += w[s+size*t]*hs[ia+ja*m_iNx];
						//}
					}
				}
				hs[idx] = hws;
				if(fabs(hs[idx]) < 1.0e-8) non_zero = false;
			}
		}
		if(non_zero) break;
	}

}


/*!
 * ハイトフィールドデータの入力
 */
void rxWave::InputHeightField(const string &fn)
{
	inputHeightField(fn, m_vH, m_iNx, m_iNy);

	// メッシュ作成
	Vec3 c1(-m_fScale/2, 0.0, -m_fScale/2);
	Vec3 c2(m_fScale/2, 0.0, m_fScale/2);
	generateMesh(c1, c2);

	// ハイトフィールドメッシュの更新
	updateMesh(m_vH);

	// 法線再計算
	CalVertexNormals(m_Mesh);

}
void rxWave::inputHeightField(const string &fn, vector<double> &hf, int &nx, int &ny)
{
		// ファイル入力
	ifstream fs;
	fs.open(fn.c_str(), ios::in);
	if(fs){
		fs >> nx >> ny;

		if(hf.size() != nx*ny) hf.resize(nx*ny);

		for(int j = 0; j < ny; ++j){
			for(int i = 0; i < nx; ++i){
				int idx = IDX(i, j);
				fs >> hf[idx];
				hf[idx] += m_fWaterDepth;
			}
		}
		fs.close();

		cout << "Height field have been read from " << fn << endl;
	}
	else{
		cout << "There is no file : " << fn << endl;
	}
}


/*!
* ハイトフィールドデータの出力
*/
void rxWave::OutputHeightField(const string &fn)
{
	outputHeightField(fn, m_vH, m_iNx, m_iNy);
}
/*!
* ハイトフィールドデータの出力
*/
void rxWave::outputHeightField(const string &fn, const vector<double> &hf, int nx, int ny)
{
	// ファイル出力
	ofstream fs;
	fs.open(fn.c_str(), ios::out);
	fs << nx << " " << ny << endl;
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i < nx; ++i){
			int idx = IDX(i, j);
			fs << hf[idx]-m_fWaterDepth << endl;
		}
	}
	fs.close();
	cout << "Height field have been written to " << fn << endl;
}



/*!
 * n×nの頂点を持つメッシュ生成(x-z平面)
 * @param[in] c1,c2 2端点座標
 */
void rxWave::generateMesh(Vec3 c1, Vec3 c2)
{
	if(!m_Mesh.vertices.empty()){
		m_Mesh.vertices.clear();
		m_Mesh.faces.clear();
	}

	// 頂点座標生成
	double dx = (c2[0]-c1[0])/m_iNx;
	double dz = (c2[2]-c1[2])/m_iNy;
	m_Mesh.vertices.resize(m_iNx*m_iNy);
	for(int k = 0; k < m_iNy; ++k){
		for(int i = 0; i < m_iNx; ++i){
			Vec3 pos;
			pos[0] = c1[0]+i*dx;
			pos[1] = c1[1];
			pos[2] = c1[2]+k*dz;
			m_Mesh.vertices[IDX(i, k)] = pos;
		}
	}

	m_fDx = dx;
	m_fDy = dz;

	// メッシュ作成
	for(int k = 0; k < m_iNy-1; ++k){
		for(int i = 0; i < m_iNx-1; ++i){
			rxFace face;
			face.resize(3);

			face[0] = IDX(i, k);
			face[1] = IDX(i+1, k+1);
			face[2] = IDX(i+1, k);
			m_Mesh.faces.push_back(face);

			face[0] = IDX(i, k);
			face[1] = IDX(i, k+1);
			face[2] = IDX(i+1, k+1);
			m_Mesh.faces.push_back(face);
		}
	}
}


/*!
 * ハイトフィールドに従ってメッシュ頂点のy座標値を更新
 * @param[in] h ハイトフィールド(m_iNx*m_iNy)
 */
void rxWave::updateMesh(const vector<double> &h)
{
	for(int k = 0; k < m_iNy; ++k){
		for(int i = 0; i < m_iNx; ++i){
			int idx = IDX(i, k);
			double y = h[idx];
			m_Mesh.vertices[idx][1] = y-m_fWaterDepth;
		}
	}
}
