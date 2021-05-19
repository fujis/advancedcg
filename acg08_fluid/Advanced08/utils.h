/*! 
  @file rx_utils.h

  @brief 各種便利な関数や定数など
 
  @author Makoto Fujisawa
  @date 2020-07
*/
// FILE --rx_utils.h--

#ifndef _RX_UTILS_H_
#define _RX_UTILS_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// C標準
#include <ctime>
#include <cmath>
//#include <cctype>
#include <cstdio>
#include <cassert>

#include <direct.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

// STL
#include <vector>
#include <string>
#include <map>
#include <bitset>
#include <algorithm>

// OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>


// ユーティリティ
#include "rx_utility.h"

#include "rx_material.h"
#include "rx_matrix.h"

#include "rx_timer.h"

using namespace std;


//-----------------------------------------------------------------------------
// MARK:定義
//-----------------------------------------------------------------------------
typedef unsigned char byte;

#define STR(x) #x

#define RXFRMT(a) (boost::format("%.3f") % (a))
#define RXFRMTE(a) (boost::format("%.3e") % (a))

// 計算結果出力のデフォルトフォルダ
const string RX_DEFAULT_RESULT_DIR = "result/";
const string RX_DEFAULT_IMAGE_DIR  = RX_DEFAULT_RESULT_DIR+"images/";
const string RX_DEFAULT_MESH_DIR   = RX_DEFAULT_RESULT_DIR+"mesh/";
const string RX_DEFAULT_DATA_DIR   = RX_DEFAULT_RESULT_DIR+"data/";



// グローバルな時間計測用(RX_USE_TIMERの定義の有無で時間計測の有無を変更できる)
static rxTimerAvg g_Time;

#ifdef RX_USE_TIMER
#define RXTIMER_CLEAR g_Time.ClearTime()
#define RXTIMER_RESET g_Time.ResetTime()
#define RXTIMER(x) g_Time.Split(x)
#define RXTIMER_STOP g_Time.StopTime()
#define RXTIMER_PRINT g_Time.Print()
#define RXTIMER_STRING(x) g_Time.PrintToString(x)
#define RXTIMER_FPS(x) g_Time.GetFPSString(x)
#else
#define RXTIMER_CLEAR
#define RXTIMER_RESET
#define RXTIMER_STOP
#define RXTIMER(x) 
#define RXTIMER_PRINT
#define RXTIMER_STRING(x)
#define RXTIMER_FPS
#endif


//-----------------------------------------------------------------------------
// 関数の定義
//-----------------------------------------------------------------------------
//! グラデーション色の生成
template<class T> 
inline void RX_COLOR_RAMP(T t, T *r)
{
	const int ncolors = 7;
	T c[ncolors][3] = {	{ 0.0, 0.0, 1.0 },
						{ 0.0, 0.5, 1.0 },
						{ 0.0, 1.0, 1.0 },
						{ 0.0, 1.0, 0.0 },
						{ 1.0, 1.0, 0.0 },
						{ 1.0, 0.0, 0.0 },
						{ 1.0, 0.0, 1.0 } };
	//T c[ncolors][3] = { { 1.0, 0.0, 0.0 },
	//					{ 1.0, 0.5, 0.0 },
	//					{ 1.0, 1.0, 0.0 },
	//					{ 0.0, 1.0, 0.0 },
	//					{ 0.0, 1.0, 1.0 },
	//					{ 0.0, 0.0, 1.0 },
	//					{ 1.0, 0.0, 1.0 } };
	t = t*(ncolors-1);
	int i = (int)t;
	T u = t-floor(t);
	r[0] = RX_LERP(c[i][0], c[i+1][0], u);
	r[1] = RX_LERP(c[i][1], c[i+1][1], u);
	r[2] = RX_LERP(c[i][2], c[i+1][2], u);
}



inline void EulerToMatrix(double *m, double pitch, double yaw, double roll)
{
	yaw   = RX_TO_RADIANS(yaw);
	pitch = RX_TO_RADIANS(pitch);
	roll  = RX_TO_RADIANS(roll);
 
	double cy = cos(yaw); 
	double sy = sin(yaw); 
	double cp = cos(pitch); 
	double sp = sin(pitch); 
	double cr = cos(roll);
	double sr = sin(roll);
 
	double cc = cy*cr; 
	double cs = cy*sr; 
	double sc = sy*cr; 
	double ss = sy*sr;
 
	m[0]  = cc+sp*ss;
	m[1]  = cs-sp*sc;
	m[2]  = -sy*cp;
	m[3]  = 0.0;
	
	m[4]  = -cp*sr;
	m[5]  = cp*cr; 
	m[6]  = -sp;
	m[7]  = 0.0;
	
	m[8]  = sc-sp*cs;
	m[9]  = ss+sp*cc;
	m[10] = cy*cp;
	m[11] = 0.0;
 
	m[12] = 0.0;
	m[13] = 0.0;
	m[14] = 0.0;
	m[15] = 1.0;
}

/*!
 * 3x3行列からOpenGLの変換行列を取得
 * @param[in] 
 * @param[out] 
 * @return 
 */
inline void GetGLMatrix(double mat[9], GLfloat glmat[16])
{
	glmat[0]  = mat[0];
	glmat[1]  = mat[3]; 
	glmat[2]  = mat[6]; 
	glmat[3]  = 0.0f; 
	glmat[4]  = mat[1]; 
	glmat[5]  = mat[4]; 
	glmat[6]  = mat[7]; 
	glmat[7]  = 0.0f; 
	glmat[8]  = mat[2]; 
	glmat[9]  = mat[5]; 
	glmat[10] = mat[8]; 
	glmat[11] = 0.0f; 
	glmat[12] = 0.0f; 
	glmat[13] = 0.0f; 
	glmat[14] = 0.0f; 
	glmat[15] = 1.0f;
}


/*!
 * 青->緑->赤->白と変化するグラデーション色生成
 * @param[out] col 生成された色
 * @param[in] x 値
 * @param[in] xmin 最小値
 * @param[in] xmax 最大値
 */
inline void Gradation(double col[3], double x, const double xmin = 0.0, const double xmax = 1.0)
{
	double l = xmax-xmin;
	if(fabs(l) < 1e-10) return;
	
	const int ncolors = 7;
	double base[ncolors][3] = { {0.0, 0.0, 0.0},
								{0.0, 0.0, 1.0},
								{0.0, 1.0, 1.0},
								{0.0, 1.0, 0.0},
								{1.0, 1.0, 0.0},
								{1.0, 0.0, 0.0},
								{1.0, 1.0, 1.0} };
	x = RX_CLAMP<double>(((x-xmin)/l), 0.0, 1.0)*(ncolors-1);
	int i = (int)x;
	double dx = x-floor(x);
	col[0] = RX_LERP(base[i][0], base[i+1][0], dx);
	col[1] = RX_LERP(base[i][1], base[i+1][1], dx);
	col[2] = RX_LERP(base[i][2], base[i+1][2], dx);
}

static inline void CalExtendedBBox(Vec3 &minp, Vec3 &maxp, double expansion_rate = 0.01)
{
	// 元のBBoxの各辺の長さと原点(最小座標)
	Vec3 size   = maxp-minp;
	Vec3 origin = minp;
	
	// 辺の長さを拡大し，それに合わせて原点位置を修正(元のBBoxの中心位置を保つようにする)
	size   *= (1.0+expansion_rate);
	origin -= 0.5*expansion_rate*size;

	// 拡張した新しいBBox
	minp = origin;
	maxp = origin+size;

}


inline int IDX4(int row, int col){ return (row | (col<<2)); }

/*!
 * 4x4行列とベクトルの掛け算
 *  | m[0]  m[1]  m[2]  m[3]  |
 *  | m[4]  m[5]  m[6]  m[7]  |
 *  | m[8]  m[9]  m[10] m[11] |
 *  | m[12] m[13] m[14] m[15] |
 * @param[in] m 元の行列
 * @param[in] v ベクトル
 * @return 積の結果のベクトル(m*v)
 */
inline Vec3 MulMat4Vec3(GLfloat* m, const Vec3& v)
{
	Vec3 d(0.0);
	d[0] = (v[0]*m[IDX4(0, 0)]+v[1]*m[IDX4(0, 1)]+v[2]*m[IDX4(0, 2)]);
	d[1] = (v[0]*m[IDX4(1, 0)]+v[1]*m[IDX4(1, 1)]+v[2]*m[IDX4(1, 2)]);
	d[2] = (v[0]*m[IDX4(2, 0)]+v[1]*m[IDX4(2, 1)]+v[2]*m[IDX4(2, 2)]);
	return d;
}


/*!
 * 任意ベクトル周りの回転
 * @param[in] pos 元の座標値
 * @param[in] axis 回転軸
 * @param[in] ang 回転角度(rad)
 * @return 回転した座標値
 */
inline Vec3 Rotate(const Vec3& pos, const Vec3& axis, const double& ang)
{
	Vec3 pos1;    // 回転後の座標値

	double c = cos(ang);
	double s = sin(ang);
	double x, y, z;
	x = axis[0]; y = axis[1]; z = axis[2];

	// | xx(1-c)+c	xy(1-c)-zs	xz(1-c)+ys	0 |
	// | yx(1-c)-zs	yy(1-c)+c	yz(1-c)-xs	0 |
	// | xz(1-c)-ys	yz(1-c)+xs	zz(1-c)+c	0 |
	// | 0			0			0			1 |
	pos1[0] = (x*x*(1.0-c)+c)*pos[0] +(x*y*(1.0-c)-z*s)*pos[1] +(x*z*(1.0-c)+y*s)*pos[2];
	pos1[1] = (y*x*(1.0-c)-z*s)*pos[0]   +(y*y*(1.0-c)+c)*pos[1] +(y*z*(1.0-c)-x*s)*pos[2];
	pos1[2] = (x*z*(1.0-c)-y*s)*pos[0] +(y*z*(1.0-c)+x*s)*pos[1]   +(z*z*(1.0-c)+c)*pos[2];

	return pos1;
}



/*!
 * "\n"が含まれるstringを複数のstringに分解する
 * @param[in] org 元の文字列
 * @param[in] div 分解後の文字列配列
 */
static inline void divideString(const string& org, vector<string>& div)
{
	size_t pos0 = 0, pos1 = 0;
	while(pos1 != string::npos){
		pos1 = org.find("\n", pos0);

		div.push_back(org.substr(pos0, pos1-pos0));

		pos0 = pos1+1;
	}
}


//-----------------------------------------------------------------------------
// カスタムストリーム出力
//-----------------------------------------------------------------------------
class rxLog
{
	fstream m_ofLog;
public:

	rxLog(const char *filename)
	{
		m_ofLog.open(filename, ios::out);
		if(!m_ofLog || !m_ofLog.is_open() || m_ofLog.bad() || m_ofLog.fail()){
			return;
		}
	}
	~rxLog()
	{
		if(m_ofLog && m_ofLog.is_open()) m_ofLog.close();
	}

	//! <<オペレータを設定
	template<typename T>
	rxLog& operator<<(const T &a)
	{
		std::cout << a;
		if(m_ofLog) m_ofLog << RX_TO_STRING(a);
		return *this;
	}

	// std::coutの型
	typedef std::basic_ostream<char, std::char_traits<char> > TypeCout;

	// std::endlのためのオペレータ<<を定義
	// (std::endlはstd::coutを引数としてとる関数)
	rxLog& operator<<(TypeCout& (*manip)(TypeCout&))
	{
		manip(std::cout);
		if(m_ofLog) m_ofLog << std::endl;
		return *this;
	}
};

static rxLog RXCOUT("log/sph.log");
static rxLog RXDOUT("log/sph.dat");

//#define RXPRINT(x) RXCOUT << boost::format(x)





//-----------------------------------------------------------------------------
// コンテナ関連
//-----------------------------------------------------------------------------

/*!
 * vectorのidx番目の要素を削除(0スタート)
 * @param[in] src vectorコンテナ
 * @param[in] idx 削除要素インデックス
 * @return 削除の可否
 */
template<class T> 
inline bool EraseSTLVectori(vector<T> &src, int idx)
{
	vector<T>::iterator iter = src.begin();
	int i = 0;
	while(iter != src.end()){
		if(i == idx){
			src.erase(iter++);
			return true;
		}
		else{
			++iter;
		}

		i++;
	}

	return false;
}

/*!
 * vectorの特定の要素を削除
 * @param[in] src vectorコンテナ
 * @param[in] comp_func 削除条件関数
 * @return 削除の可否
 */
template<class T> 
inline bool EraseSTLVector(vector<T> &src, bool (*comp_func)(T))
{
	vector<T>::iterator iter = src.begin();
	while(iter != src.end()){
		if(comp_func(*iter)){
			src.erase(iter++);
			return true;
		}
		else{
			++iter;
		}
	}

	return false;
}

/*!
 * mapの特定の要素を削除
 * @param[in] src マップ
 * @param[in] comp_func 削除条件関数
 * @return 削除の可否
 */
template<class T1, class T2> 
inline bool EraseSTLMap(map<T1, T2> &src, bool (*comp_func)(pair<T1, T2>))
{
	map<T1, T2>::iterator iter = src.begin();
	while(iter != src.end()){
		if(comp_func(*iter)){
			src.erase(iter++);
		}
		else{
			++iter;
		}
	}

	return (iter != src.end());
}



/*!
 * "(x, y, z)"の形式の文字列からVec3型へ変換
 *  - (x)となっていたら(x, x, x)とする．
 * @param[in] s 文字列
 * @param[out] v 値
 * @return 要素記述数
 */
inline int StringToVec3(const string &s, Vec3 &v)
{
	int vcount = 0;
	size_t pos;
	v = Vec3(0.0);
	if((pos = s.find('(')) != string::npos){
		while(pos != string::npos && vcount < 3){
			size_t pos1 = pos;
			if((pos1 = s.find(',', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				pos = pos1;
			}
			else if((pos1 = s.find(')', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				break;
			}
			else{
				break;
			}
		}
	}
	if(vcount < 3){
		for(int i = vcount; i < 3; ++i){
			v[i] = v[vcount-1];
		}
	}

	return vcount;
}

/*!
* "(x, y, z, w)"の形式の文字列からVec4型へ変換
*  - (x)となっていたら(x, x, x, x)とする．
* @param[in] s 文字列
* @param[out] v 値
* @return 要素記述数
*/
inline int StringToVec4(const string &s, Vec4 &v)
{
	int vcount = 0;
	size_t pos;
	v = Vec4(0.0);
	if((pos = s.find('(')) != string::npos){
		while(pos != string::npos && vcount < 4){
			size_t pos1 = pos;
			if((pos1 = s.find(',', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				pos = pos1;
			} else if((pos1 = s.find(')', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				break;
			} else{
				break;
			}
		}
	}
	if(vcount < 4){
		for(int i = vcount; i < 4; ++i){
			v[i] = v[vcount-1];
		}
	}

	return vcount;
}


inline Vec3 Vec3f(const float* v)
{
	Vec3 x;
	x.data[0] = v[0];
	x.data[1] = v[1];
	x.data[2] = v[2];
	return x;
}


//-----------------------------------------------------------------------------
// サーモグラフ
//-----------------------------------------------------------------------------
/*!
 * col0->col1と変化するサーモグラフ用の色生成
 * @param 
 * @return 
 */
template<class Type> 
inline void Gradation(Type x, Type minx, Type maxx, 
					  Type col0[3], Type col1[3], 
					  Type &r, Type &g, Type &b)
{
	Type l = maxx-minx;
	
	if(fabs(l) < 1e-8){
		r = col0[0]; g = col0[1]; b = col0[2];
		return;
	}

	if(x <= minx){
		r = col0[0]; g = col0[1]; b = col0[2];
	}
	else if(x >= maxx){
		r = col1[0]; g = col1[1]; b = col1[2];
	}
	else{
		Type col[3];
		for(int i = 0; i < 3; ++i){
			col[i] = col0[i]+(col1[i]-col0[i])*((x-minx)/l);
		}
		r = col[0]; g = col[1]; b = col[2];
	}
}

/*!
 * col0->col1と変化するサーモグラフ用の色生成
 * @param 
 * @return 
 */
template<class Type> 
inline Vec3 Gradation(Type x, Type minx, Type maxx, Type col0[3], Type col1[3])
{
	Type l = maxx-minx;
	Vec3 c;
	
	if(fabs(l) < 1e-8){
		c[0] = col0[0]; c[1] = col0[1]; c[2] = col0[2];
		return c;
	}

	if(x <= minx){
		c[0] = col0[0]; c[1] = col0[1]; c[2] = col0[2];
	}
	else if(x >= maxx){
		c[0] = col1[0]; c[1] = col1[1]; c[2] = col1[2];
	}
	else{
		Type col[3];
		for(int i = 0; i < 3; ++i){
			col[i] = col0[i]+(col1[i]-col0[i])*((x-minx)/l);
		}
		c[0] = col[0]; c[1] = col[1]; c[2] = col[2];
	}

	return c;
}


//-----------------------------------------------------------------------------
// MARK:データファイル出力
//-----------------------------------------------------------------------------
template<class Type> 
static void Output(int w, int h, Type *data, char *name)
{
	// HACK:output
	if(data == NULL) return;

	ofstream fo;
	fo.open(name);

	int size = w*h;
	Type min_ = data[0], max_ = data[0];
	for(int i = 0; i < size; ++i){
		if(data[i] < min_) min_ = data[i];
		if(data[i] > max_) max_ = data[i];
	}
	fo << "[" << min_ << ", " << max_ << "]" << endl;

	for(int j = 0; j < h; ++j){
		for(int i = 0; i < w; ++i){
			int idx = j*w+i;
			fo << data[idx] << " ";
		}
		fo << endl;
	}

	fo.close();
}



template<class Type> 
static void Output_v2(int w, int h, Type *x, Type *y, char *name)
{
	if(x == NULL || y == NULL) return;

	ofstream fo;
	fo.open(name);

	int size = w*h;
	Type minx = x[0], maxx = x[0];
	Type miny = y[0], maxy = y[0];
	for(int i = 0; i < size; ++i){
		if(x[i] < minx) minx = x[i];
		if(x[i] > maxx) maxx = x[i];
		if(y[i] < miny) miny = y[i];
		if(y[i] > maxy) maxy = y[i];
	}
	fo << "[ (" << minx << ", " << miny << "), (" << maxx << ", " << maxy << ") ]" << endl;

	for(int j = 0; j < h; ++j){
		for(int i = 0; i < w; ++i){
			int idx = j*w+i;
			fo << "(" << x[idx] << ", " << y[idx] << ") ";
		}
		fo << endl;
	}

	fo.close();
}

/*!
 * データのファイル出力
 * @param[in] fn ファイル名
 * @param[in] label データラベル
 * @param[in] n データ数
 * @param[in] d データの次元
 * @param[in] rw trueなら既存のファイルに追記
 * @param[in] label データラベル(主に追記の場合に使用)
 */
template<class Type> 
static void Dump(string fn, Type *data, int n, int d, bool rw = false, string label = "")
{
	ofstream fout;
	ios_base::openmode mode = ios::out | (rw ? (ios::app | ios::ate) : 0);
	fout.open(fn.c_str(), mode);
	if(!fout){
		RXCOUT << fn << " couldn't open (DumpData)." << endl;
		return;
	}

	if(!label.empty()) fout << "[" << label << "]" << endl;
	for(int i = 0; i < n; ++i){
		if(d >= 2) fout << "(";
		for(int j = 0; j < d; ++j){
			fout << data[d*i+j] << (j == d-1 ? "" : ", ");
		}
		if(d >= 2) fout << ")";
		fout << endl;
	}
	fout << endl;

	fout.close();
}




//-----------------------------------------------------------------------------
// MARK:ランダムソート
//-----------------------------------------------------------------------------
/*!
 * [0, 1]の範囲の実数乱数の生成
 * @return 実数乱数
 */
inline double Random()
{
	return (double)rand()/(double)RAND_MAX;
}

/*!
 * 指定した範囲の実数乱数の生成
 * @param _max, _min 範囲
 * @return 実数乱数
 */
inline double Random(const double &_max, const double &_min)
{
	return (_max-_min)*Random()+_min;
}






//-----------------------------------------------------------------------------
// MARK:ファイル・フォルダ処理
//-----------------------------------------------------------------------------
/*!
 * ファイル，フォルダの存在確認
 * @param[in] path ファイル・フォルダパス
 */
inline int ExistFile(const string fn)
{
	FILE *fp;
 
	if( (fp = fopen(fn.c_str(), "r")) == NULL ){
		return 0;
	}
 
	fclose(fp);
	return 1;
}

/*!
 * フォルダ区切りの検索
 * @param[in] str ファイル・フォルダパス
 * @param[out] pos 見つかった位置
 */
inline bool FindPathBound(const string &str, string::size_type &pos)
{
	if((pos = str.find_last_of("/")) == string::npos){
		if((pos = str.find_last_of("\\")) == string::npos){
			return false;
		}
	}
	
	return true;
}

/*!
 * ファイル名比較関数(拡張子)
 * @param[in] fn 比較したいファイル名
 * @param[in] ext 拡張子
 * @return fnの拡張子がextと同じならtrue
 */
inline bool SearchCompExt(const string &fn, const string &ext)
{
	return (fn.find(ext, 0) != string::npos);
}


/*!
 * ファイル名生成
 * @param head : 基本ファイル名
 * @param ext  : 拡張子
 * @param n    : 連番
 * @param d    : 連番桁数
 * @return 生成したファイル名
 */
inline string CreateFileName(const string &head, const string &ext, int n, const int &d)
{
	string file_name = head;
	int dn = d-1;
	if(n > 0){
		dn = (int)(log10((double)n))+1;
	}
	else if(n == 0){
		dn = 1;
	}
	else{
		n = 0;
		dn = 1;
	}

	for(int i = 0; i < d-dn; ++i){
		file_name += "0";
	}

	file_name += RX_TO_STRING(n);
	if(!ext.empty() && ext[0] != '.') file_name += ".";
	file_name += ext;

	return file_name;
}


/*!
 * 整数値の下一桁を返す
 * @param[in] x 整数値
 * @return xの下一桁
 */
inline int LowOneDigit(const int &x)
{
	int x1 = (x < 0) ? -x : x;
	double a = 10;

	//INT_MAX = 2147483647
	for(int i = 10; i >= 1; i--){
		a = pow(10.0, (double)i);
		while(x1 > a){
			x1 -= (int)a;
		}
	}

	return x1;
}

/*!
 * 0付きの数字を生成
 * @param[in] n 数字
 * @param[in] d 桁数
 * @return 0付きの数字(string)
 */
inline string GenZeroNo(int n, const int &d)
{
	string zero_no = "";
	int dn = d-1;
	if(n > 0){
		dn = (int)(log10((double)n))+1;
	}
	else if(n == 0){
		dn = 1;
	}
	else{
		n = 0;
		dn = 1;
	}

	for(int i = 0; i < d-dn; ++i){
		zero_no += "0";
	}

	zero_no += RX_TO_STRING(n);

	return zero_no;
}


/*!
 * 秒数を hh:mm:ss の形式に変換
 * @param[in] sec 秒数
 * @param[in] use_msec ミリ秒まで含める(hh:mm:ss.msec)
 * @return hh:mm:ss形式の文字列
 */
inline string GenTimeString(double sec, bool use_msec = false)
{
	long value = (int)(1000*sec+0.5);	// ミリ秒

	unsigned int h = (unsigned int)(value/3600000);	// 時間
	value -= h*3600000;
	unsigned int m = (unsigned int)(value/60000);		// 分
	value -= m*60000;
	unsigned int s = (unsigned int)(value/1000);		// 秒
	value -= s*1000;
	unsigned int ms = (unsigned int)(value);			// ミリ秒

	stringstream ss;
	if(h > 0) ss << GenZeroNo(h, 2) << ":";
	ss << GenZeroNo(m, 2) << ":";
	ss << GenZeroNo(s, 2);
	if(use_msec) ss << "." << GenZeroNo(ms, 3);

	return ss.str();
}

/*!
 * 時刻を hh:mm:ss の形式に変換
 * @param[in] h,m,s 時,分,秒
 * @return hh:mm:ss形式の文字列
 */
inline string GenTimeString(int h, int m, int s)
{
	stringstream ss;
	if(h > 0) ss << GenZeroNo(h, 2) << ":";
	ss << GenZeroNo(m, 2) << ":";
	ss << GenZeroNo(s, 2);
	return ss.str();
}
/*!
 * パスからファイル名のみ取り出す
 * @param[in] path パス
 * @return ファイル名
 */
inline string GetFileName(const string &path)
{
	size_t pos1;
 
	pos1 = path.rfind('\\');
	if(pos1 != string::npos){
		return path.substr(pos1+1, path.size()-pos1-1);
	}
 
	pos1 = path.rfind('/');
	if(pos1 != string::npos){
		return path.substr(pos1+1, path.size()-pos1-1);
	}
 
	return path;
}

/*!
 * パスから拡張子を小文字にして取り出す
 * @param[in] path ファイルパス
 * @return (小文字化した)拡張子
 */
inline string GetExtension(const string &path)
{
	string ext;
	size_t pos1 = path.rfind('.');
	if(pos1 != string::npos){
		ext = path.substr(pos1+1, path.size()-pos1);
		string::iterator itr = ext.begin();
		while(itr != ext.end()){
			*itr = tolower(*itr);
			itr++;
		}
		itr = ext.end()-1;
		while(itr != ext.begin()){	// パスの最後に\0やスペースがあったときの対策
			if(*itr == 0 || *itr == 32){
				ext.erase(itr--);
			}
			else{
				itr--;
			}
		}
	}
 
	return ext;
}

/*!
 * ファイルストリームを開く
 * @param[out] file ファイルストリーム
 * @param[in] path  ファイルパス
 * @param[in] rw    入出力フラグ (1:読込, 2:書込, 4:追記)
 * @return ファイルオープン成功:1, 失敗:0
 */
static inline int OpenFileStream(fstream &file, const string &path, int rw = 1)
{
	file.open(path.c_str(), (rw & 0x01 ? ios::in : 0)|(rw & 0x02 ? ios::out : 0)|(rw & 0x04 ? ios::app : 0));
	if(!file || !file.is_open() || file.bad() || file.fail()){
		return 0;
	}
	return 1;
}

/*!
 * ディレクトリ作成(多階層対応)
 * @param[in] dir 作成ディレクトリパス
 * @return 成功で1,失敗で0 (ディレクトリがすでにある場合も1を返す)
 */
static int MkDir(string dir)
{
	if(_mkdir(dir.c_str()) != 0){
		char cur_dir[512];
		_getcwd(cur_dir, 512);	// カレントフォルダを確保しておく
		if(_chdir(dir.c_str()) == 0){	// chdirでフォルダ存在チェック
			cout << "MkDir : " << dir << " is already exist." << endl;
			_chdir(cur_dir);	// カレントフォルダを元に戻す
			return 1;
		}
		else{
			size_t pos = dir.find_last_of("\\/");
			if(pos != string::npos){	// 多階層の可能性有り
				int parent = MkDir(dir.substr(0, pos+1));	// 親ディレクトリを再帰的に作成
				if(parent){
					if(_mkdir(dir.c_str()) == 0){
						return 1;
					}
					else{
						return 0;
					}
				}
			}
			else{
				return 0;
			}
		}
	}

	return 1;
}



//-----------------------------------------------------------------------------
// 衝突処理/距離計算用関数
//-----------------------------------------------------------------------------
//! AABBの法線方向(内向き)
const Vec3 RX_AABB_NORMALS[6] = { Vec3(1.0,  0.0,  0.0),	// x-
								  Vec3(-1.0,  0.0,  0.0),	// x+
								  Vec3(0.0,  1.0,  0.0),	// y-
								  Vec3(0.0, -1.0,  0.0),	// y+
								  Vec3(0.0,  0.0,  1.0),	// z-
								  Vec3(0.0,  0.0, -1.0) };	// z+

/*!
 * AABBと点の距離
 * @param[in] spos 立方体の中心を原点とした相対座標値
 * @param[in] r    半径(球の場合)
 * @param[in] sgn  立方体の内で距離が正:1,外で正:-1
 * @param[in] vMin 立方体の最小座標値(相対座標)
 * @param[in] vMax 立方体の最大座標値(相対座標)
 * @param[out] d   符号付距離値
 * @param[out] n   最近傍点の法線方向
 */
static inline bool AABB_point_dist(const Vec3& spos, const double& r, const int& sgn,
								   const Vec3& vMin, const Vec3& vMax,	double& d, Vec3& n)
{
	int bout = 0;
	double d0[6];
	int idx0 = -1;

	// 各軸ごとに最小と最大境界外になっていないか調べる
	int c = 0;
	for(int i = 0; i < 3; ++i){
		int idx = 2*i;
		if((d0[idx] = (spos[i]-r)-vMin[i]) < 0.0){
			bout |= (1 << idx); c++;
			idx0 = idx;
		}
		idx = 2*i+1;
		if((d0[idx] = vMax[i]-(spos[i]+r)) < 0.0){
			bout |= (1 << idx); c++;
			idx0 = idx;
		}
	}

	// AABB内(全軸で境界内)
	if(bout == 0){
		double min_d = 1e10;
		int idx1 = -1;
		for(int i = 0; i < 6; ++i){
			if(d0[i] <= min_d){
				min_d = d0[i];
				idx1 = i;
			}
		}

		d = sgn*min_d;
		n = (idx1 != -1) ? sgn*RX_AABB_NORMALS[idx1] : Vec3(0.0);
		return true;
	}

	// AABB外
	Vec3 x(0.0);
	for(int i = 0; i < 3; ++i){
		if(bout & (1 << (2*i))){
			x[i] = d0[2*i];
		} else if(bout & (1 << (2*i+1))){
			x[i] = -d0[2*i+1];
		}
	}

	// sgn = 1:箱，-1:オブジェクト
	if(c == 1){
		// 平面近傍
		d = sgn*d0[idx0];
		n = sgn*RX_AABB_NORMALS[idx0];
	} else{
		// エッジ/コーナー近傍
		d = -sgn*norm(x);
		n = sgn*(-Unit(x));
	}

	return false;
}


/*!
 * 平面の陰関数値計算
 * @param[in] pos 陰関数値を計算する位置
 * @param[out] d,n 平面までの距離(陰関数値)と法線
 * @param[out] v 平面の速度(今のところ常に0)
 * @param[in] pn,pq 平面の法線と平面上の点
 */
inline bool GetImplicitPlane(const Vec3& pos, double& d, Vec3& n, Vec3& v, const Vec3& pn, const Vec3& pq)
{
	d = dot(pq-pos, pn);
	n = pn;
	v = Vec3(0.0);
	return true;
}


/*!
 * 上が空いたボックス形状の陰関数値計算
 *  - 形状の中心を原点とした座標系での計算
 * @param[in] spos 立方体の中心を原点とした相対座標値
 * @param[in] r    半径(球の場合)
 * @param[in] sgn  立方体の内で距離が正:1,外で正:-1
 * @param[in] sl0  ボックスの内側サイズ(side length)
 * @param[in] sl1  ボックスの外側サイズ(side length)
 * @param[out] d   符号付距離値
 * @param[out] n   最近傍点の法線方向
 */
static inline bool GetImplicitOpenBox(const Vec3& spos, const double& r, const int& sgn,
									  const Vec3& sl0, const Vec3& sl1, double& d, Vec3& n)
{
	int t = 2;
	d = RX_FEQ_INF;

	double td;
	Vec3 tn;

	// 底
	Vec3 m0, m1;
	m0 = -sl1;
	m1 = sl1;
	m1[t] = sl0[t];
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 -x
	m0 = Vec3(-sl1[0], -sl1[1], -sl0[2]);
	m1 = Vec3(-sl0[0], sl1[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 +x
	m0 = Vec3(sl0[0], -sl1[1], -sl0[2]);
	m1 = Vec3(sl1[0], sl1[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 -y
	m0 = Vec3(-sl0[0], -sl1[1], -sl0[2]);
	m1 = Vec3(sl0[0], -sl0[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 +y
	m0 = Vec3(-sl0[0], sl0[1], -sl0[2]);
	m1 = Vec3(sl0[0], sl1[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}
}


/*!
 * 光線(レイ,半直線)と球の交差判定
 * @param[in] p,d レイの原点と方向
 * @param[in] c,r 球の中心と半径
 * @param[out] t1,t2 pから交点までの距離
 * @return 交点数
 */
inline int ray_sphere(const Vec3& p, const Vec3& d, const Vec3& sc, const double r, double& t1, double& t2)
{
	Vec3 q = p-sc;	// 球中心座標系での光線原点座標

	double a = norm2(d);
	double b = 2*dot(q, d);
	double c = norm2(q)-r*r;

	// 判別式
	double D = b*b-4*a*c;

	if(D < 0.0){ // 交差なし
		return 0;
	} else if(D < RX_FEQ_EPS){ // 交点数1
		t1 = -b/(2*a);
		t2 = -1;
		return 1;
	} else{ // 交点数2
		double sqrtD = sqrt(D);
		t1 = (-b-sqrtD)/(2*a);
		t2 = (-b+sqrtD)/(2*a);
		return 2;
	}

}

/*!
 * 三角形と球の交差判定
 * @param[in] v0,v1,v2	三角形の頂点
 * @param[in] n			三角形の法線
 * @param[in] p			最近傍点
 * @return
 */
inline bool triangle_sphere(const Vec3& v0, const Vec3& v1, const Vec3& v2, const Vec3& n,
	const Vec3& c, const double& r, double& dist, Vec3& ipoint)
{
	// ポリゴンを含む平面と球中心の距離
	double d = dot(v0, n);
	double l = dot(n, c)-d;

	dist = l;
	if(l > r) return false;

	// 平面との最近傍点座標
	Vec3 p = c-l*n;

	// 近傍点が三角形内かどうかの判定
	Vec3 n1 = cross((v0-c), (v1-c));
	Vec3 n2 = cross((v1-c), (v2-c));
	Vec3 n3 = cross((v2-c), (v0-c));

	ipoint = p;
	dist = l;
	if(dot(n1, n2) > 0 && dot(n2, n3) > 0){		// 三角形内
		return true;
	} else{		// 三角形外
		// 三角形の各エッジと球の衝突判定
		for(int e = 0; e < 3; ++e){
			Vec3 va0 = (e == 0 ? v0 : (e == 1 ? v1 : v2));
			Vec3 va1 = (e == 0 ? v1 : (e == 1 ? v2 : v0));

			double t1, t2;
			int n = ray_sphere(va0, Unit(va1-va0), c, r, t1, t2);

			if(n){
				double le2 = norm2(va1-va0);
				if((t1 >= 0.0 && t1*t1 < le2) || (t2 >= 0.0 && t2*t2 < le2)){
					return true;
				}
			}
		}
		return false;
	}
}

/*!
 * 線分(を含む直線)と点の距離
 * @param[in] v0,v1 線分の両端点座標
 * @param[in] p 点の座標
 * @return 距離
 */
inline double segment_point_dist(const Vec3& v0, const Vec3& v1, const Vec3& p)
{
	Vec3 v = Unit(v1-v0);
	Vec3 vp = p-v0;
	Vec3 vh = dot(vp, v)*v;
	return norm(vp-vh);
}

/*!
 * 線分と球の交差判定
 * @param[in] s0,s1	線分の端点
 * @param[in] sc,r   球の中心座標と半径
 * @param[out] d2 線分との距離の二乗
 * @return 交差ありでtrue
 */
inline bool segment_sphere(const Vec3& s0, const Vec3& s1, const Vec3& sc, const double& r, double& d2)
{
	Vec3 v = s1-s0;
	Vec3 c = sc-s0;

	double vc = dot(v, c);
	if(vc < 0){		// 球の中心が線分の始点s0の外にある
		d2 = norm2(c);
		return (d2 < r* r);	// 球中心と始点s0の距離で交差判定
	} else{
		double v2 = norm2(v);
		if(vc > v2){	// 球の中心が線分の終点s1の外にある
			d2 = norm2(s1-sc);
			return (d2 < r* r);	// 球中心と終点s1の距離で交差判定
		} else{			// 球がs0とs1の間にある
			d2 = norm2((vc*v)/norm2(v)-c);
			return (d2 < r* r);	// 直線と球中心の距離で交差判定
		}
	}

	return false;
}





//-----------------------------------------------------------------------------
// 異方性カーネル計算用
//-----------------------------------------------------------------------------
inline float RxPythag(const float a, const float b)
{
	float absa = abs(a), absb = abs(b);
	return (absa > absb ? absa*(float)sqrt((double)(1.0+(absb/absa)*(absb/absa))) :
		   (absb == 0.0 ? 0.0 : absb*(float)sqrt((double)(1.0+(absa/absb)*(absa/absb)))));
}

static void RxSVDecomp3(float w[3], float u[9], float v[9], float eps)
{
	bool flag;
	int i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z;
	float rv1[3];
	g = scale = anorm = 0.0;
	for(i = 0; i < 3; ++i){
		l = i+2;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		for(k = i; k < 3; ++k) scale += abs(u[k*3+i]);
		if(scale != 0.0){
			for(k = i; k < 3; ++k){
				u[k*3+i] /= scale;
				s += u[k*3+i]*u[k*3+i];
			}
			f = u[i*3+i];
			g = -RX_SIGN2(sqrt(s), f);
			h = f*g-s;
			u[i*3+i] = f-g;
			for(j = l-1; j < 3; ++j){
				for(s = 0.0, k = i; k < 3; ++k) s += u[k*3+i]*u[k*3+j];
				f = s/h;
				for(k = i; k < 3; ++k) u[k*3+j] += f*u[k*3+i];
			}
			for(k = i; k < 3; ++k) u[k*3+i] *= scale;
		}

		w[i] = scale*g;
		g = s = scale = 0.0;
		if(i+1 <= 3 && i+1 != 3){
			for(k = l-1; k < 3; ++k) scale += abs(u[i*3+k]);
			if(scale != 0.0){
				for(k = l-1; k < 3; ++k){
					u[i*3+k] /= scale;
					s += u[i*3+k]*u[i*3+k];
				}
				f = u[i*3+l-1];
				g = -RX_SIGN2(sqrt(s), f);
				h = f*g-s;
				u[i*3+l-1] = f-g;
				for(k = l-1; k < 3; ++k) rv1[k] = u[i*3+k]/h;
				for(j = l-1; j < 3; ++j){
					for(s = 0.0,k = l-1; k < 3; ++k) s += u[j*3+k]*u[i*3+k];
					for(k = l-1; k < 3; ++k) u[j*3+k] += s*rv1[k];
				}
				for(k = l-1; k < 3; ++k) u[i*3+k] *= scale;
			}
		}
		anorm = RX_MAX(anorm, (abs(w[i])+abs(rv1[i])));
	}
	for(i = 2; i >= 0; --i){
		if(i < 2){
			if(g != 0.0){
				for(j = l; j < 3; ++j){
					v[j*3+i] = (u[i*3+j]/u[i*3+l])/g;
				}
				for(j = l; j < 3; ++j){
					for(s = 0.0, k = l; k < 3; ++k) s += u[i*3+k]*v[k*3+j];
					for(k = l; k < 3; ++k) v[k*3+j] += s*v[k*3+i];
				}
			}
			for(j = l; j < 3; ++j) v[i*3+j] = v[j*3+i] = 0.0;
		}
		v[i*3+i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for(i = 2; i >= 0; --i){
		l = i+1;
		g = w[i];
		for(j = l; j < 3; ++j) u[i*3+j] = 0.0;
		if(g != 0.0){
			g = 1.0/g;
			for(j = l; j < 3; ++j){
				for(s = 0.0, k = l; k < 3; ++k) s += u[k*3+i]*u[k*3+j];
				f = (s/u[i*3+i])*g;
				for(k = i; k < 3; ++k) u[k*3+j] += f*u[k*3+i];
			}
			for(j = i; j < 3; ++j) u[j*3+i] *= g;
		}
		else{
			for(j = i; j < 3; ++j) u[j*3+i] = 0.0;
		}
		++u[i*3+i];
	}
	for(k = 2; k >= 0; --k){
		for(its = 0; its < 30; ++its){
			flag = true;
			for(l = k; l >= 0; --l){
				nm = l-1;
				if(l == 0 || abs(rv1[l]) <= eps*anorm){
					flag = false;
					break;
				}
				if(abs(w[nm]) <= eps*anorm) break;
			}
			if(flag){
				c = 0.0;
				s = 1.0;
				for(i = l; i < k+1; ++i){
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if(abs(f) <= eps*anorm) break;
					g = w[i];
					h = RxPythag(f, g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for(j = 0; j < 3; ++j){
						y = u[j*3+nm];
						z = u[j*3+i];
						u[j*3+nm] = y*c+z*s;
						u[j*3+i] = z*c-y*s;
					}
				}
			}
			z = w[k];
			if(l == k){
				if(z < 0.0){
					w[k] = -z;
					for(j = 0; j < 3; ++j) v[j*3+k] = -v[j*3+k];
				}
				break;
			}
			if(its == 29){
				printf("no convergence in 30 svdcmp iterations");
				return;
			}
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = RxPythag(f, 1.0f);
			f = ((x-z)*(x+z)+h*((y/(f+RX_SIGN2(g, f)))-h))/x;
			c = s = 1.0;
			for(j = l; j <= nm; ++j){
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = RxPythag(f, h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for(jj = 0; jj < 3; ++jj){
					x = v[jj*3+j];
					z = v[jj*3+i];
					v[jj*3+j] = x*c+z*s;
					v[jj*3+i] = z*c-x*s;
				}
				z = RxPythag(f, h);
				w[j] = z;
				if(z){
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for(jj = 0; jj < 3; ++jj){
					y = u[jj*3+j];
					z = u[jj*3+i];
					u[jj*3+j] = y*c+z*s;
					u[jj*3+i] = z*c-y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}

	// reorder
	int inc = 1;
	float sw;
	float su[3], sv[3];

	do{
		inc *= 3;
		inc++; 
	}while(inc <= 3);

	do{
		inc /= 3;
		for(i = inc; i < 3; ++i){
			sw = w[i];
			for(k = 0; k < 3; ++k) su[k] = u[k*3+i];
			for(k = 0; k < 3; ++k) sv[k] = v[k*3+i];
			j = i;
			while (w[j-inc] < sw){
				w[j] = w[j-inc];
				for(k = 0; k < 3; ++k) u[k*3+j] = u[k*3+j-inc];
				for(k = 0; k < 3; ++k) v[k*3+j] = v[k*3+j-inc];
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for(k = 0; k < 3; ++k) u[k*3+j] = su[k];
			for(k = 0; k < 3; ++k) v[k*3+j] = sv[k];

		}
	}while(inc > 1);

	for(k = 0; k < 3; ++k){
		s = 0;
		for(i = 0; i < 3; ++i) if(u[i*3+k] < 0.) s++;
		for(j = 0; j < 3; ++j) if(v[j*3+k] < 0.) s++;
		if(s > 3){
			for(i = 0; i < 3; ++i) u[i*3+k] = -u[i*3+k];
			for(j = 0; j < 3; ++j) v[j*3+k] = -v[j*3+k];
		}
	}
}



//-----------------------------------------------------------------------------
// 数値計算
//-----------------------------------------------------------------------------
/*!
* 2x2行列の行列式の計算
*  | m[0] m[1] |
*  | m[2] m[3] |
* @param[in] m 元の行列
* @return 行列式の値
*/
static inline double CalDetMat2x2(const double *m)
{
	return m[0]*m[3]-m[1]*m[2];
}

/*!
* 2x2行列の逆行列の計算
*  | m[0] m[1] |
*  | m[2] m[3] |
* @param[in] m 元の行列
* @param[out] invm 逆行列
* @return 逆行列の存在
*/
static inline bool CalInvMat2x2(const double *m, double *invm)
{
	double det = CalDetMat2x2(m);
	if(fabs(det) < RX_FEQ_EPS){
		return false;
	}
	else{
		double inv_det = 1.0/det;
		invm[0] =  inv_det*m[3];
		invm[1] = -inv_det*m[1];
		invm[2] = -inv_det*m[2];
		invm[3] =  inv_det*m[0];
		return true;
	}
}

/*!
* 3x3行列の行列式の計算
*  | m[0] m[1] m[2] |
*  | m[3] m[4] m[5] |
*  | m[6] m[7] m[8] |
* @param[in] m 元の行列
* @return 行列式の値
*/
static inline double CalDetMat3x3(const double *m)
{
	// a11a22a33+a21a32a13+a31a12a23-a11a32a23-a31a22a13-a21a12a33
	return m[0]*m[4]*m[8]+m[3]*m[7]*m[2]+m[6]*m[1]*m[5]
		-m[0]*m[7]*m[5]-m[6]*m[4]*m[2]-m[3]*m[1]*m[8];
}

/*!
* 3x3行列の逆行列の計算
*  | m[0] m[1] m[2] |
*  | m[3] m[4] m[5] |
*  | m[6] m[7] m[8] |
* @param[in] m 元の行列
* @param[out] invm 逆行列
* @return 逆行列の存在
*/
static inline bool CalInvMat3x3(const double *m, double *invm)
{
	double det = CalDetMat3x3(m);
	if(fabs(det) < RX_FEQ_EPS){
		return false;
	}
	else{
		double inv_det = 1.0/det;

		invm[0] = inv_det*(m[4]*m[8]-m[5]*m[7]);
		invm[1] = inv_det*(m[2]*m[7]-m[1]*m[8]);
		invm[2] = inv_det*(m[1]*m[5]-m[2]*m[4]);

		invm[3] = inv_det*(m[5]*m[6]-m[3]*m[8]);
		invm[4] = inv_det*(m[0]*m[8]-m[2]*m[6]);
		invm[5] = inv_det*(m[2]*m[3]-m[0]*m[5]);

		invm[6] = inv_det*(m[3]*m[7]-m[4]*m[6]);
		invm[7] = inv_det*(m[1]*m[6]-m[0]*m[7]);
		invm[8] = inv_det*(m[0]*m[4]-m[1]*m[3]);

		return true;
	}
}



#endif // #ifndef _RX_SPH_COMMON_H_