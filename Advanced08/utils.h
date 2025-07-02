/*!
  @file rx_utils.h

  @brief 各種便利な関数や定数など

  @author Makoto Fujisawa
  @date 2020-07
*/
// FILE --rx_utils.h--

#ifndef _RX_UTILS_H_
#define _RX_UTILS_H_

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// C標準
#include <ctime>
#include <cmath>
//#include <cctype>
#include <cstdio>
#include <cassert>

#ifdef WIN32
#include <direct.h>
#endif

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

// glm
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/ext.hpp"	// for glm::to_string()

// 画像読み書き
#include "rx_bitmap.h"

// OpenGLのシェーダ
#include "rx_shaders.h"

using namespace std;


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
const double RX_PI = 3.14159265358979323846;	// 円周率
const double RX_FEQ_EPS = 1.0e-8;				// 許容誤差
const double RX_FEQ_INF = 1.0e10;				// 大きい数の初期値設定用
const double RX_DEGREES_TO_RADIANS = 0.0174532925199432957692369076848;	//! degree -> radian の変換係数(pi/180.0)
const double RX_RADIANS_TO_DEGREES = 57.295779513082320876798154814114;	//! radian -> degree の変換係数(180.0/pi)

//! 許容誤差を含めた等値判定
template<class T>
inline bool RX_FEQ(const T &a, const T &b){ return (fabs(a-b) < RX_FEQ_EPS); }

//! 値のクランプ(クランプした値を返す)
template<class T>
inline T RX_CLAMP(const T &x, const T &a, const T &b){ return ((x < a) ? a : (x > b) ? b : x); }

//! 1次元線型補間
template<class T>
inline T RX_LERP(const T &a, const T &b, const T &t){ return a + t*(b-a); }

//! 乱数
inline double RX_RAND(const double &_min, const double &_max){ return (_max-_min)*(double(rand())/(1.0+RAND_MAX))+_min; }

//! スワップ
template<class T>
inline void RX_SWAP(T &a, T &b){ T c; c = a; a = b; b = c; }


//-----------------------------------------------------------------------------
// VAO関連
//-----------------------------------------------------------------------------
/*!
 * 頂点配列オブジェクトの作成
 * - 最低限頂点情報があれば良い．必要のない引数は0を入れておく．
 * @param[in] vrts,nvrts 頂点座標配列と頂点数
 * @param[in] dim 次元(2or3)
 * @param[in] tris,ntris 三角形ポリゴンを構成する頂点インデックス列を格納した配列と三角形ポリゴン数
 * @param[in] nrms,nnrms 各頂点の法線情報配列と法線数(=頂点数)
 * @param[in] cols,ncols 各頂点の色情報配列と色情報数(=頂点数)
 * @param[in] tcds,ntcds 各頂点のテクスチャ座標情報配列とテクスチャ座標数(=頂点数)
 * @return VAOオブジェクト
 */
static inline GLuint CreateVAO(
	const GLfloat *vrts, GLuint nvrts, int dim = 3,
	const GLuint  *tris = 0, GLuint ntris = 0,
	const GLfloat *nrms = 0, GLuint nnrms = 0,
	const GLfloat *cols = 0, GLuint ncols = 0,
	const GLfloat *tcds = 0, GLuint ntcds = 0)
{
	// VAOの生成
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// VBO:頂点座標
	GLuint vvbo;
	glGenBuffers(1, &vvbo);
	glBindBuffer(GL_ARRAY_BUFFER, vvbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*dim*nvrts, vrts, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, dim, GL_FLOAT, GL_FALSE, 0, 0);

	// VBO:頂点テクスチャ座標
	if(tcds){
		GLuint tvbo;
		glGenBuffers(1, &tvbo);
		glBindBuffer(GL_ARRAY_BUFFER, tvbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*2*ntcds, tcds, GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
	}

	// VBO:頂点法線
	if(nrms){
		GLuint nvbo;
		glGenBuffers(1, &nvbo);
		glBindBuffer(GL_ARRAY_BUFFER, nvbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*dim*nnrms, nrms, GL_STATIC_DRAW);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, dim, GL_FLOAT, GL_FALSE, 0, 0);
	}

	// VBO:頂点描画色
	if(cols){
		GLuint cvbo;
		glGenBuffers(1, &cvbo);
		glBindBuffer(GL_ARRAY_BUFFER, cvbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*3*ncols, cols, GL_STATIC_DRAW);
		glEnableVertexAttribArray(3);
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);
	}

	// VBO:三角形ポリゴンインデックス
	if(tris){
		GLuint pvbo;
		glGenBuffers(1, &pvbo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pvbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int)*ntris*3, tris, GL_STATIC_DRAW);
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	return vao;
}

/*!
 * VAOに関連づけられたデータの更新
 * @param[in] vao VAOオブジェクト
 * @param[in] index glVertexAttribPointerで関連づけられたインデックス
 * @param[in] data データ配列
 * @param[in] n データ数
 * @return VBOオブジェクト
 */
static inline GLuint UpdateDataVAO(GLuint vao, GLuint index, const GLfloat *data, GLuint n)
{
	glBindVertexArray(vao);
	GLuint vbo;
	glGetVertexAttribIuiv(index, GL_VERTEX_ATTRIB_ARRAY_BUFFER_BINDING, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*n, data);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	return vbo;
}

/*!
 * VAOの破棄
 *  - リンクされているVBOなどは破棄されないようなのでそれらは個別に行う必要あり
 * @param[in] vao VAOオブジェクト
 */
static inline void DeleteVAO(GLuint vao)
{
	glDeleteVertexArrays(1, &vao);
}



/*!
 * 球体形状のポリゴンメッシュを生成
 * - 球の中心は原点(0,0,0)
 * @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
 * @param[in] rad 球の半径
 * @param[in] slices,stacks 緯度方向(360度)と傾度方向(180度)のポリゴン分割数
 * @return 生成したVAOオブジェクト
 */
static inline int MakeSphere(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float rad, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();
	nvrts = 0;
	ntris = 0;
	vrts.clear();
	nrms.clear();
	tris.clear();

	for(int j = 0; j <= stacks; ++j){
		float t = float(j)/float(stacks);
		float y = rad*cos(pi*t);
		float rj = rad*sin(pi*t);	// 高さyでの球の断面円半径
		for(int i = 0; i <= slices; ++i){
			float s = float(i)/float(slices);
			float x = rj*sin(2*pi*s);
			float z = rj*cos(2*pi*s);
			vrts.push_back(glm::vec3(x, y, z));
			nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
		}
	}
	nvrts = static_cast<int>(vrts.size());
	// メッシュ作成
	int nx = slices+1;
	for(int j = 0; j < stacks; ++j){
		for(int i = 0; i < slices; ++i){
			tris.push_back((i)+(j)*nx);
			tris.push_back((i+1)+(j+1)*nx);
			tris.push_back((i+1)+(j)*nx);

			tris.push_back((i)+(j)*nx);
			tris.push_back((i)+(j+1)*nx);
			tris.push_back((i+1)+(j+1)*nx);
		}
	}
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}

/*!
 * 円筒形状のポリゴンメッシュを生成してVAOとして登録
 * - 円筒の中心は原点(0,0,0)
 * - 円筒の軸方向はz軸方向(0,0,1) - gluCylinderに合わせている
 * @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
 * @param[in] rad,len 円筒の半径と長さ
 * @param[in] slices 円筒の円に沿ったポリゴン分割数
 * @return 生成したVAOオブジェクト
 */
static inline int MakeCylinder(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float rad, float len, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();

	for(int i = 0; i <= slices; ++i){
		float t = float(i)/float(slices);
		float x = rad*cos(2*pi*t);
		float y = rad*sin(2*pi*t);
		vrts.push_back(glm::vec3(x, y, -0.5*len));
		nrms.push_back(glm::normalize(glm::vec3(x, y, 0.0)));
		vrts.push_back(glm::vec3(x, y, 0.5*len));
		nrms.push_back(glm::normalize(glm::vec3(x, y, 0.0)));

	}
	nvrts = static_cast<int>(vrts.size());

	// メッシュ作成
	for(int i = 0; i < 2*slices; i += 2){
		tris.push_back(i);
		tris.push_back((i+2 >= 2*slices ? 0 : i+2));
		tris.push_back(i+1);

		tris.push_back(i+1);
		tris.push_back((i+2 >= 2*slices ? 0 : i+2));
		tris.push_back((i+2 >= 2*slices ? 1 : i+3));
	}
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}


/*!
 * カプセル形状(円筒の両端が半球)のポリゴンメッシュを生成してVAOとして登録
 * - 形状の中心は原点(0,0,0)
 * - カプセル形状の軸方向はz軸方向(0,0,1) - gluCylinderに合わせている
 * @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
 * @param[in] rad,len 円筒の半径と長さ
 * @param[in] slices 円筒の円に沿ったポリゴン分割数
 * @return 生成したVAOオブジェクト
 */
static inline int MakeCapsule(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float rad, float len, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();
	int offset = 0;

	// 球体の中心断面(赤道部分)に沿った頂点を2重にして，
	// その部分を円筒長さ分z方向に伸ばすことでカプセル形状を生成
	for(int j = 0; j <= stacks; ++j){
		float t = float(j)/float(stacks);
		float z = rad*cos(pi*t);
		float rj = rad*sin(pi*t);	// 高さyでの球の断面円半径
		float zlen = (j < stacks/2 ? 0.5*len : -0.5*len);	// z方向のオフセット量

		if(j == stacks/2){	// 中心断面(赤道部分)に頂点を追加
			for(int i = 0; i <= slices; ++i){
				float s = float(i)/float(slices);
				float x = rj*sin(2*pi*s);
				float y = rj*cos(2*pi*s);
				vrts.push_back(glm::vec3(x, y, z-zlen));
				nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
			}
		}
		for(int i = 0; i <= slices; ++i){
			float s = float(i)/float(slices);
			float x = rj*sin(2*pi*s);
			float y = rj*cos(2*pi*s);
			vrts.push_back(glm::vec3(x, y, z+zlen));
			nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
		}
	}

	// メッシュ作成
	int nx = slices+1;
	for(int j = 0; j < stacks+1; ++j){
		for(int i = 0; i < slices; ++i){
			tris.push_back((i)+(j)*nx+offset);
			tris.push_back((i+1)+(j)*nx+offset);
			tris.push_back((i+1)+(j+1)*nx+offset);

			tris.push_back((i)+(j)*nx+offset);
			tris.push_back((i+1)+(j+1)*nx+offset);
			tris.push_back((i)+(j+1)*nx+offset);
		}
	}

	nvrts = static_cast<int>(vrts.size());
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}



/*!
 * VAOによるプリミティブ描画
 * @param[in] obj 頂点数,ポリゴン数情報を含むVAO
 * @param[in] draw 描画フラグ
 */
inline static void DrawPrimitive(const GLuint vao, const int nvrts, const int ntris, int draw)
{
	// エッジ描画における"stitching"をなくすためのオフセットの設定
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	// 図形の描画
	glDisable(GL_LIGHTING);
	glBindVertexArray(vao);
	if(draw & 0x01){
		glColor3d(1.0, 0.3, 0.3);
		glPointSize(5.0);
		glDrawArrays(GL_POINTS, 0, nvrts);
	}
	if(draw & 0x02){
		glColor3d(0.5, 0.9, 0.9);
		glLineWidth(4.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, 0);
	}
	if(draw & 0x04){
		glEnable(GL_LIGHTING);
		//glDisable(GL_CULL_FACE);
		glEnable(GL_AUTO_NORMAL);
		glEnable(GL_NORMALIZE);
		glColor3d(0.1, 0.5, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, 0);
	}
	glBindVertexArray(0);
}






//-----------------------------------------------------------------------------
// glm関係の関数
//-----------------------------------------------------------------------------
inline int IDX4(int row, int col){ return (row | (col<<2)); }

/*!
 * 1次元配列に格納された4x4行列とglm::vec3ベクトルの掛け算
 *  | m[0]  m[1]  m[2]  m[3]  |
 *  | m[4]  m[5]  m[6]  m[7]  |
 *  | m[8]  m[9]  m[10] m[11] |
 *  | m[12] m[13] m[14] m[15] |
 * @param[in] m 元の行列
 * @param[in] v ベクトル
 * @return 積の結果のベクトル(m*v)
 */
inline glm::vec3 MulMat4Vec3(GLfloat* m, const glm::vec3& v)
{
	glm::vec3 d(0.0);
	d[0] = (v[0]*m[IDX4(0, 0)]+v[1]*m[IDX4(0, 1)]+v[2]*m[IDX4(0, 2)]);
	d[1] = (v[0]*m[IDX4(1, 0)]+v[1]*m[IDX4(1, 1)]+v[2]*m[IDX4(1, 2)]);
	d[2] = (v[0]*m[IDX4(2, 0)]+v[1]*m[IDX4(2, 1)]+v[2]*m[IDX4(2, 2)]);
	return d;
}



//-----------------------------------------------------------------------------
// 文字列処理
//-----------------------------------------------------------------------------
/*!
 * 様々な型のstringへの変換(stringstreamを使用)
 * @param[in] x 入力
 * @return string型に変換したもの
 */
template<typename T>
inline std::string RX_TO_STRING(const T &x)
{
	std::stringstream ss;
	ss << x;
	return ss.str();
}

//! string型に<<オペレータを設定
template<typename T>
inline std::string &operator<<(std::string &cb, const T &a)
{
	cb += RX_TO_STRING(a);
	return cb;
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

/*!
 * "(x, y, z)"の形式の文字列からglm::vec3型へ変換
 *  - (x)となっていたら(x, x, x)とする．
 * @param[in] s 文字列
 * @param[out] v 値
 * @return 要素記述数
 */
inline int StringToVec3(const string &s, glm::vec3 &v)
{
	int vcount = 0;
	size_t pos;
	v = glm::vec3(0.0);
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
 * 整数値の下一桁を返す
 * @param[in] x 整数値
 * @return xの下一桁
 */
inline int LowOneDigit(const int &x)
{
	int x1 = (x < 0) ? -x : x;
	float a = 10;

	//INT_MAX = 2147483647
	for(int i = 10; i >= 1; i--){
		a = pow(10.0, (float)i);
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
		dn = (int)(log10((float)n))+1;
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
inline string GenTimeString(float sec, bool use_msec = false)
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
	if(i == string::npos) {
		sub = src.substr(pos, string::npos);
		return (int)string::npos;
	}
	else {
		int cnt = 1;
		while(src[i + cnt] == ' ') {    // 区切り文字の後のスペースを消す
			cnt++;
		}
		sub = src.substr(pos, i - pos);
		return (int)(i + cnt >= src.size() ? (int)string::npos : i + cnt);
	}
}


//-----------------------------------------------------------------------------
// ファイル・フォルダ処理
//-----------------------------------------------------------------------------
/*!
 * ファイル，フォルダの存在確認
 * @param[in] path ファイル・フォルダパス
 */
inline int ExistFile(const string fn)
{
	FILE *fp;

	if((fp = fopen(fn.c_str(), "r")) == NULL){
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
		dn = (int)(log10((float)n))+1;
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
 * ディレクトリ作成(多階層対応) - Windows only
 * @param[in] dir 作成ディレクトリパス
 * @return 成功で1,失敗で0 (ディレクトリがすでにある場合も1を返す)
 */
static int MkDir(string dir)
{
#ifdef WIN32
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
#else
	return 0;
#endif
}

//! ファイルパス検索用(made by 金森先生)
// [How to use]
// PathFinder p;
// p.addSearchPath("bin");
// p.addSearchPath("../bin");
// p.addSearchPath("../../bin");
// std::string filename = p.find("sample.bmp");
class PathFinder
{
public:
	void addSearchPath(const std::string& dirname) { m_SearchDirs.push_back(dirname); }

	std::string find(const std::string& filename, const std::string& separator = "/")
	{
		FILE* fp = 0;
		std::string path_name;

		for(unsigned int i = 0; i < m_SearchDirs.size(); i++)
		{
			path_name = m_SearchDirs[i] + separator + filename;
#ifdef _WIN32
			::fopen_s(&fp, path_name.c_str(), "r");
#else
			fp = fopen(path_name.c_str(), "r");
#endif

			if(fp != 0)
			{
				fclose(fp);
				return path_name;
			}
		}

		path_name = filename;
#ifdef _WIN32
		::fopen_s(&fp, path_name.c_str(), "r");
#else
		fp = fopen(path_name.c_str(), "r");
#endif

		if(fp != 0)
		{
			fclose(fp);
			return path_name;
		}

		printf("file not found: %s\n", filename.c_str());

		return "";
	}

private:
	std::vector<std::string> m_SearchDirs;
};

//! OpenGLのエラーチェック(made by 金森先生)
// see https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glGetError.xhtml
static inline void CheckGLError(const char *func, const char *file, int line)
{
	GLenum errCode = glGetError();

	if(errCode != GL_NO_ERROR)
	{
		std::cerr << func << ": OpenGL Error: ";

		switch(errCode)
		{
		case GL_INVALID_ENUM:
			std::cerr << "GL_INVALID_ENUM: An unacceptable value is specified for an enumerated argument";
			break;
		case GL_INVALID_VALUE:
			std::cerr << "GL_INVALID_VALUE: A numeric argument is out of range";
			break;
		case GL_INVALID_OPERATION:
			std::cerr << "GL_INVALID_OPERATION: The specified operation is not allowed in the current state";
			break;
		case GL_STACK_OVERFLOW:
			std::cerr << "GL_STACK_OVERFLOW: This command would cause a stack overflow";
			break;
		case GL_STACK_UNDERFLOW:
			std::cerr << "GL_STACK_UNDERFLOW: This command would cause a stack underflow";
			break;
		case GL_OUT_OF_MEMORY:
			std::cerr << "GL_OUT_OF_MEMORY: There is not enough memory left to execute the command";
			break;
		case GL_TABLE_TOO_LARGE:
			std::cerr << "GL_TABLE_TOO_LARGE: The specified table exceeds the implementation's maximum supported table size";
			break;
		default:
			std::cerr << "Unknown command, error code = " << std::showbase << std::hex << errCode;
			break;
		}

		std::cerr << " (file: " << file << " at line " << line << std::endl;
	}
}
#define CHECK_GL_ERROR CheckGLError(__FUNCTION__, __FILE__, __LINE__)


//-----------------------------------------------------------------------------
// テクスチャ関連
//-----------------------------------------------------------------------------
/*!
 * 画像読込み -> OpenGLテクスチャ登録
 * @param[in] fn ファイル名
 * @param[inout] tex_name テクスチャ名(0なら新たに生成)
 * @param[in] mipmap ミップマップ使用フラグ
 * @param[in] compress テクスチャ圧縮使用フラグ
 */
static int loadTexture(const string& fn, GLuint& tex_name, bool mipmap, bool compress)
{
	// 画像読み込み
	int w, h, c, wstep;
	unsigned char* pimg;
	pimg = ReadBitmapFile(fn, w, h, c, wstep, false, true);
	if(!pimg){
		return 0;
	}

	GLuint iformat, format;

	// 画像フォーマット
	format = GL_RGBA;
	if(c == 1){
		format = GL_LUMINANCE;
	}
	else if(c == 3){
		format = GL_RGB;
	}

	// OpenGL内部の格納フォーマット
	if(compress){
		iformat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT;
		if(c == 1){
			iformat = GL_COMPRESSED_LUMINANCE_ARB;
		}
		else if(c == 3){
			iformat = GL_COMPRESSED_RGB_S3TC_DXT1_EXT;
		}
	}
	else{
		iformat = GL_RGBA;
		if(c == 1){
			iformat = GL_LUMINANCE;
		}
		else if(c == 3){
			iformat = GL_RGB;
		}
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// テクスチャ作成
	if(tex_name == 0){
		glGenTextures(1, &tex_name);

		glBindTexture(GL_TEXTURE_2D, tex_name);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR));
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		if(mipmap){
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 6);
		}

		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}
	else{
		glBindTexture(GL_TEXTURE_2D, tex_name);
		//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, format, GL_UNSIGNED_BYTE, pimg);
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}

	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] pimg;

	return 1;
}

/*!
 * 画像読込み -> OpenGLテクスチャ登録
 * @param[in] fn ファイル名
 * @param[inout] tex_name テクスチャ名(0なら新たに生成)
 * @param[in] mipmap ミップマップ使用フラグ
 * @param[in] compress テクスチャ圧縮使用フラグ
 */
static int makeCheckerBoardTexture(GLuint& tex_name, const int size = 72, bool mipmap = false, bool compress = false)
{
	// ピクセルデータ生成
	unsigned char* pimg = new unsigned char[size*size*4];
	if(!pimg) return 0;

	int cw = size/36;

	// チェッカーボードパターン生成
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			if((i+j%cw)%cw == 0){
				// dark
				pimg[i*size*4 + j*4 + 0] = static_cast<GLubyte>(195);
				pimg[i*size*4 + j*4 + 1] = static_cast<GLubyte>(188);
				pimg[i*size*4 + j*4 + 2] = static_cast<GLubyte>(207);
				pimg[i*size*4 + j*4 + 3] = static_cast<GLubyte>(255);
			}
			else {
				// light
				pimg[i*size*4 + j*4 + 0] = static_cast<GLubyte>(220);
				pimg[i*size*4 + j*4 + 1] = static_cast<GLubyte>(213);
				pimg[i*size*4 + j*4 + 2] = static_cast<GLubyte>(232);
				pimg[i*size*4 + j*4 + 3] = static_cast<GLubyte>(255);
			}
		}
	}


	GLuint iformat, format;

	// 画像フォーマット
	format = GL_RGBA;

	// OpenGL内部の格納フォーマット
	if(compress){
		iformat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT;
	}
	else{
		iformat = GL_RGBA;
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	// テクスチャ作成
	if(tex_name == 0){
		glGenTextures(1, &tex_name);
		glBindTexture(GL_TEXTURE_2D, tex_name);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (mipmap ? GL_NEAREST_MIPMAP_NEAREST : GL_NEAREST));

		if(mipmap){
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 6);
		}

		glTexImage2D(GL_TEXTURE_2D, 0, iformat, size, size, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}
	else{
		glBindTexture(GL_TEXTURE_2D, tex_name);
		//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size, size, format, GL_UNSIGNED_BYTE, pimg);
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, size, size, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}

	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] pimg;

	return 1;
}

/*!
 * フレームバッファのRGB情報を一時的なバッファに確保
 * @retval true 保存成功
 * @retval false 保存失敗
 */
static inline bool saveFrameBuffer(string fn, int w, int h)
{
	string ext = GetExtension(fn);
	if(ext != "bmp") fn += ".bmp";

	int c = 3;
	int wstep = (((w+1)*c)/4)*4;
	vector<unsigned char> imm_buf(wstep*h);

	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &imm_buf[0]);
	WriteBitmapFile(fn, &imm_buf[0], w, h, c, RX_BMP_WINDOWS_V3, wstep, false, true);
	return true;
}




//-----------------------------------------------------------------------------
// 各種オブジェクト描画関数
//-----------------------------------------------------------------------------
/*!
 * 床描画
 *  - GLSLでテクスチャ付き＆柔らかいスポットライトで照らされているような床を描画
 * @param[in] light_pos,light_color 光源位置と色
 * @param[in] y,s 床の高さと水平方向の長さ
 * @return
 */
static inline void drawFloor(glm::vec3 light_pos, glm::vec3 light_color, float y = -1.0, float s = 20.0)
{
	static bool init = true;
	static rxGLSL glslFloor;			//!< GLSLを使った描画
	static GLuint texFloor = 0;				//!< 床のテクスチャ
	if(init){
#ifdef WIN32
		// シェーダの初期化(最初の呼び出し時だけ実行)
		glslFloor = CreateGLSLFromFile("shaders/floor.vp", "shaders/floor.fp", "floor");
#endif
		init = false;

		// 床テクスチャ読み込み
		glActiveTexture(GL_TEXTURE0);
		//loadTexture("floortile.bmp", texFloor, true, false);
		makeCheckerBoardTexture(texFloor, 72, true, false);
		glBindTexture(GL_TEXTURE_2D, texFloor);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16.0f);

	}

	glEnable(GL_TEXTURE_2D);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
#ifdef WIN32
	// 視点座標系のモデルビュー行列取得と視点座標系での光源位置計算
	GLfloat eye_modelview[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, eye_modelview);
	glm::vec3 light_pos_eye = MulMat4Vec3(eye_modelview, light_pos);

	// GLSLシェーダをセット
	glUseProgram(glslFloor.Prog);

	glUniform3f(glGetUniformLocation(glslFloor.Prog, "v3LightPosEye"), light_pos_eye[0], light_pos_eye[1], light_pos_eye[2]);
	glUniform3f(glGetUniformLocation(glslFloor.Prog, "v3LightColor"), light_color[0], light_color[1], light_color[2]);
#endif 

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texFloor);
	glUniform1i(glGetUniformLocation(glslFloor.Prog, "texFloor"), 0);

	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);	// 裏面をカリング

	// 床用ポリゴン描画
	glColor3d(1.0, 1.0, 1.0);
	glNormal3d(0.0, 1.0, 0.0);
	glBegin(GL_QUADS);
	{
		glTexCoord2d(0.0, 0.0);	glVertex3d(-s, y, -s);
		glTexCoord2d(0.0, 1.0);	glVertex3d(-s, y, s);
		glTexCoord2d(1.0, 1.0);	glVertex3d(s, y, s);
		glTexCoord2d(1.0, 0.0);	glVertex3d(s, y, -s);
	}
	glEnd();

	glUseProgram(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}

/*!
 * xyz軸描画(x軸:赤,y軸:緑,z軸:青)
 * @param[in] len 軸の長さ
 */
static inline int drawAxis(float len, float line_width)
{
	glLineWidth((GLfloat)line_width);

	// x軸
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(len, 0.0, 0.0);
	glEnd();

	// y軸
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, len, 0.0);
	glEnd();

	// z軸
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, len);
	glEnd();

	return 1;
}


/*!
 * 直方体のワイヤーフレーム描画
 * @param[in] min 最小座標値
 * @param[in] max 最大座標値
 * @param[in] color 描画色
 */
static inline void drawWireCuboid(const glm::vec3& min, const glm::vec3& max, const glm::vec3& color, float line_width = 1.0)
{
	glLineWidth((GLfloat)line_width);

	glPushMatrix();
	glColor3d(color[0], color[1], color[2]);

	glBegin(GL_LINES);

	// x軸平行
	glVertex3f(min[0], min[1], min[2]);	glVertex3f(max[0], min[1], min[2]);
	glVertex3f(min[0], min[1], max[2]); glVertex3f(max[0], min[1], max[2]);
	glVertex3f(min[0], max[1], min[2]);	glVertex3f(max[0], max[1], min[2]);
	glVertex3f(min[0], max[1], max[2]);	glVertex3f(max[0], max[1], max[2]);

	// z軸平行
	glVertex3f(min[0], min[1], min[2]);	glVertex3f(min[0], min[1], max[2]);
	glVertex3f(min[0], max[1], min[2]);	glVertex3f(min[0], max[1], max[2]);
	glVertex3f(max[0], min[1], min[2]);	glVertex3f(max[0], min[1], max[2]);
	glVertex3f(max[0], max[1], min[2]);	glVertex3f(max[0], max[1], max[2]);

	// z軸平行
	glVertex3f(min[0], min[1], min[2]);	glVertex3f(min[0], max[1], min[2]);
	glVertex3f(min[0], min[1], max[2]);	glVertex3f(min[0], max[1], max[2]);
	glVertex3f(max[0], min[1], min[2]);	glVertex3f(max[0], max[1], min[2]);
	glVertex3f(max[0], min[1], max[2]);	glVertex3f(max[0], max[1], max[2]);

	glEnd();

	glPopMatrix();
}



/*!
 * 原点中心の円のワイヤーフレーム描画
 * @param rad 円の半径
 * @param n 分割数
 */
static inline void DrawWireCircle(const float& rad, const int& n)
{
	float t = 0.0;
	float dt = 2.0*RX_PI/(float)n;

	glBegin(GL_LINE_LOOP);
	do{
		glVertex3f(rad*cos(t), rad*sin(t), 0.0);
		t += dt;
	} while(t < 2.0*RX_PI);
	glEnd();
}

/*!
 * 原点中心の円のワイヤーフレーム描画(XZ平面)
 * @param rad 円の半径
 * @param n 分割数
 */
static inline void DrawWireCircleXZ(const float& rad, const int& n)
{
	float t = 0.0;
	float dt = 2.0*RX_PI/(float)n;

	glBegin(GL_LINE_LOOP);
	do{
		glVertex3f(rad*cos(t), 0.0, rad*sin(t));
		t += dt;
	} while(t < 2.0*RX_PI);
	glEnd();
}

/*!
 * 球のワイヤーフレーム描画
 * @param cen 球の中心
 * @param rad 球の半径
 * @param col 描画色
 */
static inline void DrawWireSphere(const glm::vec3& cen, const float& rad, const glm::vec3& col)
{
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glTranslatef(cen[0], cen[1], cen[2]);
	glRotatef(90, 1.0, 0.0, 0.0);
	glColor3f(col[0], col[1], col[2]);

	// 緯度(x-y平面に平行)
	float z, dz;
	dz = 2.0*rad/8.0f;
	z = -(rad-dz);
	do{
		glPushMatrix();
		glTranslatef(0.0, 0.0, z);
		DrawWireCircle(sqrt(rad*rad-z*z), 32);
		glPopMatrix();
		z += dz;
	} while(z < rad);

	// 経度(z軸まわりに回転)
	float t, dt;
	t = 0.0f;
	dt = 180.0/8.0;
	do{
		glPushMatrix();
		glRotatef(t, 0.0, 0.0, 1.0);
		DrawWireCircleXZ(rad, 32);
		glPopMatrix();

		t += dt;
	} while(t < 180);

	//glutWireSphere(rad, 10, 5);
	glPopMatrix();
}

/*!
 * 球描画
 * @param cen 球の中心
 * @param rad 球の半径
 * @param col 描画色
 */
static inline void drawSolidSphere(const glm::vec3& cen, const float& rad, const glm::vec3& col)
{
	glPushMatrix();
	glTranslatef(cen[0], cen[1], cen[2]);
	glRotatef(90, 1.0, 0.0, 0.0);
	glColor3f(col[0], col[1], col[2]);
	//glutSolidSphere(rad, 20, 10);
	glPopMatrix();
}

/*!
 * 直方体の8頂点座標を中心と辺の長さの1/2(side length)から計算
 * @param[in] cn,sl 立方体の中心と辺の長さの1/2(side length)
 * @param[out] v 8頂点座標
 */
inline void SetVerticesBox(const glm::vec3& cn, const glm::vec3& sl, glm::vec3 v[8])
{
	v[0] = cn+glm::vec3(-sl[0], -sl[1], -sl[2]);
	v[1] = cn+glm::vec3(-sl[0], sl[1], -sl[2]);
	v[2] = cn+glm::vec3(-sl[0], sl[1], sl[2]);
	v[3] = cn+glm::vec3(-sl[0], -sl[1], sl[2]);

	v[4] = cn+glm::vec3(sl[0], -sl[1], -sl[2]);
	v[5] = cn+glm::vec3(sl[0], sl[1], -sl[2]);
	v[6] = cn+glm::vec3(sl[0], sl[1], sl[2]);
	v[7] = cn+glm::vec3(sl[0], -sl[1], sl[2]);
}

/*!
 * 上が空いたボックス形状(壁面厚みあり)を表すポリゴンデータの生成
 * @param[in] sl0,sl1 ボックスの内側と外側のサイズ(side length)
 * @param[in] d 空ける方向
 * @param[out] vrts,idxs ポリゴン頂点座標と接続関係
 */
inline void CreateBoxPolygon(const glm::vec3& sl0, const glm::vec3& sl1, const int& d,
	vector<glm::vec3>& vrts, vector< vector<int> >& idxs)
{
	if(d < 0 || d > 2) return;

	float h = sl1[d]-sl0[d];
	glm::vec3 cn(0.0);

	vrts.resize(16);

	// 外側の頂点
	SetVerticesBox(cn, sl1, &vrts[0]);

	// 内側の頂点
	cn[d] += h;
	SetVerticesBox(cn, sl0, &vrts[8]);


	int idxs0[5][4] = { {0, 3, 2, 1},
						{1, 2, 6, 5},
						{5, 6, 7, 4},
						{4, 7, 3, 0},
						{0, 1, 5, 4} };

	int idxs1[4][4] = { {2, 3, 11, 10},
						{3, 7, 15, 11},
						{7, 6, 14, 15},
						{6, 2, 10, 14} };

	// 三角形作成
	idxs.resize(28);
	for(int i = 0; i < 28; ++i) idxs[i].resize(3);

	int c = 0;

	// 外側の箱
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][j];
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][((j+2 > 3) ? 0 : j+2)];
		}
		c++;
	}

	// 内側の箱
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][2-j]+8;
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][(((2-j)+2 > 3) ? 0 : (2-j)+2)]+8;
		}
		c++;
	}

	// 上部
	for(int i = 0; i < 4; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs1[i][j];
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs1[i][((j+2 > 3) ? 0 : j+2)];
		}
		c++;
	}

}


//-----------------------------------------------------------------------------
// 衝突処理/距離計算用関数
//-----------------------------------------------------------------------------
//! AABBの法線方向(内向き)
const glm::vec3 RX_AABB_NORMALS[6] = { glm::vec3(1.0,  0.0,  0.0),	// x-
								  glm::vec3(-1.0,  0.0,  0.0),	// x+
								  glm::vec3(0.0,  1.0,  0.0),	// y-
								  glm::vec3(0.0, -1.0,  0.0),	// y+
								  glm::vec3(0.0,  0.0,  1.0),	// z-
								  glm::vec3(0.0,  0.0, -1.0) };	// z+

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
static inline bool AABB_point_dist(const glm::vec3& spos, const float& r, const int& sgn,
	const glm::vec3& vMin, const glm::vec3& vMax, float& d, glm::vec3& n)
{
	int bout = 0;
	float d0[6];
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
		float min_d = 1e10;
		int idx1 = -1;
		for(int i = 0; i < 6; ++i){
			if(d0[i] <= min_d){
				min_d = d0[i];
				idx1 = i;
			}
		}

		d = sgn*min_d;
		n = (idx1 != -1) ? float(sgn)*RX_AABB_NORMALS[idx1] : glm::vec3(0.0);
		return true;
	}

	// AABB外
	glm::vec3 x(0.0);
	for(int i = 0; i < 3; ++i){
		if(bout & (1 << (2*i))){
			x[i] = d0[2*i];
		}
		else if(bout & (1 << (2*i+1))){
			x[i] = -d0[2*i+1];
		}
	}

	// sgn = 1:箱，-1:オブジェクト
	if(c == 1){
		// 平面近傍
		d = sgn*d0[idx0];
		n = float(sgn)*RX_AABB_NORMALS[idx0];
	}
	else{
		// エッジ/コーナー近傍
		d = -sgn*glm::length(x);
		n = float(sgn)*(-glm::normalize(x));
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
inline bool GetImplicitPlane(const glm::vec3& pos, float& d, glm::vec3& n, glm::vec3& v, const glm::vec3& pn, const glm::vec3& pq)
{
	d = glm::dot(pq-pos, pn);
	n = pn;
	v = glm::vec3(0.0);
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
static inline void GetImplicitOpenBox(const glm::vec3& spos, const float& r, const int& sgn,
	const glm::vec3& sl0, const glm::vec3& sl1, float& d, glm::vec3& n)
{
	int t = 2;
	d = RX_FEQ_INF;

	float td;
	glm::vec3 tn;

	// 底
	glm::vec3 m0, m1;
	m0 = -sl1;
	m1 = sl1;
	m1[t] = sl0[t];
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 -x
	m0 = glm::vec3(-sl1[0], -sl1[1], -sl0[2]);
	m1 = glm::vec3(-sl0[0], sl1[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 +x
	m0 = glm::vec3(sl0[0], -sl1[1], -sl0[2]);
	m1 = glm::vec3(sl1[0], sl1[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 -y
	m0 = glm::vec3(-sl0[0], -sl1[1], -sl0[2]);
	m1 = glm::vec3(sl0[0], -sl0[1], sl1[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 +y
	m0 = glm::vec3(-sl0[0], sl0[1], -sl0[2]);
	m1 = glm::vec3(sl0[0], sl1[1], sl1[2]);
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
inline int ray_sphere(const glm::vec3& p, const glm::vec3& d, const glm::vec3& sc, const float r, float& t1, float& t2)
{
	glm::vec3 q = p-sc;	// 球中心座標系での光線原点座標

	float a = glm::length(d);
	float b = 2*glm::dot(q, d);
	float c = glm::length2(q)-r*r;

	// 判別式
	float D = b*b-4*a*c;

	if(D < 0.0){ // 交差なし
		return 0;
	}
	else if(D < RX_FEQ_EPS){ // 交点数1
		t1 = -b/(2*a);
		t2 = -1;
		return 1;
	}
	else{ // 交点数2
		float sqrtD = sqrt(D);
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
inline bool triangle_sphere(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n,
	const glm::vec3& c, const float& r, float& dist, glm::vec3& ipoint)
{
	// ポリゴンを含む平面と球中心の距離
	float d = glm::dot(v0, n);
	float l = glm::dot(n, c)-d;

	dist = l;
	if(l > r) return false;

	// 平面との最近傍点座標
	glm::vec3 p = c-l*n;

	// 近傍点が三角形内かどうかの判定
	glm::vec3 n1 = glm::cross((v0-c), (v1-c));
	glm::vec3 n2 = glm::cross((v1-c), (v2-c));
	glm::vec3 n3 = glm::cross((v2-c), (v0-c));

	ipoint = p;
	dist = l;
	if(glm::dot(n1, n2) > 0 && glm::dot(n2, n3) > 0){		// 三角形内
		return true;
	}
	else{		// 三角形外
	 // 三角形の各エッジと球の衝突判定
		for(int e = 0; e < 3; ++e){
			glm::vec3 va0 = (e == 0 ? v0 : (e == 1 ? v1 : v2));
			glm::vec3 va1 = (e == 0 ? v1 : (e == 1 ? v2 : v0));

			float t1, t2;
			int n = ray_sphere(va0, glm::normalize(va1-va0), c, r, t1, t2);

			if(n){
				float le2 = glm::length2(va1-va0);
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
inline float segment_point_dist(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& p)
{
	glm::vec3 v = glm::normalize(v1-v0);
	glm::vec3 vp = p-v0;
	glm::vec3 vh = glm::dot(vp, v)*v;
	return glm::length(vp-vh);
}

/*!
 * 線分と球の交差判定
 * @param[in] s0,s1	線分の端点
 * @param[in] sc,r   球の中心座標と半径
 * @param[out] d2 線分との距離の二乗
 * @return 交差ありでtrue
 */
inline bool segment_sphere(const glm::vec3& s0, const glm::vec3& s1, const glm::vec3& sc, const float& r, float& d2)
{
	glm::vec3 v = s1-s0;
	glm::vec3 c = sc-s0;

	float vc = glm::dot(v, c);
	if(vc < 0){		// 球の中心が線分の始点s0の外にある
		d2 = glm::length2(c);
		return (d2 < r* r);	// 球中心と始点s0の距離で交差判定
	}
	else{
		float v2 = glm::length2(v);
		if(vc > v2){	// 球の中心が線分の終点s1の外にある
			d2 = glm::length2(s1-sc);
			return (d2 < r* r);	// 球中心と終点s1の距離で交差判定
		}
		else{			// 球がs0とs1の間にある
			d2 = glm::length2((vc*v)/glm::length2(v)-c);
			return (d2 < r* r);	// 直線と球中心の距離で交差判定
		}
	}

	return false;
}





#endif // #ifndef _RX_SPH_COMMON_H_