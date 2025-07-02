/*!
  @file rx_msh.h
	
  @brief MSH File Input
         Gmsh(http://gmsh.info//)の出力ファイル(*.msh)のインポータ
		 - 表面三角形ポリゴンと四面体メッシュの読み込み
		 - 表面が三角形で構成され，内部に四面体メッシュが定義されているMSHファイルのみを対象としている
		 - 四角形ポリゴン+六面体メッシュなどは対応してないので注意
 
  @author Makoto Fujisawa
  @date   2021
*/

#ifndef _RX_MSH_H_
#define _RX_MSH_H_


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------
#include <fstream>
#include <sstream>
#include "rx_mesh.h"

//-----------------------------------------------------------------------------
// Name Space
//-----------------------------------------------------------------------------
using namespace std;



//-----------------------------------------------------------------------------
// rxMSHクラスの宣言 - OBJ形式の読み込み
//-----------------------------------------------------------------------------
class rxMSH
{
public:
	//! コンストラクタ
	rxMSH(){}

	//! デストラクタ
	~rxMSH(){}

	/*!
	 * MSHファイル読み込み
	 * @param[in] file_name ファイル名(フルパス)
	 * @param[out] vrts 頂点座標
	 * @param[out] tets 四面体情報
	 */
	bool Read(string file_name, vector<glm::vec3> &vrts, vector< vector<int> > &tris, vector< vector<int> > &tets);

};

/*!
 * 最後の空白(スペース，タブ)を削除
 * @param[inout] buf 処理文字列
 */
inline void DeleteLastSpace(std::string &buf)
{
	if(buf.empty()) return;
	std::size_t pos;
	while((pos = buf.find_last_of(" 　\t")) == buf.size()-1){
		buf.erase(pos);
		if(buf.empty()) break;
	}
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


/*!
* OBJファイル読み込み
* @param[in] file_name ファイル名(フルパス)
* @param[out] vrts 頂点座標
* @param[out] vnms 頂点法線
* @param[out] poly ポリゴン
* @param[out] mats 材質情報
* @param[in] triangle ポリゴンの三角形分割フラグ
*/
inline bool rxMSH::Read(string file_name, 
						vector<glm::vec3> &vrts, 
						vector< vector<int> > &tris, 
						vector< vector<int> > &tets)
{
	ifstream file;

	file.open(file_name.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "rxMSH::Read : Invalid file specified" << endl;
		return false;
	}

	int div = 0;
	int line = 0;
	int num_nodes = 0;
	int num_elems = 0;

	vector<glm::vec3> nodes;

	int offset = static_cast<int>(vrts.size());
	int cnt = 0;
	int elem_id = 0;
	cout << "rxMSH::Read : " << file_name << endl;

	string buf;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);

		// '#'以降はコメントとして無視
		if((comment_start = buf.find('#')) != string::size_type(-1))
			buf = buf.substr(0, comment_start);

		// 行頭,行尾のスペース，タブを削除
		DeleteHeadSpace(buf);
		DeleteLastSpace(buf);

		// 空行は無視
		if(buf.empty())
			continue;


		// ノード座標(頂点座標)
		if(buf.find("$Nodes") != std::string::npos){
			div = 1;
			line = 0;
		}
		if(buf.find("$EndNodes") != std::string::npos){
			div = 0;
			for(int i = 0; i < num_nodes; ++i) vrts.push_back(nodes[i]);
			cout << "  num. nodes : " << nodes.size() << endl;

		}
		if(div == 1){
			if(line == 1){
				int tmp1, tmp2, tmp3;
				std::stringstream s(buf);
				s >> tmp1 >> tmp2 >> tmp3 >> num_nodes;
				if(num_nodes == 0) break;
			}
			else if(nodes.size() < num_nodes){
				int num_space = CountString(buf, 0, " ");
				if(num_space == 2){
					glm::vec3 p;
					if(StringToVec3s(buf, "", p) == 3){
						nodes.push_back(p);
					}
				}
			}

			line++;
		}

		// 四面体
		if(buf.find("$Elements") != std::string::npos){
			div = 2;
			line = 0;
			cnt = 0;
			elem_id = 0;
		}
		if(buf.find("$EndElements") != std::string::npos){
			div = 0;

			cout << "  num. elems : " << num_elems;
			cout << " (tris: " << tris.size() << ", tets: " << tets.size() << ")" << endl;
		}
		if(div == 2){
			if(line == 1){
				// エレメント数(ここでは三角形ポリゴンと四面体の数の合計を想定)
				int tmp1, tmp2, tmp3;
				if(sscanf(&buf[0], "%d %d %d %d", &tmp1, &tmp2, &tmp3, &num_elems) != 4){
					continue;
				}
				if(num_elems == 0) break;
			}
			else if(line >= 3){
				string sub;
				size_t pos;

				vector<int> t(3, 0);
				pos = GetNextString(buf, sub, 0, " ");
				int id = atoi(sub.c_str());
				if(id != elem_id+1){
					line++; continue;
				}
				elem_id++;

				pos = GetNextString(buf, sub, pos, " "); t[0] = atof(sub.c_str());
				pos = GetNextString(buf, sub, pos, " "); t[1] = atof(sub.c_str());
				pos = GetNextString(buf, sub, pos, " "); t[2] = atof(sub.c_str());
				//if(pos == string::npos){	// 三角形 (こちらだとmac環境で上手くいかなかった)
				if(pos >= buf.size()){	// 三角形
					for(int i = 0; i < 3; ++i) t[i] += offset-1;
					tris.push_back(t);
				}
				else{	// bufに更に文字列が含まれる場合(もう一つ要素がある場合は)は四面体として処理
					pos = GetNextString(buf, sub, pos, " "); t.push_back(atof(sub.c_str()));
					for(int i = 0; i < 4; ++i) t[i] += offset-1;
					if(t[3] < 0){ // 最後の要素が-1の場合(atofが0を返した場合)は三角形として登録
						t.pop_back();
						tris.push_back(t);
					}
					else{	// 四面体の登録
						tets.push_back(t);
					}
				}

				cnt++;
			}

			line++;
		}

	}

	file.close();

	return true;
}


#endif // _RX_MSH_H_
