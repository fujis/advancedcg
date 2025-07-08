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

	int offset = vrts.size();
	int cnt = 0;
	int elem_id = 0;

	string buf;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);

		// '#'以降はコメントとして無視
		if((comment_start = buf.find('#')) != string::size_type(-1))
			buf = buf.substr(0, comment_start);

		// 行頭のスペース，タブを削除
		DeleteHeadSpace(buf);

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

			//cout << "num. verts : " << num_nodes << ", " << nodes.size() << endl;
			//for(int i = 0; i < num_nodes; ++i) cout << glm::to_string(nodes[i]) << endl;
		}
		if(div == 1){
			if(line == 1){
				int tmp1, tmp2, tmp3;
				if(sscanf(&buf[0], "%d %d %d %d", &tmp1, &tmp2, &tmp3, &num_nodes) != 4){
					break;
				}
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

			cout << "num. elems : " << num_elems << endl;
			cout << "num. tris : " << tris.size() << endl;
			cout << "num. tets : " << tets.size() << endl;
			//for(int i = 0; i < num_nodes; ++i) cout << glm::to_string(nodes[i]) << endl;
		}
		if(div == 2){
			if(line == 1){
				// エレメント数(ここでは三角形ポリゴンと四面体の数の合計を想定)
				int tmp1, tmp2, tmp3;
				if(sscanf(&buf[0], "%d %d %d %d", &tmp1, &tmp2, &tmp3, &num_elems) != 4){
					break;
				}
				if(num_elems == 0) break;
			}
			else if(line >= 3){
				int id;
				if(sscanf(&buf[0], "%d", &id) != 1){
					break;
				}
				if(id != elem_id+1) continue;
				elem_id++;
				
				int num_space = CountString(buf, 0, " ");
				int tmp;
				if(num_space == 3){	// 三角形
					vector<int> t(3, 0);
					if(sscanf(&buf[0], "%d %d %d %d", &tmp, &t[0], &t[1], &t[2]) != 4){
						break;
					}
					for(int i = 0; i < 3; ++i) t[i] += offset-1;
					tris.push_back(t);
				}
				else if(num_space == 4){	// 四面体
					vector<int> t(4, 0);
					if(sscanf(&buf[0], "%d %d %d %d %d", &tmp, &t[0], &t[1], &t[2], &t[3]) != 5){
						break;
					}
					for(int i = 0; i < 4; ++i) t[i] += offset-1;
					tets.push_back(t);
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
