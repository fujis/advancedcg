/*!
  @file shading.fp

  @brief GLSLフラグメントシェーダ

  @author
  @date 2021
*/

// mac環境でテクスチャ座標を確実に読み込むためにlayoutの指定が可能なver4.1にしています．
// Mac mini(Late 2012)(mac os 10.15.7 Catalina)で動作チェック済み
#version 410 core

//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
// GLから設定される定数(uniform)
uniform vec3 color;	// 描画色
uniform sampler2D tex;	// 表面テクスチャ
uniform float alpha;	// テクスチャ色と描画色の混合割合

// 頂点シェーダからの入力
in vec2 vtexcoord;	// テクスチャ座標

// 出力
out vec4 fragcolor;

//-----------------------------------------------------------------------------
// エントリ関数
//-----------------------------------------------------------------------------
void main(void)
{
	vec3 texcolor = texture(tex, vtexcoord).rgb;
	fragcolor = vec4(alpha*texcolor+(1.0-alpha)*color, 1.0);
}

