/*!
  @file shading.fp
	
  @brief GLSLフラグメントシェーダ
 
  @author Makoto Fujisawa
  @date 2011
*/
#version 120


//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
// GLから設定される定数(uniform)
uniform vec3 v3LightPosEye;	// カメラ座標系での光源位置
uniform vec3 v3LightColor;	// 光源色
uniform sampler2D texFloor;	// 表面テクスチャ
//uniform sampler2D shadowTex;	// 影用テクスチャ

// 頂点シェーダから渡される変数(varying)
varying vec4 v4VertexPosEye;  // カメラ座標系での頂点位置
varying vec3 v3NormalEye;		// カメラ座標系での頂点法線

//-----------------------------------------------------------------------------
// 関数
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// エントリ関数
//-----------------------------------------------------------------------------
void main(void)
{	
	//vec4 shadowPos = gl_TextureMatrix[0] * vertexPosEye;
	//vec4 colorMap = vec4(0.9, 0.9, 0.9, 0.0);
	vec4 colorMap = texture2D(texFloor, gl_TexCoord[0].xy);

	vec3 N = normalize(v3NormalEye);
	vec3 L = normalize(v3LightPosEye - v4VertexPosEye.xyz);
	float diffuse = max(0.0, dot(N, L));

	vec3 shadow = vec3(1.0);
	//vec3 shadow = vec3(1.0) -texture2DProj(shadowTex, shadowPos.xyw).xyz;
	//if(shadowPos.w < 0.0) shadow = lightColor;	// avoid back projections

	gl_FragColor = vec4(gl_Color.xyz *colorMap.xyz *diffuse * shadow, 1.0);
	//gl_FragColor = vec4(gl_Color.xyz*diffuse, 1.0);
}

