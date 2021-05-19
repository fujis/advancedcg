/*!
  @file pointsprite.fp
	
  @brief GLSLフラグメントシェーダ
		 - PointSpriteによる点描画
 
  @author Makoto Fujisawa
  @date 2011
*/
// FILE --pointsprite.fp--
#version 120


//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
varying float vValid;
uniform vec3 lightDir;	// 光源方向

//-----------------------------------------------------------------------------
// エントリ関数
//-----------------------------------------------------------------------------

void main(void)
{
	if(vValid < 0.5) discard;
	//const vec3 lightDir = vec3(0.577, 0.577, 0.577);

	// テクスチャ座標から法線を計算(球として描画)
	vec3 N;
	N.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0)+vec2(-1.0, 1.0);

	float mag = dot(N.xy, N.xy);	// 中心からの2乗距離
	if(mag > 1.0) discard;   // 円の外のピクセルは破棄
	N.z = sqrt(1.0-mag);

	// 光源方向と法線から表面色を計算
	float diffuse = max(0.0, dot(lightDir, N));

	gl_FragColor = gl_Color*diffuse;
}


