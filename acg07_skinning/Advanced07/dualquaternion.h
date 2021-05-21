/*! @file rx_dualquaternion.h
	
	@brief デュアルクオータニオン
			- Dual Number[Clifford1882]を四元数に拡張したもの
			- キャラクターモーションでの剛体移動変換を表すのに用いる
			 
	@author Makoto Fujisawa
	@date  
*/

#ifndef _DUALQUATERNION_H_
#define _DUALQUATERNION_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// glm
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

//-----------------------------------------------------------------------------
//! Dual Quaternionクラス
//  - 回転を表す四元数(Real)と平行移動を表す四元数(Dual)の組み合わせ
//  - DQS(Dual Quaternion Skinning)用
//  - 四元数はglm::quatを使用(参考: https://glm.g-truc.net/0.9.0/api/a00135.html )
//-----------------------------------------------------------------------------
class DualQuaternion
{
public:
	glm::quat m_real;	// Real part (回転を表す四元数)
	glm::quat m_dual;	// Dual part (平行移動を表す四元数=スカラー部が0+ベクトル部が平行移動ベクトル)

public:
	//! デフォルトコンストラクタ
	DualQuaternion()
	{
		m_real = glm::quat(1, 0, 0, 0);	// 回転なし(θ=0)で四元数を初期化(glm::quatはw,x,y,zの順番)
		m_dual = glm::quat(0, 0, 0, 0);	// 平行移動なし(0,0,0)で初期化
	}

	//! コンストラクタ:四元数2つで初期化
	DualQuaternion(glm::quat r, glm::quat d)
	{
		m_real = r; m_dual = d;
	}

	//! コンストラクタ:回転を表す四元数と平行移動ベクトルで初期化
	DualQuaternion(glm::quat r, glm::vec3 t)
	{
		m_real = r; 
		glm::quat tq = glm::quat(0.0, t[0], t[1], t[2]);
		m_dual = 0.5*tq*r;
	}

	//! 正規化
	void normalize()
	{
		float r2 = glm::length(m_real);
		if(fabs(r2) < 1.0e-6){
			return;
		}
		m_real *= 1.0/r2;
		m_dual *= 1.0/r2;
	}

	//! 共役
	DualQuaternion& conjugate()
	{
		m_real = glm::conjugate(m_real);
		m_dual = glm::conjugate(m_dual);
		return *this;
	}

	//! 3次元座標のDQによる変換
	glm::vec3 transform(const glm::vec3 &p) const
	{
		glm::quat q_trans = 2.0f*m_dual*glm::conjugate(m_real);
		glm::vec3 trans(q_trans.x, q_trans.y, q_trans.z);
		// glmでの四元数(glm::quat)とベクトル(glm::vec3)の掛け算(*)はそのベクトルを四元数で回転させる(qvq*)
		return m_real*p+trans;
	}

	glm::vec3 rotate(const glm::vec3 &p) const
	{
		glm::quat tmp = m_real;
		tmp = glm::normalize(tmp);
		return tmp*p;
	}

	//! 4x4変換行列からDQへの変換
	void setMat(const glm::mat4 &m)
	{
		glm::quat q_real(m);		// mat4から回転を表す四元数を取り出す
		glm::vec3 v_trans(m[3]);	// mat4から平行移動ベクトルを取り出す
		DualQuaternion dq(q_real, v_trans);
		*this = dq;
	}

	//! 回転四元数の取得
	glm::quat getRotation()
	{
		return m_real;
	}

	//! 平行移動ベクトルの取得
	glm::vec3 getTranslation()
	{
		glm::quat rc = m_real;
		glm::quat t = 2.0f*m_dual*glm::conjugate(rc);
		return glm::vec3(t.x, t.y, t.z);
	}

	//
	// オペレータ
	// 
	DualQuaternion& operator+=(const DualQuaternion &dq)
	{
		m_real = m_real + dq.m_real;
		m_dual = m_dual + dq.m_dual;
		return *this;
	}
	DualQuaternion& operator-=(const DualQuaternion &dq)
	{
		m_real = m_real + (-dq.m_real);
		m_dual = m_dual + (-dq.m_dual);
		return *this;
	}
	DualQuaternion& operator*=(const DualQuaternion &dq)
	{
		glm::quat r, d;
		r = m_real*dq.m_real;
		d = m_dual*dq.m_real + m_real*dq.m_dual;
		m_real = r; m_dual = d;
		return *this;
	}
	DualQuaternion& operator/=(const DualQuaternion &dq)
	{
		float r2 = glm::length(dq.m_real);
		if(fabs(r2) > 1.0e-6){
			glm::quat r, d;
			r = m_real*dq.m_real/r2;
			d = (m_dual*dq.m_real + (-m_real*dq.m_dual))/r2;
			m_real = r; m_dual = d;
		}
		return *this;
	}
};

// 算術演算子はC++だとクラスのメンバ関数にするとa+10には対応できるが10+aに対応できないのでグローバル関数として定義する
inline DualQuaternion operator+(const DualQuaternion &dq1, const DualQuaternion &dq2)
{	
	DualQuaternion dq = dq1; 
	dq += dq2; 
	return dq; 
}
inline DualQuaternion operator-(const DualQuaternion &dq1, const DualQuaternion &dq2)
{	
	DualQuaternion dq = dq1; 
	dq -= dq2; 
	return dq; 
}
inline DualQuaternion operator*(const DualQuaternion &dq1, const DualQuaternion &dq2)
{	
	DualQuaternion dq = dq1; 
	dq *= dq2; 
	return dq; 
}
inline DualQuaternion operator/(const DualQuaternion &dq1, const DualQuaternion &dq2)
{	
	DualQuaternion dq = dq1; 
	dq /= dq2; 
	return dq; 
}
inline DualQuaternion operator*(const DualQuaternion &dq1, const float &s)
{	
	DualQuaternion dq = dq1; 
	dq.m_real = s*dq.m_real;
	dq.m_dual = s*dq.m_dual;
	return dq; 
}
inline DualQuaternion operator*(const float &s, const DualQuaternion &dq1)
{	
	DualQuaternion dq = dq1; 
	dq.m_real = s*dq.m_real;
	dq.m_dual = s*dq.m_dual;
	return dq; 
}




#endif// _DUAL_QUATERNION_H_
