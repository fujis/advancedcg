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
	DualQuaternion conjugate() const
	{
		DualQuaternion dq_conj;
		dq_conj.m_real = glm::conjugate(m_real);
		dq_conj.m_dual = glm::conjugate(m_dual);
		return dq_conj;
	}

	//! 回転四元数の取得
	glm::quat getRotation()
	{
		return m_real;
	}

	//! 平行移動ベクトルの取得
	glm::quat getTranslation()
	{
		return m_dual;
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
