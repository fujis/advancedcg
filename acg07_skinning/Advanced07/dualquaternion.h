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
#include <memory.h>

// glm
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"


//-----------------------------------------------------------------------------
//! クオータニオンクラス
//-----------------------------------------------------------------------------
class rxDualQuaternion
{
public:
	rxQuaternion m_real;	// Real part (回転を表す四元数)
	rxQuaternion m_dual;	// Dual part (平行移動を表す四元数=スカラー部が0+ベクトル部が平行移動ベクトル)

public:
	//! デフォルトコンストラクタ
	rxDualQuaternion()
	{
		m_real = rxQuaternion(0, 0, 0, 1);	// 回転なし(θ=0)で四元数を初期化
		m_dual = rxQuaternion(0, 0, 0, 0);	// 平行移動なし(0,0,0)で初期化
	}

	//! コンストラクタ:四元数2つで初期化
	rxDualQuaternion(rxQuaternion r, rxQuaternion d)
	{
		m_real = r; m_dual = d;
	}

	//! コンストラクタ:回転を表す四元数と平行移動ベクトルで初期化
	rxDualQuaternion(rxQuaternion r, Vec3 t)
	{
		m_real = r; 
		rxQuaternion tq = rxQuaternion(t[0], t[1], t[2], 0.0);
		m_dual = 0.5*tq*r;
	}

	//! 正規化
	void normalize()
	{
		real r2 = norm(m_real);
		if(fabs(r2) < 1.0e-6){
			return;
		}
		m_real *= 1.0/r2;
		m_dual *= 1.0/r2;
	}

	//! 共役
	rxDualQuaternion& conjugate()
	{
		m_real.conjugate();
		m_dual.conjugate();
		return *this;
	}

	//! 3次元座標のDQによる変換
	Vec3 transform(const Vec3 &p) const
	{
		double norm = m_real.norm();
		rxQuaternion q_real = m_real/norm;
		rxQuaternion q_dual = m_dual/norm;

		Vec3 v_real = q_real.GetVector();
		Vec3 v_dual = q_dual.GetVector();
		Vec3 trans = (v_dual*q_real.w_() - v_real*q_dual.w_() + cross(v_real, v_dual))*2.0;

		return q_real.rotate(p)+trans;
	}

	Vec3 rotate(const Vec3 &p) const
	{
		rxQuaternion tmp = m_real;
		tmp.normalize();
		return tmp.rotate(p);
	}

	//! 4x4変換行列からDQへの変換
	void setMat(const Mat4x4 &m)
	{
		rxQuaternion q_real(m);
		Vec3 v_trans(m(0,3), m(1,3), m(2,3));
		rxDualQuaternion dq(q_real, v_trans);
		*this = dq;
	}


	//! DQから4x4変換行列への変換
	Mat4x4 getMat(void) const
	{
		Mat4x4 m;

		return m;
	}



	//! 回転四元数の取得
	rxQuaternion getRotation()
	{
		return m_real;
	}

	//! 平行移動ベクトルの取得
	Vec3 getTranslation()
	{
		rxQuaternion rc = m_real;
		rc.conjugate();
		rxQuaternion t = (m_dual*2.0)*rc;
		return t.GetVector();
	}

	//
	// オペレータ
	// 
	rxDualQuaternion& operator+=(const rxDualQuaternion &dq)
	{
		m_real += dq.m_real;
		m_dual += dq.m_dual;
		return *this;
	}
	rxDualQuaternion& operator-=(const rxDualQuaternion &dq)
	{
		m_real -= dq.m_real;
		m_dual -= dq.m_dual;
		return *this;
	}
	rxDualQuaternion& operator*=(const rxDualQuaternion &dq)
	{
		rxQuaternion r, d;
		r = m_real*dq.m_real;
		d = m_dual*dq.m_real + m_real*dq.m_dual;
		m_real = r; m_dual = d;
		return *this;
	}
	rxDualQuaternion& operator/=(const rxDualQuaternion &dq)
	{
		real r2 = norm(dq.m_real);
		if(fabs(r2) > 1.0e-6){
			rxQuaternion r, d;
			r = m_real*dq.m_real/r2;
			d = (m_dual*dq.m_real - m_real*dq.m_dual)/r2;
			m_real = r; m_dual = d;
		}
		return *this;
	}
};

// 算術演算子はC++だとクラスのメンバ関数にするとa+10には対応できるが10+aに対応できないのでグローバル関数として定義する
inline rxDualQuaternion operator+(const rxDualQuaternion &dq1, const rxDualQuaternion &dq2)
{	
	rxDualQuaternion dq = dq1; 
	dq += dq2; 
	return dq; 
}
inline rxDualQuaternion operator-(const rxDualQuaternion &dq1, const rxDualQuaternion &dq2)
{	
	rxDualQuaternion dq = dq1; 
	dq -= dq2; 
	return dq; 
}
inline rxDualQuaternion operator*(const rxDualQuaternion &dq1, const rxDualQuaternion &dq2)
{	
	rxDualQuaternion dq = dq1; 
	dq *= dq2; 
	return dq; 
}
inline rxDualQuaternion operator/(const rxDualQuaternion &dq1, const rxDualQuaternion &dq2)
{	
	rxDualQuaternion dq = dq1; 
	dq /= dq2; 
	return dq; 
}
inline rxDualQuaternion operator*(const rxDualQuaternion &dq1, const real &s)
{	
	rxDualQuaternion dq = dq1; 
	dq.m_real *= s; 
	dq.m_dual *= s; 
	return dq; 
}
inline rxDualQuaternion operator*(const real &s, const rxDualQuaternion &dq1)
{	
	rxDualQuaternion dq = dq1; 
	dq.m_real *= s; 
	dq.m_dual *= s; 
	return dq; 
}




#endif// _DUAL_QUATERNION_H_
