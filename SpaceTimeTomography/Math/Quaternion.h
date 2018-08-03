#ifndef _Quaternion_h
#define _Quaternion_h

#include "Matrix4x4.h"

class  Quaternion
{
public:

	Quaternion();	
	Quaternion(float ax, float ay, float az, float angle);	
	Quaternion&	operator=(const Quaternion &q)	
	{
		x = q.x;
		y = q.y;
		z = q.z;
		w = q.w;
		return(*this);	
	}
	Quaternion(const Quaternion &q)                    
	{
		x = q.x;
		y = q.y;
		z = q.z;
		w = q.w;
	}

	void Identity()
	{
		x	= 0.0f;
		y	= 0.0f;
		z	= 0.0f;
		w	= 1.0f;
	}
	void Normalize()
	{
		float normInv = 1.0f/sqrtf(x*x+y*y+z*z+w*w);
		x *= normInv;
		y *= normInv;
		z *= normInv;
		w *= normInv;	
		return;
	}

	Quaternion& operator *= (const Quaternion&);	
	Matrix4x4 ToMatrix () const
	{
		Matrix4x4 mat;
		float *m = mat;
		m[0]  = 1.0f - 2.0f * (y*y + z*z);
		m[1]  = 2.0f * (x*y + z*w);
		m[2]  = 2.0f * (z*x - y*w);
		m[3]  = 0.0f;

		m[4]  = 2.0f * (x*y - z*w);
		m[5]  = 1.0f - 2.0f * (z*z + x*x);
		m[6]  = 2.0f * (y*z + x*w);
		m[7]  = 0.0f;

		m[8]  = 2.0f * (z*x + y*w);
		m[9]  = 2.0f * (y*z - x*w);
		m[10] = 1.0f - 2.0f * (x*x + y*y);
		m[11] = 0.0f;

		m[12] = 0.0f;
		m[13] = 0.0f;
		m[14] = 0.0f;
		m[15] = 1.0f;	
		return mat;
	}

	float	x,y,z,w;			

};

inline Quaternion operator * (const Quaternion &p, const Quaternion &q)
{
	Quaternion res;

	res.w	= p.w*q.w - p.x*q.x - p.y*q.y - p.z*q.z;
	res.x	= p.w*q.x + p.x*q.w + p.y*q.z - p.z*q.y;
	res.y	= p.w*q.y + p.y*q.w + p.z*q.x - p.x*q.z;
	res.z	= p.w*q.z + p.z*q.w + p.x*q.y - p.y*q.x;

	return res;
}
inline Quaternion& Quaternion::operator *= (const Quaternion& q)
{
	*this=*this*q;
	return (*this);	
}


#endif

