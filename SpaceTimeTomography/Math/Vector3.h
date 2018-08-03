#ifndef _Vector3_h
#define _Vector3_h

#include <math.h>

class  Vector3
{
public:
	float ele[3];
	Vector3();	
	Vector3(const Vector3 &v);
	Vector3(float x, float y, float z);
	Vector3(float srcVector[3]);
	operator float* () { return ele; }	
	operator const float* () const { return ele; }

	Vector3& operator = (const Vector3 &a)
	{
		ele[0] = a.ele[0];
		ele[1] = a.ele[1];
		ele[2] = a.ele[2];
		return *this; 
	}
	Vector3& operator *= (float);
	Vector3& operator += (const Vector3 &);
	Vector3& operator -= (const Vector3 &);

	void Normalize();
	float Length();
};

inline float operator * (const Vector3& A, const Vector3& B) 
{
	return A.ele[0]*B.ele[0] + A.ele[1]*B.ele[1] + A.ele[2]*B.ele[2];
}
//----------------------------------------------------------------------------
inline Vector3 operator % (const Vector3& A, const Vector3& B) 
{
	Vector3 Res;
	Res.ele[0] = A.ele[1] * B.ele[2] - A.ele[2] * B.ele[1]; 
	Res.ele[1] = A.ele[2] * B.ele[0] - A.ele[0] * B.ele[2]; 
	Res.ele[2] = A.ele[0] * B.ele[1] - A.ele[1] * B.ele[0]; 
	return Res;
}
//----------------------------------------------------------------------------
inline Vector3 operator * (const Vector3 &V, float s) 
{
	Vector3 Res;
	Res.ele[0] = V.ele[0] * s;
	Res.ele[1] = V.ele[1] * s;
	Res.ele[2] = V.ele[2] * s;
	return Res;
}
//----------------------------------------------------------------------------
inline Vector3 operator * (float s, const Vector3 &V) 
{
	Vector3 Res;
	Res.ele[0] = V.ele[0] * s;
	Res.ele[1] = V.ele[1] * s;
	Res.ele[2] = V.ele[2] * s;
	return Res;
}
//----------------------------------------------------------------------------
inline Vector3& Vector3::operator *= (float s) 
{
	ele[0] *= s;
	ele[1] *= s;
	ele[2] *= s;
	return *this;
}
//----------------------------------------------------------------------------
inline Vector3 operator + (const Vector3& A, const Vector3& B) 
{
	Vector3 Res;
	Res.ele[0] = A.ele[0] + B.ele[0];
	Res.ele[1] = A.ele[1] + B.ele[1];
	Res.ele[2] = A.ele[2] + B.ele[2];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector3 operator - (const Vector3& A, const Vector3& B) 
{
	Vector3 Res;
	Res.ele[0] = A.ele[0] - B.ele[0];
	Res.ele[1] = A.ele[1] - B.ele[1];
	Res.ele[2] = A.ele[2] - B.ele[2];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector3 operator - (const Vector3& A)
{
	Vector3 Res;
	Res.ele[0] = - A.ele[0];
	Res.ele[1] = - A.ele[1];
	Res.ele[2] = - A.ele[2];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector3 & Vector3::operator += (const Vector3 &B) 
{
	ele[0] += B.ele[0];
	ele[1] += B.ele[1];
	ele[2] += B.ele[2];
	return *this;
}
//----------------------------------------------------------------------------
inline Vector3 & Vector3::operator -= (const Vector3 &B) 
{
	ele[0] -= B.ele[0];
	ele[1] -= B.ele[1];
	ele[2] -= B.ele[2];
	return *this;
}
//----------------------------------------------------------------------------
inline void Vector3::Normalize() 
{
	float vecLenInv = 1.0f / sqrtf(ele[0]*ele[0] + ele[1]*ele[1] + ele[2]*ele[2]);
	ele[0] *= vecLenInv;
	ele[1] *= vecLenInv;
	ele[2] *= vecLenInv;
}
//----------------------------------------------------------------------------
inline float Vector3::Length()
{
	return sqrtf(ele[0]*ele[0] + ele[1]*ele[1] + ele[2]*ele[2]);
}


#endif
