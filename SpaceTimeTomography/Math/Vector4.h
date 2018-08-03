#ifndef _Vector4_h
#define _Vector4_h

#include <math.h>

class  Vector4
{
public:
	float ele[4];
	Vector4();	
	Vector4(const Vector4 &v);
	Vector4(float x, float y, float z, float w = 1.0f);
	Vector4(float srcVector[4]);
	operator float* () { return ele; }	
	operator const float* () const { return ele; }

	Vector4& operator = (const Vector4 &a)
	{
		ele[0] = a.ele[0];
		ele[1] = a.ele[1];
		ele[2] = a.ele[2];
		ele[3] = a.ele[3];
		return *this; 
	}
	Vector4& operator *= (float);
	Vector4& operator += (const Vector4 &);
	Vector4& operator -= (const Vector4 &);

	void Normalize();
	float Length();
};

inline float operator * (const Vector4& A, const Vector4& B) 
{
	return A.ele[0]*B.ele[0] + A.ele[1]*B.ele[1] + A.ele[2]*B.ele[2];
}
//----------------------------------------------------------------------------
inline Vector4 operator % (const Vector4& A, const Vector4& B) 
{
	Vector4 Res;
	Res.ele[0] = A.ele[1] * B.ele[2] - A.ele[2] * B.ele[1]; 
	Res.ele[1] = A.ele[2] * B.ele[0] - A.ele[0] * B.ele[2]; 
	Res.ele[2] = A.ele[0] * B.ele[1] - A.ele[1] * B.ele[0]; 
	Res.ele[3] = 0.0f;
	return Res;
}
//----------------------------------------------------------------------------
inline Vector4 operator * (const Vector4 &V, float s) 
{
	Vector4 Res;
	Res.ele[0] = V.ele[0] * s;
	Res.ele[1] = V.ele[1] * s;
	Res.ele[2] = V.ele[2] * s;
	Res.ele[3] = V.ele[3] * s;	
	return Res;
}
//----------------------------------------------------------------------------
inline Vector4 operator * (float s, const Vector4 &V) 
{
	Vector4 Res;
	Res.ele[0] = V.ele[0] * s;
	Res.ele[1] = V.ele[1] * s;
	Res.ele[2] = V.ele[2] * s;
	Res.ele[3] = V.ele[3] * s;	
	return Res;
}
//----------------------------------------------------------------------------
inline Vector4& Vector4::operator *= (float s) 
{
	ele[0] *= s;
	ele[1] *= s;
	ele[2] *= s;
	ele[3] *= s;	
	return *this;
}
//----------------------------------------------------------------------------
inline Vector4 operator + (const Vector4& A, const Vector4& B) 
{
	Vector4 Res;
	Res.ele[0] = A.ele[0] + B.ele[0];
	Res.ele[1] = A.ele[1] + B.ele[1];
	Res.ele[2] = A.ele[2] + B.ele[2];
	Res.ele[3] = A.ele[3] + B.ele[3];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector4 operator - (const Vector4& A, const Vector4& B) 
{
	Vector4 Res;
	Res.ele[0] = A.ele[0] - B.ele[0];
	Res.ele[1] = A.ele[1] - B.ele[1];
	Res.ele[2] = A.ele[2] - B.ele[2];
	Res.ele[3] = A.ele[3] - B.ele[3];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector4 operator - (const Vector4& A)
{
	Vector4 Res;
	Res.ele[0] = - A.ele[0];
	Res.ele[1] = - A.ele[1];
	Res.ele[2] = - A.ele[2];
	Res.ele[3] = - A.ele[3];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector4 & Vector4::operator += (const Vector4 &B) 
{
	ele[0] += B.ele[0];
	ele[1] += B.ele[1];
	ele[2] += B.ele[2];
	ele[3] += B.ele[3];
	return *this;
}
//----------------------------------------------------------------------------
inline Vector4 & Vector4::operator -= (const Vector4 &B) 
{
	ele[0] -= B.ele[0];
	ele[1] -= B.ele[1];
	ele[2] -= B.ele[2];
	ele[3] -= B.ele[3];
	return *this;
}
//----------------------------------------------------------------------------
inline void Vector4::Normalize() 
{
	float vecLenInv = 1.0f / sqrtf(ele[0]*ele[0] + ele[1]*ele[1] + ele[2]*ele[2]);
	ele[0] *= vecLenInv;
	ele[1] *= vecLenInv;
	ele[2] *= vecLenInv;
}
//----------------------------------------------------------------------------
inline float Vector4::Length()
{
	return sqrtf(ele[0]*ele[0] + ele[1]*ele[1] + ele[2]*ele[2]);
}


#endif
