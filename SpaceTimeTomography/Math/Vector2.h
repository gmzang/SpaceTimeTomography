#ifndef _Vector2_h
#define _Vector2_h

#include <math.h>

class  Vector2
{
public:
	float ele[2];
	Vector2();	
	Vector2(const Vector2 &v);
	Vector2(float x, float y);
	Vector2(float srcVector[2]);
	operator float* () { return ele; }	
	operator const float* () const { return ele; }

	Vector2& operator = (const Vector2 &a)
	{
		ele[0] = a.ele[0];
		ele[1] = a.ele[1];
		return *this; 
	}
	Vector2& operator *= (float);
	Vector2& operator += (const Vector2 &);
	Vector2& operator -= (const Vector2 &);

	void Normalize();
	float Length();
};

inline float operator * (const Vector2& A, const Vector2& B) 
{
	return A.ele[0]*B.ele[0] + A.ele[1]*B.ele[1];
}
//----------------------------------------------------------------------------
inline Vector2 operator * (const Vector2 &V, float s) 
{
	Vector2 Res;
	Res.ele[0] = V.ele[0] * s;
	Res.ele[1] = V.ele[1] * s;
	return Res;
}
//----------------------------------------------------------------------------
inline Vector2 operator * (float s, const Vector2 &V) 
{
	Vector2 Res;
	Res.ele[0] = V.ele[0] * s;
	Res.ele[1] = V.ele[1] * s;
	return Res;
}
//----------------------------------------------------------------------------
inline Vector2& Vector2::operator *= (float s) 
{
	ele[0] *= s;
	ele[1] *= s;
	return *this;
}
//----------------------------------------------------------------------------
inline Vector2 operator + (const Vector2& A, const Vector2& B) 
{
	Vector2 Res;
	Res.ele[0] = A.ele[0] + B.ele[0];
	Res.ele[1] = A.ele[1] + B.ele[1];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector2 operator - (const Vector2& A, const Vector2& B) 
{
	Vector2 Res;
	Res.ele[0] = A.ele[0] - B.ele[0];
	Res.ele[1] = A.ele[1] - B.ele[1];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector2 operator - (const Vector2& A)
{
	Vector2 Res;
	Res.ele[0] = - A.ele[0];
	Res.ele[1] = - A.ele[1];
	return Res;
}
//----------------------------------------------------------------------------
inline Vector2 & Vector2::operator += (const Vector2 &B) 
{
	ele[0] += B.ele[0];
	ele[1] += B.ele[1];
	return *this;
}
//----------------------------------------------------------------------------
inline Vector2 & Vector2::operator -= (const Vector2 &B) 
{
	ele[0] -= B.ele[0];
	ele[1] -= B.ele[1];
	return *this;
}
//----------------------------------------------------------------------------
inline void Vector2::Normalize() 
{
	float vecLenInv = 1.0f / sqrtf(ele[0]*ele[0] + ele[1]*ele[1]);
	ele[0] *= vecLenInv;
	ele[1] *= vecLenInv;
}
//----------------------------------------------------------------------------
inline float Vector2::Length()
{
	return sqrtf(ele[0]*ele[0] + ele[1]*ele[1] );
}


#endif
