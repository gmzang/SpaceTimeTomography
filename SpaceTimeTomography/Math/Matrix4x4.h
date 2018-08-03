#ifndef _Matrix4x4_h
#define _Matrix4x4_h

#include "Vector4.h"

class  Matrix4x4
{
public:
	float ele[16];

	Matrix4x4();
	Matrix4x4(const Matrix4x4 &m);
	Matrix4x4(float e0, float e1, float e2, float e3,
		float e4, float e5, float e6, float e7,
		float e8, float e9, float e10, float e11,
		float e12, float e13, float e14, float e15);

	Matrix4x4& operator = (const Matrix4x4 &m) 
	{
		ele[0] =  m.ele[0]; 
		ele[1] =  m.ele[1]; 
		ele[2] =  m.ele[2]; 
		ele[3] =  m.ele[3]; 
		ele[4] =  m.ele[4]; 
		ele[5] =  m.ele[5]; 
		ele[6] =  m.ele[6]; 
		ele[7] =  m.ele[7]; 
		ele[8] =  m.ele[8]; 
		ele[9] =  m.ele[9]; 
		ele[10] = m.ele[10];
		ele[11] = m.ele[11];
		ele[12] = m.ele[12];
		ele[13] = m.ele[13];
		ele[14] = m.ele[14];
		ele[15] = m.ele[15];
		return *this;
	}

	Matrix4x4& operator *= (const Matrix4x4 &);
	Matrix4x4& operator *= (float);
	Matrix4x4& operator += (const Matrix4x4 &);
	Matrix4x4& operator -= (const Matrix4x4 &);

	operator const float* () const { return ele; }	
	operator float* () { return ele; }

	void IdentityMatrix();

	void Transpose();	
	float Inverse();	

	static Matrix4x4 TranslateMatrix(float tx,float ty,float tz);
	static Matrix4x4 ScaleMatrix(float sx,float sy,float sz);
	static Matrix4x4 RotateMatrix(float angle,float x,float y,float z);

};

inline void Matrix4x4::IdentityMatrix() 
{
	ele[0] = ele[5] = ele[10] = ele[15] = 1.0f;
	ele[1] = ele[2] = ele[3] = ele[4] = 0.0f;
	ele[6] = ele[7] = ele[8] = ele[9] = 0.0f;
	ele[11] = ele[12] = ele[13] = ele[14] = 0.0f;
}
//----------------------------------------------------------------------------
inline Matrix4x4 operator * (const Matrix4x4 &A, const Matrix4x4 &B) 
{
	Matrix4x4 Res;

	Res.ele[0]   = A.ele[0]*B.ele[0] + A.ele[4]*B.ele[1] + A.ele[8]*B.ele[2] + A.ele[12]*B.ele[3]; 
	Res.ele[1]   = A.ele[1]*B.ele[0] + A.ele[5]*B.ele[1] + A.ele[9]*B.ele[2] + A.ele[13]*B.ele[3]; 
	Res.ele[2]   = A.ele[2]*B.ele[0] + A.ele[6]*B.ele[1] + A.ele[10]*B.ele[2] + A.ele[14]*B.ele[3]; 
	Res.ele[3]   = A.ele[3]*B.ele[0] + A.ele[7]*B.ele[1] + A.ele[11]*B.ele[2] + A.ele[15]*B.ele[3]; 

	Res.ele[4]   = A.ele[0]*B.ele[4] + A.ele[4]*B.ele[5] + A.ele[8]*B.ele[6] + A.ele[12]*B.ele[7]; 
	Res.ele[5]   = A.ele[1]*B.ele[4] + A.ele[5]*B.ele[5] + A.ele[9]*B.ele[6] + A.ele[13]*B.ele[7]; 
	Res.ele[6]   = A.ele[2]*B.ele[4] + A.ele[6]*B.ele[5] + A.ele[10]*B.ele[6] + A.ele[14]*B.ele[7]; 
	Res.ele[7]   = A.ele[3]*B.ele[4] + A.ele[7]*B.ele[5] + A.ele[11]*B.ele[6] + A.ele[15]*B.ele[7]; 

	Res.ele[8]   = A.ele[0]*B.ele[8] + A.ele[4]*B.ele[9] + A.ele[8]*B.ele[10] + A.ele[12]*B.ele[11]; 
	Res.ele[9]   = A.ele[1]*B.ele[8] + A.ele[5]*B.ele[9] + A.ele[9]*B.ele[10] + A.ele[13]*B.ele[11]; 
	Res.ele[10]   = A.ele[2]*B.ele[8] + A.ele[6]*B.ele[9] + A.ele[10]*B.ele[10] + A.ele[14]*B.ele[11]; 
	Res.ele[11]   = A.ele[3]*B.ele[8] + A.ele[7]*B.ele[9] + A.ele[11]*B.ele[10] + A.ele[15]*B.ele[11]; 

	Res.ele[12]   = A.ele[0]*B.ele[12] + A.ele[4]*B.ele[13] + A.ele[8]*B.ele[14] + A.ele[12]*B.ele[15]; 	
	Res.ele[13]   = A.ele[1]*B.ele[12] + A.ele[5]*B.ele[13] + A.ele[9]*B.ele[14] + A.ele[13]*B.ele[15]; 	
	Res.ele[14]   = A.ele[2]*B.ele[12] + A.ele[6]*B.ele[13] + A.ele[10]*B.ele[14] + A.ele[14]*B.ele[15]; 	
	Res.ele[15]   = A.ele[3]*B.ele[12] + A.ele[7]*B.ele[13] + A.ele[11]*B.ele[14] + A.ele[15]*B.ele[15]; 		
	return Res;	

}
//----------------------------------------------------------------------------
inline Matrix4x4& Matrix4x4::operator *= (const Matrix4x4 &B) 
{
	*this=(*this)*B;
	return *this;
}
//----------------------------------------------------------------------------
inline Matrix4x4 operator * (const Matrix4x4 &A, float s) 
{
	Matrix4x4 Res;
	Res.ele[0] = A.ele[0] * s;
	Res.ele[1] = A.ele[1] * s;
	Res.ele[2] = A.ele[2] * s;
	Res.ele[3] = A.ele[3] * s;
	Res.ele[4] = A.ele[4] * s;
	Res.ele[5] = A.ele[5] * s;
	Res.ele[6] = A.ele[6] * s;
	Res.ele[7] = A.ele[7] * s;
	Res.ele[8] = A.ele[8] * s;
	Res.ele[9] = A.ele[9] * s;
	Res.ele[10] = A.ele[10] * s;
	Res.ele[11] = A.ele[11] * s;
	Res.ele[12] = A.ele[12] * s;
	Res.ele[13] = A.ele[13] * s;
	Res.ele[14] = A.ele[14] * s;
	Res.ele[15] = A.ele[15] * s;
	return Res;	
}
//----------------------------------------------------------------------------
inline Matrix4x4 operator * (float s, const Matrix4x4 &A) 
{
	Matrix4x4 Res;
	Res.ele[0] = A.ele[0] * s;
	Res.ele[1] = A.ele[1] * s;
	Res.ele[2] = A.ele[2] * s;
	Res.ele[3] = A.ele[3] * s;
	Res.ele[4] = A.ele[4] * s;
	Res.ele[5] = A.ele[5] * s;
	Res.ele[6] = A.ele[6] * s;
	Res.ele[7] = A.ele[7] * s;
	Res.ele[8] = A.ele[8] * s;
	Res.ele[9] = A.ele[9] * s;
	Res.ele[10] = A.ele[10] * s;
	Res.ele[11] = A.ele[11] * s;
	Res.ele[12] = A.ele[12] * s;
	Res.ele[13] = A.ele[13] * s;
	Res.ele[14] = A.ele[14] * s;
	Res.ele[15] = A.ele[15] * s;
	return Res;	

}
//----------------------------------------------------------------------------
inline Matrix4x4& Matrix4x4::operator *= (float s) 
{
	ele[0] *= s;
	ele[1] *= s;
	ele[2] *= s;
	ele[3] *= s;
	ele[4] *= s;
	ele[5] *= s;
	ele[6] *= s;
	ele[7] *= s;
	ele[8] *= s;
	ele[9] *= s;
	ele[10] *= s;
	ele[11] *= s;
	ele[12] *= s;
	ele[13] *= s;
	ele[14] *= s;
	ele[15] *= s;
	return *this;
}
//----------------------------------------------------------------------------
inline Matrix4x4 operator + (const Matrix4x4 &A, const Matrix4x4 &B) 
{
	Matrix4x4 Res;

	Res.ele[0] = A.ele[0] + B.ele[0];
	Res.ele[1] = A.ele[1] + B.ele[1];
	Res.ele[2] = A.ele[2] + B.ele[2];
	Res.ele[3] = A.ele[3] + B.ele[3];
	Res.ele[4] = A.ele[4] + B.ele[4];
	Res.ele[5] = A.ele[5] + B.ele[5];
	Res.ele[6] = A.ele[6] + B.ele[6];
	Res.ele[7] = A.ele[7] + B.ele[7];	
	Res.ele[8] = A.ele[8] + B.ele[8];
	Res.ele[9] = A.ele[9] + B.ele[9];
	Res.ele[10] = A.ele[10] + B.ele[10];
	Res.ele[11] = A.ele[11] + B.ele[11];
	Res.ele[12] = A.ele[12] + B.ele[12];
	Res.ele[13] = A.ele[13] + B.ele[13];
	Res.ele[14] = A.ele[14] + B.ele[14];
	Res.ele[15] = A.ele[15] + B.ele[15];
	return Res;

}
//----------------------------------------------------------------------------
inline Matrix4x4 & Matrix4x4::operator += (const Matrix4x4 &B) 
{
	ele[0] += B.ele[0];
	ele[1] += B.ele[1];
	ele[2] += B.ele[2];
	ele[3] += B.ele[3];
	ele[4] += B.ele[4];
	ele[5] += B.ele[5];
	ele[6] += B.ele[6];
	ele[7] += B.ele[7];
	ele[8] += B.ele[8];
	ele[9] += B.ele[9];
	ele[10] += B.ele[10];
	ele[11] += B.ele[11];
	ele[12] += B.ele[12];
	ele[13] += B.ele[13];
	ele[14] += B.ele[14];
	ele[15] += B.ele[15];
	return *this;
}
//----------------------------------------------------------------------------
inline Matrix4x4 operator - (const Matrix4x4 &A, const Matrix4x4 &B) 
{
	Matrix4x4 Res;

	Res.ele[0] = A.ele[0] - B.ele[0];
	Res.ele[1] = A.ele[1] - B.ele[1];
	Res.ele[2] = A.ele[2] - B.ele[2];
	Res.ele[3] = A.ele[3] - B.ele[3];
	Res.ele[4] = A.ele[4] - B.ele[4];
	Res.ele[5] = A.ele[5] - B.ele[5];
	Res.ele[6] = A.ele[6] - B.ele[6];
	Res.ele[7] = A.ele[7] - B.ele[7];	
	Res.ele[8] = A.ele[8] - B.ele[8];
	Res.ele[9] = A.ele[9] - B.ele[9];
	Res.ele[10] = A.ele[10] - B.ele[10];
	Res.ele[11] = A.ele[11] - B.ele[11];
	Res.ele[12] = A.ele[12] - B.ele[12];
	Res.ele[13] = A.ele[13] - B.ele[13];
	Res.ele[14] = A.ele[14] - B.ele[14];
	Res.ele[15] = A.ele[15] - B.ele[15];
	return Res;

}
//----------------------------------------------------------------------------
inline Matrix4x4 operator - (const Matrix4x4 &A)
{
	Matrix4x4 Res;

	Res.ele[0] = - A.ele[0];
	Res.ele[1] = - A.ele[1];
	Res.ele[2] = - A.ele[2];
	Res.ele[3] = - A.ele[3];
	Res.ele[4] = - A.ele[4];
	Res.ele[5] = - A.ele[5];
	Res.ele[6] = - A.ele[6];
	Res.ele[7] = - A.ele[7];	
	Res.ele[8] = - A.ele[8];
	Res.ele[9] = - A.ele[9];
	Res.ele[10] = - A.ele[10];
	Res.ele[11] = - A.ele[11];
	Res.ele[12] = - A.ele[12];
	Res.ele[13] = - A.ele[13];
	Res.ele[14] = - A.ele[14];
	Res.ele[15] = - A.ele[15];
	return Res;

}
//----------------------------------------------------------------------------
inline Matrix4x4 & Matrix4x4::operator -= (const Matrix4x4 &B) 
{

	ele[0] -= B.ele[0];
	ele[1] -= B.ele[1];
	ele[2] -= B.ele[2];
	ele[3] -= B.ele[3];
	ele[4] -= B.ele[4];
	ele[5] -= B.ele[5];
	ele[6] -= B.ele[6];
	ele[7] -= B.ele[7];
	ele[8] -= B.ele[8];
	ele[9] -= B.ele[9];
	ele[10] -= B.ele[10];
	ele[11] -= B.ele[11];
	ele[12] -= B.ele[12];
	ele[13] -= B.ele[13];
	ele[14] -= B.ele[14];
	ele[15] -= B.ele[15];
	return *this;

}
//----------------------------------------------------------------------------
inline Vector4 operator * (const Vector4& Vec, const Matrix4x4& Mat) 
{
	Vector4 Res;


	Res.ele[0] = Mat.ele[0]*Vec.ele[0] + Mat.ele[1]*Vec.ele[1] + Mat.ele[2]*Vec.ele[2] + Mat.ele[3]*Vec.ele[3]; 
	Res.ele[1] = Mat.ele[4]*Vec.ele[0] + Mat.ele[5]*Vec.ele[1] + Mat.ele[6]*Vec.ele[2] + Mat.ele[7]*Vec.ele[3]; 
	Res.ele[2] = Mat.ele[8]*Vec.ele[0] + Mat.ele[9]*Vec.ele[1] + Mat.ele[10]*Vec.ele[2] + Mat.ele[11]*Vec.ele[3]; 
	Res.ele[3] = Mat.ele[12]*Vec.ele[0] + Mat.ele[13]*Vec.ele[1] + Mat.ele[14]*Vec.ele[2] + Mat.ele[15]*Vec.ele[3]; 
	return Res;	

}
//----------------------------------------------------------------------------
inline Vector4 operator * (const Matrix4x4& Mat, const Vector4& Vec) 
{
	Vector4 Res;

	Res.ele[0] = Mat.ele[0]*Vec.ele[0] + Mat.ele[4]*Vec.ele[1] + Mat.ele[8]*Vec.ele[2] + Mat.ele[12]*Vec.ele[3]; 
	Res.ele[1] = Mat.ele[1]*Vec.ele[0] + Mat.ele[5]*Vec.ele[1] + Mat.ele[9]*Vec.ele[2] + Mat.ele[13]*Vec.ele[3]; 
	Res.ele[2] = Mat.ele[2]*Vec.ele[0] + Mat.ele[6]*Vec.ele[1] + Mat.ele[10]*Vec.ele[2] + Mat.ele[14]*Vec.ele[3]; 
	Res.ele[3] = Mat.ele[3]*Vec.ele[0] + Mat.ele[7]*Vec.ele[1] + Mat.ele[11]*Vec.ele[2] + Mat.ele[15]*Vec.ele[3]; 
	return Res;	

}

#endif
