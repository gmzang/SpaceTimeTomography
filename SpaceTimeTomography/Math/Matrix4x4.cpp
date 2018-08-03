#include "Matrix4x4.h"
#include "Quaternion.h"

Matrix4x4::Matrix4x4() 
{
	IdentityMatrix();
}
Matrix4x4::Matrix4x4(const Matrix4x4 &m)
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
}

Matrix4x4::Matrix4x4(float e0, float e1, float e2, float e3,
	float e4, float e5, float e6, float e7,
	float e8, float e9, float e10, float e11,
	float e12, float e13, float e14, float e15)
{
	ele[0] = e0;
	ele[1] = e1;
	ele[2] = e2;
	ele[3] = e3;
	ele[4] = e4;
	ele[5] = e5;
	ele[6] = e6;
	ele[7] = e7;
	ele[8] = e8;
	ele[9] = e9;
	ele[10] = e10;
	ele[11] = e11;
	ele[12] = e12;
	ele[13] = e13;
	ele[14] = e14;
	ele[15] = e15;
}	

void Matrix4x4::Transpose()
{
	float tmp;

	tmp = ele[1];
	ele[1] = ele[4];
	ele[4] = tmp;

	tmp = ele[2];
	ele[2] = ele[8];
	ele[8] = tmp;

	tmp = ele[3];
	ele[3] = ele[12];
	ele[12] = tmp;

	tmp = ele[6];
	ele[6] = ele[9];
	ele[9] = tmp;

	tmp = ele[7];
	ele[7] = ele[13];
	ele[13] = tmp;

	tmp = ele[11];
	ele[11] = ele[14];
	ele[14] = tmp;	
}
#define Det3x3(a, b, c, d, e, f, g, h, i) (a*e*i + b*f*g + c*d*h - c*e*g - a*f*h - b*d*i)
float Matrix4x4::Inverse()
{
	float a1, a2, a3, a4, b1, b2, b3, b4;
	float c1, c2, c3, c4, d1, d2, d3, d4;		

	a1 = ele[0];
	b1 = ele[1];
	c1 = ele[2];
	d1 = ele[3];

	a2 = ele[4];
	b2 = ele[5];
	c2 = ele[6];
	d2 = ele[7];

	a3 = ele[8];
	b3 = ele[9];
	c3 = ele[10];
	d3 = ele[11];

	a4 = ele[12];
	b4 = ele[13];
	c4 = ele[14];
	d4 = ele[15];

	float det1 = Det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
	float det2 = Det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
	float det3 = Det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
	float det4 = Det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);
	float det = a1*det1 - b1*det2 + c1*det3 - d1*det4;	
	if(det == 0.0f)    return 0.0f;	

	float invdet = 1.0f / det;

	ele[0] =  det1*invdet;
	ele[4] = -det2*invdet;
	ele[8] =  det3*invdet;
	ele[12] = -det4*invdet;

	ele[1] = -Det3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4)*invdet;
	ele[5] =  Det3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4)*invdet;
	ele[9] = -Det3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4)*invdet;
	ele[13] =  Det3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4)*invdet;

	ele[2] =  Det3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4)*invdet;
	ele[6] = -Det3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4)*invdet;
	ele[10]=  Det3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4)*invdet;
	ele[14]= -Det3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4)*invdet;

	ele[3]= -Det3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3)*invdet;
	ele[7]=  Det3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3)*invdet;
	ele[11]= -Det3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3)*invdet;
	ele[15]=  Det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3)*invdet;	

	return det;
}

Matrix4x4 Matrix4x4::TranslateMatrix(float tx,float ty,float tz)
{
	Matrix4x4 tMat;
	tMat[12]=tx; tMat[13]=ty; tMat[14]=tz;
	return tMat;
}
Matrix4x4 Matrix4x4::ScaleMatrix(float sx,float sy,float sz)
{
	Matrix4x4 sMat;
	sMat[0]=sx; sMat[5]=sy; sMat[10]=sz;
	return sMat;
}
Matrix4x4 Matrix4x4::RotateMatrix(float angle,float x,float y,float z)
{
	Quaternion q(x,y,z,angle);
	return q.ToMatrix();
}
