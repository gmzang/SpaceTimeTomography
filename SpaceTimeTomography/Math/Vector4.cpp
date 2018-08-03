#include "Vector4.h"
Vector4::Vector4()
{
	ele[0] = 0.0f;
	ele[1] = 0.0f;
	ele[2] = 0.0f;
	ele[3] = 1.0f;
}

Vector4::Vector4(const Vector4 &v)
{
	ele[0] = v.ele[0];
	ele[1] = v.ele[1];
	ele[2] = v.ele[2];
	ele[3] = v.ele[3];
}

Vector4::Vector4(float x, float y, float z, float w)
{
	ele[0] = x;
	ele[1] = y;
	ele[2] = z;
	ele[3] = w;
}

Vector4::Vector4(float srcVector[4])
{
	ele[0] = srcVector[0];
	ele[1] = srcVector[1];
	ele[2] = srcVector[2];
	ele[3] = srcVector[3];
}	
