#include "Vector3.h"
Vector3::Vector3()
{
	ele[0] = 0.0f;
	ele[1] = 0.0f;
	ele[2] = 0.0f;
}

Vector3::Vector3(const Vector3 &v)
{
	ele[0] = v.ele[0];
	ele[1] = v.ele[1];
	ele[2] = v.ele[2];
}

Vector3::Vector3(float x, float y, float z)
{
	ele[0] = x;
	ele[1] = y;
	ele[2] = z;
}

Vector3::Vector3(float srcVector[3])
{
	ele[0] = srcVector[0];
	ele[1] = srcVector[1];
	ele[2] = srcVector[2];
}	
