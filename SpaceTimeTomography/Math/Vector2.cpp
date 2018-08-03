#include "Vector2.h"
Vector2::Vector2()
{
	ele[0] = 0.0f;
	ele[1] = 0.0f;
}

Vector2::Vector2(const Vector2 &v)
{
	ele[0] = v.ele[0];
	ele[1] = v.ele[1];
}

Vector2::Vector2(float x, float y)
{
	ele[0] = x;
	ele[1] = y;
}

Vector2::Vector2(float srcVector[2])
{
	ele[0] = srcVector[0];
	ele[1] = srcVector[1];
}	
