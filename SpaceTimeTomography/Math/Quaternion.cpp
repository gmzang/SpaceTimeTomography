#include "Quaternion.h"

Quaternion::Quaternion()
{
	Identity();
}
Quaternion::Quaternion(float ax, float ay, float az, float angle)
{
	float xx,yy,zz;			
	float rad, scale;		
	float length;

	length = sqrtf(ax*ax+ay*ay+az*az);

	if (0.0f==length)			// if axis is zero, then return Quaternion (1,0,0,0)
	{
		w	= 1.0f;
		x	= 0.0f;
		y	= 0.0f;
		z	= 0.0f;		
		return;
	}		

	float lenInv=1.0f/length;
	xx = ax *lenInv;
	yy = ay *lenInv;
	zz = az *lenInv;	    

	rad		= angle / 2;		
	w		= cosf(rad);		
	scale	= sinf(rad);

	x = xx * scale;
	y = yy * scale;
	z = zz * scale;
	//Change by GM, test rotation.
	//x=0;
	//y=0;
	//z=0;
}
