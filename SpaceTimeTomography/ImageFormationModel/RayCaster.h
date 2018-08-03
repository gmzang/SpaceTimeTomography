// raycaster.h
#ifndef _RayCaster_h
#define _RayCaster_h

#include "Matrix4x4.h"



class RayCaster
{
public:
	RayCaster();
	virtual ~RayCaster();

	// interpolating distance
	void  SetSampleDistance(float sd) {m_SampleDistance = sd;}
	float GetSampleDistance() { return m_SampleDistance;}

	void MakeReference();
	void RotateLightDirection(const Matrix4x4& matrix);

private:

	// ray direction
	Vector4 m_LightDirection;
	Vector4 m_OldLightDirection;
	
	// Sample Distance
	float m_SampleDistance;
};

#endif


