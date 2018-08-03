// RayCaster.cpp

#include "RayCaster.h"
#include <iostream>
using namespace std;
RayCaster::RayCaster()
{
	m_LightDirection[0] = 1.0f;
	m_LightDirection[1] = 0.0f;
	m_LightDirection[2] = 0.0f;
	m_LightDirection[3] = 0.0f;


//	m_LightIntensity = 1.0f;


//	m_SampleDistance = 1.0f;   
	m_SampleDistance = 1.0f;   

	

}
RayCaster::~RayCaster()
{

}



void RayCaster::MakeReference()
{
	m_OldLightDirection=m_LightDirection;
}

void RayCaster::RotateLightDirection(const Matrix4x4& matrix)
{
	m_LightDirection=matrix*m_OldLightDirection;
}


