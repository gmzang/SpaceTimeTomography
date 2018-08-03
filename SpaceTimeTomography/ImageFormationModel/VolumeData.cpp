// VolumeData.cpp

#include "VolumeData.h"
#include <malloc.h>
#include <fstream>

#include <iostream>

using namespace  std;



VolumeData::VolumeData()
{
	m_Dims[0]=128;
	m_Dims[1]=128;
	m_Dims[2]=128;
	m_Spacings[0]=m_Spacings[1]=m_Spacings[2]=1.0f;

	int _index=0 ;
	
	m_IsBackProjectMode=false;

}


 VolumeData::VolumeData(int d1,int d2, int d3, float sp1,float sp2,float sp3,float *data)
{
	int d[]={d1,d2,d3};
	this->SetDimensions(d);
	float sp[]={sp1,sp2,sp3};
	this->SetSpacings(sp);
	this->_allocate();
}



 void VolumeData::setParas(int *vols, double vs,int nframes)
 {
	 m_volX= vols[0];
	 m_volY= vols[1];
	 m_volZ= vols[2];
	 Vol_Spacing=vs;
	 m_nframes = nframes;

 }
 
VolumeData::~VolumeData()
{
	Clean();
}

void VolumeData::Clean()
{

	
}


void VolumeData::_allocate()
{
	
}



void VolumeData::GetDimensions(int *dims) const
{
	dims[0]=m_Dims[0]; dims[1]=m_Dims[1]; dims[2]=m_Dims[2];
}

void VolumeData::SetDimensions(const int *dims)
{
	m_Dims[0]=dims[0]; m_Dims[1]=dims[1]; m_Dims[2]=dims[2];
}

void VolumeData::GetSpacings(float *spacings) const
{
	spacings[0]=m_Spacings[0]; spacings[1]=m_Spacings[1]; spacings[2]=m_Spacings[2];
}

void VolumeData::SetSpacings(const float *spacings)
{
	m_Spacings[0]=spacings[0]; m_Spacings[1]=spacings[1]; m_Spacings[2]=spacings[2];
}


CImg<float> VolumeData::GetData(int iframe)
{
	return m_Volumelist[iframe];
	//return 0;
}

void  VolumeData::SetData(CImgList<float> _vol_list)
{
	//m_Data = _data;
	//return true;
	m_Volumelist= _vol_list;
}



void VolumeData::ReadStructuredFile()
{
	m_Dims[0]=m_volX;
	m_Dims[1]=m_volY;
	m_Dims[2]=m_volZ;
	m_Spacings[0]=m_Spacings[1]=m_Spacings[2]=Vol_Spacing;  

	m_Flowlist = CImgList<float>(m_nframes, m_Dims[0], m_Dims[1], m_Dims[2], 3, 0);

	char basename[1024];
	sprintf(basename, "Input");
	//rose.-000-frame-9iter
	for (int i = 0; i < m_nframes; i++)
	{
		char temp[200];
		
		sprintf(temp, "%s/rose.-%03d-frame-9iter.tiff", basename, i);

		m_Volumelist.insert(CImg<float>(temp));

	}
	// init the time steps as n+1 
	// in case potential boundary issue.
	m_Volumelist.insert(m_Volumelist[m_nframes - 1]);

	
}

