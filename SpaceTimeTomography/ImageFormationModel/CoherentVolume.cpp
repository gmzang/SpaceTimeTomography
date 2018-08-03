#include "CoherentVolume.h"
#include "VolumeData.h"
#include <memory.h>
#include <iostream>
using namespace std;

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif


class Block
{
public:
	Block()	{	}
	virtual ~Block() {	}
	virtual void GetMinMax(float& minV,float& maxV)=0;
	virtual void SetData(float *data,int xSize,int ySize,int zSize,int yInc,int zInc)=0;
	virtual float GetValue(int x,int y,int z)=0;

protected:
	int m_yInc,m_zInc;
};

//template<typename T float>
class TBlock : public Block
{
public:
	TBlock()
	{
		m_minV=0;
		m_maxV=0;
		m_data=NULL;
	}
	virtual ~TBlock()
	{
		delete[] m_data;
	}
	virtual void GetMinMax(float& minV,float& maxV)
	{
		minV=(float)m_minV;
		maxV=(float)m_maxV;
	}
	virtual void SetData(float *data,int xSize,int ySize,int zSize,int yInc,int zInc)
	{
		delete[] m_data;
		m_data=new float[xSize*ySize*zSize];
		m_yInc=xSize;
		m_zInc=xSize*ySize;
		m_maxV=m_minV=*(float*)(data);
		float* pData10=m_data;
		float* pData00=(float*)(data);
		int i,j,k;
	//	cout<<"m_yInc,m_zInc"<<m_yInc<<" "<<m_zInc <<endl;
		for (i=0;i<zSize;i++,pData00+=zInc,pData10+=m_zInc)
		{
			float* pData01=pData00;
			float* pData11=pData10;
			for (j=0;j<ySize;j++,pData01+=yInc,pData11+=m_yInc)
			{
				memcpy(pData11,pData01,sizeof(float)*xSize);
				for (k=0;k<xSize;k++) 
				{
					m_minV=min(m_minV,pData11[k]);
					m_maxV=max(m_maxV,pData11[k]);
				}
			}
		}		
		if (m_minV==m_maxV) 
		{
			delete[] m_data;
			m_data=NULL;
		}
	}
	virtual float GetValue(int x,int y,int z)
	{
		if (!m_data) return (float)m_minV;
		else return (float)m_data[z*m_zInc+y*m_yInc+x];
	}

private:
	float* m_data;
	float m_minV,m_maxV;
};

class CoherentVolume::BlockArray 
{
public:
	BlockArray(){}
	virtual ~BlockArray(){}
	virtual Block& operator[](int index)=0;
};

//template<typename T>
class TBlockArray : public CoherentVolume::BlockArray 
{
public:
	TBlockArray(int size)
	{
		m_Data=new TBlock[size];
	}
	virtual ~TBlockArray()
	{
		delete[] m_Data;
	}
	virtual TBlock& operator[](int index)
	{
		return m_Data[index];
	}
private:
	TBlock* m_Data;
};


CoherentVolume::CoherentVolume(int x, int y, int z)
{
	//m_BlockSize[0]=m_BlockSize[1]=m_BlockSize[2]=64;

	//m_BlockSize[0=this->GetDimensions()
	
	//m_BlockSize[0]=208;
	//m_BlockSize[1]=256;
	//m_BlockSize[2]=225;

	//m_BlockSize[0]=256;
	//m_BlockSize[1]=256;
	//m_BlockSize[2]=256;
	m_BlockSize[0]=x;
	m_BlockSize[1]=y;
	m_BlockSize[2]=z;
	m_BlockArray=NULL;
}

CoherentVolume::~CoherentVolume()
{
	delete m_BlockArray;
}

void CoherentVolume::Initialize()
{
	delete m_BlockArray;
	m_BlockArray=NULL;
}


//DataType CoherentVolume::GetDataType() const
//{
//	return m_DataType;
//}
//
//void CoherentVolume::SetDataType(DataType type)
//{
//	m_DataType=type;
//}

void CoherentVolume::AllocateBlocks()
{
	Initialize();
	 m_BlockArray=new TBlockArray(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
	/*if (m_DataType.isFloat) m_BlockArray=new TBlockArray<float>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
	else if (m_DataType.isSigned)
	{
		if (m_DataType.bitsPerSample==8) m_BlockArray=new TBlockArray<char>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
		else if (m_DataType.bitsPerSample==16) m_BlockArray=new TBlockArray<short>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
		else if (m_DataType.bitsPerSample==32) m_BlockArray=new TBlockArray<int>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
	}
	else
	{
		if (m_DataType.bitsPerSample==8) m_BlockArray=new TBlockArray<unsigned char>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
		else if (m_DataType.bitsPerSample==16) m_BlockArray=new TBlockArray<unsigned short>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
		else if (m_DataType.bitsPerSample==32) m_BlockArray=new TBlockArray<unsigned int>(m_BlockNum[0]*m_BlockNum[1]*m_BlockNum[2]);
	}*/
}

void CoherentVolume::SetBlockSlice(void *sliceData,int sliceID,int dimz)
{
	if (!m_BlockArray) return;
	if (dimz<0) dimz=m_BlockSize[2];

	float* pData=(float*)sliceData;

	int i,j;
	for (i=0;i<m_BlockNum[1];i++,pData+=m_BlockSize[1]*m_Dimensions[0]*(sizeof(float)))
	{
		float* pData1=pData;
		for (j=0;j<m_BlockNum[0];j++,pData1+=(sizeof(float))*m_BlockSize[0])
		{
			(*m_BlockArray)[(sliceID*m_BlockNum[1]+i)*m_BlockNum[0]+j].SetData(
				pData1,min(m_BlockSize[0],m_Dimensions[0]-m_BlockSize[0]*j),
				min(m_BlockSize[1],m_Dimensions[1]-m_BlockSize[1]*i),min(m_BlockSize[2],dimz),
				m_Dimensions[0],m_Dimensions[0]*m_Dimensions[1]);
		}

	} 

	/*for (i=0;i<m_BlockNum[1];i++,pData+=m_BlockSize[1]*m_Dimensions[0]*(m_DataType.bitsPerSample>>3))
	{
		unsigned char* pData1=pData;
		for (j=0;j<m_BlockNum[0];j++,pData1+=(m_DataType.bitsPerSample>>3)*m_BlockSize[0])
		{
			(*m_BlockArray)[(sliceID*m_BlockNum[1]+i)*m_BlockNum[0]+j].SetData(
				pData1,min(m_BlockSize[0],m_Dimensions[0]-m_BlockSize[0]*j),
				min(m_BlockSize[1],m_Dimensions[1]-m_BlockSize[1]*i),min(m_BlockSize[2],dimz),
				m_Dimensions[0],m_Dimensions[0]*m_Dimensions[1]);
		}

	} */
}

void CoherentVolume::SetData(float *data)
{
	if (!m_BlockArray) return;
	
	float* pData=(float*)data;

	int i;


	for (i=0;i<m_BlockNum[2];i++,pData+=m_BlockSize[2]*m_Dimensions[0]*m_Dimensions[1]*(sizeof(float)))
		SetBlockSlice(pData,i,min(m_BlockSize[2],m_Dimensions[2]-m_BlockSize[2]*i));

	/*for (i=0;i<m_BlockNum[2];i++,pData+=m_BlockSize[2]*m_Dimensions[0]*m_Dimensions[1]*(m_DataType.bitsPerSample>>3))
		SetBlockSlice(pData,i,min(m_BlockSize[2],m_Dimensions[2]-m_BlockSize[2]*i));
	*/
}

//void CoherentVolume::SetData(VolumeData *vol)
//{
//	if (!vol || !vol->GetData()) return;
//	vol->GetDimensions(m_Dimensions);
//	vol->GetSpacings(m_Spacings);
//
//	m_BlockNum[0]= (m_Dimensions[0]-1)/m_BlockSize[0]+1;
//	m_BlockNum[1]= (m_Dimensions[1]-1)/m_BlockSize[1]+1;
//	m_BlockNum[2]= (m_Dimensions[2]-1)/m_BlockSize[2]+1;
//
//	//this->SetDataType(vol->GetDataType());
//	AllocateBlocks();
//	SetData(vol->GetData());
//}

float CoherentVolume::GetValue(int x,int y,int z)
{
	if (!m_BlockArray) return 0.0f;

	int bx=x/m_BlockSize[0]; x=x%m_BlockSize[0];
	int by=y/m_BlockSize[1]; y=y%m_BlockSize[1];
	int bz=z/m_BlockSize[2]; z=z%m_BlockSize[2];

	return (*m_BlockArray)[bx+(by+bz*m_BlockNum[1])*m_BlockNum[0]].GetValue(x,y,z);
}

void CoherentVolume::GetMinMaxValue(double &minValue,double &maxValue)
{
	float minV,maxV;
	(*m_BlockArray)[0].GetMinMax(minV,maxV);

	int i,j,k;

	for (i=0;i<m_BlockNum[2];i++)
	{
		for (j=0;j<m_BlockNum[1];j++)
		{
			for (k=0;k<m_BlockNum[0];k++)
			{
				float bmin,bmax;
				(*m_BlockArray)[k+(j+i*m_BlockNum[1])*m_BlockNum[0]].GetMinMax(bmin,bmax);
				minV=min(minV,bmin);
				maxV=max(maxV,bmax);
				
			}
		}
	}
	minValue=minV;
	maxValue=maxV;
}
