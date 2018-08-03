#ifndef __CoherentVolume_h
#define __CoherentVolume_h

#include "MyStandardIncludes.h"



class VolumeData;

class   CoherentVolume 
{
public:
	class BlockArray;

	//CoherentVolume();
	CoherentVolume(int x, int y, int z);
	virtual ~CoherentVolume();

	virtual void Initialize();

	void SetBlockSize(const int bs[3])
	{
		m_BlockSize[0]=bs[0];
		m_BlockSize[1]=bs[1];
		m_BlockSize[2]=bs[2];
	}
	void GetBlcokSize(int bs[3])
	{
		bs[0]=m_BlockSize[0];
		bs[1]=m_BlockSize[1];
		bs[2]=m_BlockSize[2];
	}

	void SetDimensions(const int dims[3])
	{
		m_Dimensions[0] = dims[0];
		m_Dimensions[1] = dims[1];
		m_Dimensions[2] = dims[2];
		m_BlockNum[0]= (m_Dimensions[0]-1)/m_BlockSize[0]+1;
		m_BlockNum[1]= (m_Dimensions[1]-1)/m_BlockSize[1]+1;
		m_BlockNum[2]= (m_Dimensions[2]-1)/m_BlockSize[2]+1;
	}

	void GetDimensions(int dims[3]) const 
	{
		dims[0] = m_Dimensions[0];
		dims[1] = m_Dimensions[1];
		dims[2] = m_Dimensions[2];
	}

	void SetSpacings(const float s[3])
	{
		m_Spacings[0] = s[0];
		m_Spacings[1] = s[1];
		m_Spacings[2] = s[2];
	}

	void GetSpacings(float s[3]) const 
	{
		s[0] = m_Spacings[0];
		s[1] = m_Spacings[1];
		s[2] = m_Spacings[2];
	}

	void GetBlockNum(int bn[3]) 
	{
		bn[0]=m_BlockNum[0];
		bn[1]=m_BlockNum[1];
		bn[2]=m_BlockNum[2];
	}

	//DataType GetDataType() const;
	//void SetDataType(DataType type); 

	void AllocateBlocks();

	void SetBlockSlice(void *sliceData,int sliceID,int dimz=-1);

	void SetData(float *data);

	//void SetData(VolumeData *vol);

	float GetValue(int x,int y,int z);

	void GetMinMaxValue(double &minValue,double &maxValue);

private:

	BlockArray *m_BlockArray;

	int m_BlockSize[3];
	int m_Dimensions[3];
	float m_Spacings[3];
	int m_BlockNum[3];
	//DataType m_DataType;

};

#endif
