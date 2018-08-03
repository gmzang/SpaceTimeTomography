// VolumeData.h

#ifndef _VolumeData_h
#define _VolumeData_h

#include <stdio.h>


#include "CImg.h"
////#include <map>
////#include"tv.h"
////
using namespace cimg_library;
////
////namespace cl = cimg_library;
////#define cimg_use_tiff
////typedef float				real;
////typedef geom::vec3<real>	real3;
////typedef vec2<real>	real2;
////typedef cl::CImg<float>		image;
////typedef cl::CImg<int>		table;
////typedef CImgList<float>  imagelist;

class VolumeData
{
public:
	//struct Voxel
	//{

	//};
	
//	static const DataType dtFloat32,dtInt8,dtUint8,dtInt16,dtUint16,dtInt32,dtUint32;
	//GetSliceIndex
	//GetVoxelIndex
	VolumeData();
	VolumeData(int d1,int d2, int d3, float sp1,float sp2,float sp3,float *data);
	virtual ~VolumeData();

	void Clean();

	//DataType GetDataType() const;
	//void SetDataType(DataType type);

	void GetDimensions(int *dims) const;
	void SetDimensions(const int *dims);

	void GetSpacings(float *spacings) const;
	void SetSpacings(const float *spacings) ;
	
	void SetProjectionMode(bool mode)
	{
		m_IsBackProjectMode=mode;

	}
	bool GetProjectionMode(){return m_IsBackProjectMode;}
	//void* GetData() const;
	//float * GetData(int iframe) const;

	//float * GetData(int iframe) const;

	CImg<float>  GetData(int iframe) ;
	//CImg<float> GetData();
	//float * GetData()const;
	//void SetData(float *_data) ;
	void SetData(CImgList<float> _vol_list);
	void setParas(int *vols,double vs, int nframes);
	//float *GetTotalWeight() const;
	//float *GetContributions() const;
	//float * GetCorrection() const;
	//bool SetCorrection(float * corr); 
	
	//void GetMinMaxValue(double &minVal, double &maxVal);

	//void WriteRawFile(FILE *fp) const;
	//void ReadRawFile(FILE *fp);
	//void ReadRawFile( );

	
	//void WriteStructuredFile(const char* filename) const;
	void ReadStructuredFile();
	//void ReadStructuredFile(const char* filename);
	void _allocate();
	bool m_IsBackProjectMode;
	//float *m_EachVoxel_TotalWeight;
	//float *m_Voxel_Contributions;
	//float *m_correction;

	int m_Dims[3];
	float m_Spacings[3];

	//void *m_Data;
	//float *m_Data;
	int m_volX;
	int m_volY;
	int m_volZ;
	double Vol_Spacing;
	CImgList<float>  m_Volumelist;
	CImgList<float>  m_Flowlist;
	int m_nframes;
};

#endif
