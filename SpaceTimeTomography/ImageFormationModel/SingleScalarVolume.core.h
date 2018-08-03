// SingleScalarVolume.core.h

#include "VolumeData.h"
#include "TypicalRayProducer.h"
#include "Vector4.h"
#include<iostream>
using namespace std;
CImg<float> idata;
static float *pWeightForVoxels;
static float *pContributions;
static float *pCorrection;
static int m_dims[3];
static float m_spacings[3];
static int m_incs[3];
static bool m_IsBackProjectionMode;

class VolumeAccess
{
public:
	VolumeAccess(){}
	virtual ~VolumeAccess(){}

	// visit sample  points
	virtual float GetScalarValue(int x,int y,int z,float &SumWeight_Ray, int indexOfimage)=0;
	// get the interpolation position
	//virtual float GetScalarValue(Vector4 pos)=0;
	virtual float GetScalarValue(Vector4 pos,float &SumWeight_Ray,int indexOfimage)=0;
	// visit interpolation  kaiser bessel
	virtual float GetScalarValue_Kaiser(Vector4 pos, float theta)=0;
	// get normal
	virtual Vector4 GetNormal(Vector4 pos)=0;
};


//template <typename T>
class T_VolumeAccess : public VolumeAccess
{

public: 

	T_VolumeAccess(CImg<float> data, const int* dims, const float* spacings, bool IsBackProMode, float * top, float*down, float *correction)
	{
		//m_data=(float*)data;
		idata = data;
		pContributions=(float*)top;
		pWeightForVoxels=(float*)down;
		pCorrection=(float*)correction;
		m_IsBackProjectionMode=IsBackProMode;
		m_dims[0]=dims[0];m_dims[1]=dims[1];m_dims[2]=dims[2];
		m_spacings[0]=spacings[0];m_spacings[1]=spacings[1];m_spacings[2]=spacings[2];

		m_incs[0]=1;
		m_incs[1]=m_dims[0];
		m_incs[2]=m_dims[0]*m_dims[1];
	

	}
	virtual ~T_VolumeAccess(){}

	virtual float GetScalarValue(int x,int y,int z,float &SumWeight_Ray, int indexOfimage)
	//virtual float GetScalarValue(int x,int y,int z)
	{
		

		x=max(min(x,m_dims[0]-1),0);
		y=max(min(y,m_dims[1]-1),0);
		z=max(min(z,m_dims[2]-1),0);


		return idata(x,y,z);
	}

	virtual void UpdateContributionSum(int x,int y,int z,int iImag,float w)
	{
		x=max(min(x,m_dims[0]-1),0);
		y=max(min(y,m_dims[1]-1),0);
		z=max(min(z,m_dims[2]-1),0);

		pContributions[x+(y+z*m_dims[1])*m_dims[0]]+=pCorrection[iImag]*w;
		pWeightForVoxels[x+(y+z*m_dims[1])*m_dims[0]]+=w;



	}

	//trilinear interpolation
	virtual float GetScalarValue(Vector4 pos,float &SumWeight_Ray, int indexOfimage)
	{
		pos[0]=max(min(pos[0],m_dims[0]-1),0);
		pos[1]=max(min(pos[1],m_dims[1]-1),0);
		pos[2]=max(min(pos[2],m_dims[2]-1),0);

		float a,b,c,d;
		int ipos[3];
		float fpos1,fpos2;

		ipos[0]=(int)pos[0]; ipos[1]=(int)pos[1]; ipos[2]=(int)pos[2];

			float t0 = pos[0] - ipos[0];
			float t1 = pos[1] - ipos[1];
			float t2 = pos[2] - ipos[2];
			float M_t0 = 1.0f - t0;
			float M_t1 = 1.0f - t1;
			float M_t2 = 1.0f - t2;
			float w[8] = { M_t0*M_t1*M_t2, t0*M_t1*M_t2, M_t0*t1*M_t2, t0*t1*M_t2, M_t0*M_t1*t2, t0*M_t1*t2, M_t0*t1*t2, t0*t1*t2 };



		fpos2=pos[0]-ipos[0];
		fpos1=1.0f-fpos2;
		if (!m_IsBackProjectionMode)
			SumWeight_Ray+=1.0f;
		a= GetScalarValue(ipos[0],ipos[1],ipos[2],SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1],ipos[2],SumWeight_Ray,indexOfimage)*fpos2;
		b= GetScalarValue(ipos[0],ipos[1]+1,ipos[2],SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1]+1,ipos[2],SumWeight_Ray,indexOfimage)*fpos2;
		c= GetScalarValue(ipos[0],ipos[1],ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1],ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos2;
		d= GetScalarValue(ipos[0],ipos[1]+1,ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1]+1,ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos2;

		fpos2=pos[1]-ipos[1];
		fpos1=1.0f-fpos2;

		a=a*fpos1+b*fpos2;
		c=c*fpos1+d*fpos2;

		fpos2=pos[2]-ipos[2];
		fpos1=1.0f-fpos2;
//

			if (m_IsBackProjectionMode)
			{


				UpdateContributionSum(ipos[0],ipos[1],ipos[2],indexOfimage,w[0]);
				UpdateContributionSum(ipos[0]+1,ipos[1],ipos[2],indexOfimage,w[1]);
				UpdateContributionSum(ipos[0],ipos[1]+1,ipos[2],indexOfimage,w[2]);
				UpdateContributionSum(ipos[0]+1,ipos[1]+1,ipos[2],indexOfimage,w[3]);
				UpdateContributionSum(ipos[0],ipos[1],ipos[2]+1,indexOfimage,w[4]);
				UpdateContributionSum(ipos[0]+1,ipos[1],ipos[2]+1,indexOfimage,w[5]);
				UpdateContributionSum(ipos[0],ipos[1]+1,ipos[2]+1,indexOfimage,w[6]);
				UpdateContributionSum(ipos[0]+1,ipos[1]+1,ipos[2]+1,indexOfimage,w[7]);


				return 0;
			}
		return a*fpos1+c*fpos2;//a

	}


	virtual float GetScalarValue_Kaiser(Vector4 pos, float theta)
	{

		return 0;
	}

	virtual Vector4 GetNormal(Vector4 pos)
	{
		return Vector4(0,0,0,1);
	
	}

};

// static volume access
static VolumeAccess *va;


// core function, get the ray value from data
static inline float GetScalarValue(int x,int y,int z,float &SumWeight_Ray,int indexOfimage)
{
		if (va) return va->GetScalarValue(x,y,z,SumWeight_Ray, indexOfimage); else return 0.0f;
}


// in the case of interpolation
static inline float GetScalarValue(Vector4 pos,float &SumWeight_Ray,int indexOfimage)
{
	if (va) return va->GetScalarValue(pos,SumWeight_Ray, indexOfimage); else return 0.0f;
}

// in the case of interpolation

 static float GetScalarValue_Kaiser(Vector4 pos, float theta)
{

	if (va) return va->GetScalarValue_Kaiser(pos,theta); else return 0.0f;
}

// get normal
static inline Vector4 GetNormal(Vector4 pos)
{
	if (va) return va->GetNormal(pos); else return Vector4();
}

// Preparation
static bool SingleScalarVolumePrepare(VolumeData *volume,TypicalRayProducer *rayproducer, int frame)
{
	if (!volume) return false;
	//VolumeData::DataType type=volume->GetDataType();
//	float type=volume->GetDataType();
	va=0;
	//if(va) 

	float *data;
	
	int dims[3];
	float spacings[3];
	bool IsBackProMode;
	float *top;
	float *down;
	float *correction;

	top=rayproducer->GetContributions();
	down=rayproducer->GetTotalWeight();
	volume->GetDimensions(dims);
	volume->GetSpacings(spacings);
	IsBackProMode=volume->GetProjectionMode();
	correction=rayproducer->GetCorrection();
	 va = new T_VolumeAccess(volume->m_Volumelist[frame], dims, spacings, IsBackProMode, top, down, correction);
	


	return true;
}

// recycling
static void SingleScalarVolumeClear()
{
	if (va) delete va;
}







