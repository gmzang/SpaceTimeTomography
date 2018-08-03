#include "RayCaster.h"
#include "Vector4.h"
#include "Vector3.h"
#include<iostream>
using namespace std;
// for data access module.
static inline float GetScalarValue(Vector4 pos);
static inline Vector4 GetNormal(Vector4 pos);

// glboal params, to store parames temporarily
static float sampleDistance;
static float isovalue;
static Vector4 baseColor;
static Vector3 ambient;
static Vector3 diffuse;
static Vector3 specular;
static float specular_power;
static Vector4 light_direction;
static float spec_threshold;
static int bsize[3];

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

static inline void GetMinMax(int x,int y,int z,float &minV,float &maxV);
//
static inline int ConvertToVoxelIndex(int blockIndex, int dimId)
{
	return blockIndex*bsize[dimId];
}

static inline int ConvertToBlockIndex(int voxelIndex, int dimId)
{
	return voxelIndex / bsize[dimId];
}
//core functions, for generating  each ray  (from X-Ray source to measurement image) 
static inline double CastingRay(const Vector4& voxel_zero,const Vector4& voxel_unit,const Vector4& CameraDirection,float tnear,float tfar,float &SumWeight_Ray, int indexOfImage)
{

	
	
	int index=0;
	float density=0.0f;
	float t1=tnear;
	Vector4 pos1=voxel_zero+voxel_unit*t1;

	float tMaxX,tMaxY,tMaxZ;
	float tDeltaX,tDeltaY,tDeltaZ;

	//Vector4 acc(0.0f,0.0f,0.0f,0.0f);
	Vector4 voxel_step=voxel_unit*sampleDistance;



	float scalar1=GetScalarValue(pos1,SumWeight_Ray,indexOfImage);
	density+=scalar1*sampleDistance;
	int x,y,z;
	int stepX,stepY,stepZ;
	if (voxel_unit[0]==0)
	{
		x=ConvertToBlockIndex(((int)ceilf(pos1[0]-0.5f)-1),0);
		stepX=0;
		tMaxX=tfar;
	}
	else if (voxel_unit[0]<0)
	{
		x=ConvertToBlockIndex(((int)ceilf(pos1[0]-0.5f)-1),0);
		stepX=-1;
		tMaxX=(ConvertToVoxelIndex(x,0)+0.5f-voxel_zero[0])/voxel_unit[0];
		tDeltaX=-ConvertToVoxelIndex(1,0)/voxel_unit[0];
	}
	else
	{
		x=ConvertToBlockIndex((int)floorf(pos1[0]-0.5f),0);
		stepX=1;
		tMaxX=(ConvertToVoxelIndex(x+1,0)+0.5f-voxel_zero[0])/voxel_unit[0];
		tDeltaX=ConvertToVoxelIndex(1,0)/voxel_unit[0];
	}

	if (voxel_unit[1]==0)
	{
		y=ConvertToBlockIndex(((int)ceilf(pos1[1]-0.5f)-1),1);
		stepY=0;
		tMaxY=tfar;
	}
	else if (voxel_unit[1]<0)
	{
		y=ConvertToBlockIndex(((int)ceilf(pos1[1]-0.5f)-1),1);
		stepY=-1;
		tMaxY=(ConvertToVoxelIndex(y,1)+0.5f-voxel_zero[1])/voxel_unit[1];
		tDeltaY=-ConvertToVoxelIndex(1,1)/voxel_unit[1];
	}
	else
	{
		y=ConvertToBlockIndex((int)floorf(pos1[1]-0.5f),1);
		stepY=1;
		tMaxY=(ConvertToVoxelIndex(y+1,1)+0.5f-voxel_zero[1])/voxel_unit[1];
		tDeltaY=ConvertToVoxelIndex(1,1)/voxel_unit[1];
	}

	if (voxel_unit[2]==0)
	{
		z=ConvertToBlockIndex(((int)ceilf(pos1[2]-0.5f)-1),2);
		stepZ=0;
		tMaxZ=tfar;
	}
	else if (voxel_unit[2]<0)
	{
		z=ConvertToBlockIndex(((int)ceilf(pos1[2]-0.5f)-1),2);
		stepZ=-1;
		tMaxZ=(ConvertToVoxelIndex(z,2)+0.5f-voxel_zero[2])/voxel_unit[2];
		tDeltaZ=-ConvertToVoxelIndex(1,2)/voxel_unit[2];
	}
	else
	{
		z=ConvertToBlockIndex((int)floorf(pos1[2]-0.5f),2);
		stepZ=1;
		tMaxZ=(ConvertToVoxelIndex(z+1,2)+0.5f-voxel_zero[2])/voxel_unit[2];
		tDeltaZ=ConvertToVoxelIndex(1,2)/voxel_unit[2];
	}

	float t_b;

	
	while (t1<tfar)
	{
		t_b=min(min(tMaxX,tMaxY),tMaxZ);

	
		if (t_b>tfar) t_b=tfar;

		float t=ceilf((t1-tnear)/sampleDistance)*sampleDistance+tnear;
		Vector4 pos=voxel_zero + voxel_unit * t;		

		
		while(t<t_b)
		{
			//Vector4 color=GetShadedColor(pos,CameraDirection);
			density+=GetScalarValue(pos,SumWeight_Ray,indexOfImage)*sampleDistance;
			//float alpha=color[3];
			//color[3]=1.0f;
			//color*=alpha;		
			//acc+=(1.0f-acc[3])*color;				
			t += sampleDistance;
			pos += voxel_step;
			index++;
		}
		

		if (tMaxX<tMaxY)
		{
			if (tMaxX<tMaxZ) 
			{
				tMaxX+=tDeltaX;
				x+=stepX;
			}
			else 
			{
				tMaxZ+=tDeltaZ;
				z+=stepZ;
			}
		}
		else
		{
			if (tMaxY<tMaxZ) 
			{
				tMaxY+=tDeltaY;
				y+=stepY;
			}
			else 
			{
				tMaxZ+=tDeltaZ;
				z+=stepZ;
			}
		}
		t1=t_b;
	}


	
	return density;

}

// prepare for rendering, called by render_core
static bool RayCasterPrepare(RayCaster *RayCaster, int size[])
{
	if (!RayCaster) return false;

	
	sampleDistance=RayCaster->GetSampleDistance();
	
	baseColor[3]=1.0f;
	bsize[0] = size[0];
	bsize[1] = size[1];
	bsize[2] = size[2];
	float lightColor[3];

	light_direction[3]=0.0f;

	return true;
}



