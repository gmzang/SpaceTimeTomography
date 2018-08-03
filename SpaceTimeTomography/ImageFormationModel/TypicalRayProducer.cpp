// TypicalRayProducer.cpp

#include "TypicalRayProducer.h"
#include "GL/glew.h"
//#include "ImageIO.h"
#include <iostream>
using namespace std;
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif



static int IMAGE_WITH;
static int IMAGE_HEIGHT;
static int Vol_X;
static int Vol_Y;
static int Vol_Z;
static double m_ds;
bool isRaybased =false;

TypicalRayProducer::TypicalRayProducer(int vX, int vY, int vZ, int imgW, int imgH, double ds,  const char * _bp)
{
	
	m_ds=ds;

	m_ImageInUseSize[0]      = 0;
	m_ImageInUseSize[1]      = 0;
	Vol_X=vX;
	Vol_Y=vY;
	Vol_Z=vZ;
	IMAGE_WITH=imgW;
	IMAGE_HEIGHT=imgH;
	m_VolSize[0]=m_VolSize[1]=m_VolSize[2]=0;
	cout<<"BB: image: "<<imgW<<" "<<imgH<<endl;
	if(strcmp(_bp, "Raybased") == 0)
	{
		isRaybased=true;
	}
	m_ImageData=0;
	m_WeightOfImage=0;

	m_BackProjection=false;
	
	m_WeightOfImage=new float[IMAGE_WITH*IMAGE_HEIGHT];
	
	m_ImageData=new float[IMAGE_WITH*IMAGE_HEIGHT];
	m_imageraylength=new float[IMAGE_WITH*IMAGE_HEIGHT];

	m_Correction=new float[IMAGE_WITH*IMAGE_HEIGHT];

	m_EachVoxel_TotalWeight=new float[Vol_X*Vol_Y*Vol_Z];
	m_Voxel_Contributions=new float[Vol_X*Vol_Y*Vol_Z];

	memset(m_imageraylength, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);

	memset(m_ImageData, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);
	memset(m_WeightOfImage, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);

	memset(m_Correction, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);
	memset(m_EachVoxel_TotalWeight, 0, sizeof(float)*Vol_X*Vol_Y*Vol_Z);
	memset(m_Voxel_Contributions, 0, sizeof(float)*Vol_X*Vol_Y*Vol_Z);
}

TypicalRayProducer::~TypicalRayProducer()
{
	delete[] m_ImageData;
	delete[] m_WeightOfImage;
}




void TypicalRayProducer::Initial()
{

	m_BackProjection=false;


	memset(m_imageraylength, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);
	memset(m_ImageData, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);
	memset(m_WeightOfImage, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);

	memset(m_Correction, 0, sizeof(float)*IMAGE_WITH*IMAGE_HEIGHT);
	/////////
	////if voxelbased, we do not need this initialization.
	if(isRaybased)
	{
	memset(m_EachVoxel_TotalWeight, 0, sizeof(float)*Vol_X*Vol_Y*Vol_Z);
	memset(m_Voxel_Contributions, 0, sizeof(float)*Vol_X*Vol_Y*Vol_Z);
	}
}



void TypicalRayProducer::SetViewPort(int width,int height)
{
	m_ViewPort[0]=width;
	m_ViewPort[1]=height;
	cout<<"TypicalRayProducer::SetViewPort "<<width<<" "<<height<<endl;
}
void TypicalRayProducer::SetVolSize(const int volDims[3],const float volSpacings[3])
{
	m_VolDims[0]=volDims[0];
	m_VolDims[1]=volDims[1];
	m_VolDims[2]=volDims[2];

	m_VolSpacings[0]=volSpacings[0];
	m_VolSpacings[1]=volSpacings[1];
	m_VolSpacings[2]=volSpacings[2];

	m_VolSize[0]=(volDims[0]-1)*volSpacings[0];
	m_VolSize[1]=(volDims[1]-1)*volSpacings[1];
	m_VolSize[2]=(volDims[2]-1)*volSpacings[2];

}

void TypicalRayProducer::GetVolSize(int volDims[3],float volSpacings[3])
{
	volDims[0]=m_VolDims[0];
	volDims[1]=m_VolDims[1];
	volDims[2]=m_VolDims[2];

	volSpacings[0]=m_VolSpacings[0];
	volSpacings[1]=m_VolSpacings[1];
	volSpacings[2]=m_VolSpacings[2];
}

void TypicalRayProducer::SetViewMatrix(const Matrix4x4& ViewMatrix)
{
	m_ViewMatrix=ViewMatrix;
}
void TypicalRayProducer::SetProjectionMatrix(const Matrix4x4& ProjectionMatrix)
{
	m_ProjectionMatrix=ProjectionMatrix;
}

void TypicalRayProducer::_calculateMatrix()
{


	Matrix4x4 modelViewMatrix
		=m_ViewMatrix*Matrix4x4::TranslateMatrix(-m_VolSize[0]*0.5f,-m_VolSize[1]*0.5f,-m_VolSize[2]*0.5f); // 
	

	m_model_view=modelViewMatrix;

	Matrix4x4 invertModelViewMatrix=modelViewMatrix;
	invertModelViewMatrix.Inverse();

	m_CameraPos[0]=invertModelViewMatrix[12]/invertModelViewMatrix[15];
	m_CameraPos[1]=invertModelViewMatrix[13]/invertModelViewMatrix[15];
	m_CameraPos[2]=invertModelViewMatrix[14]/invertModelViewMatrix[15];

	// model to world
	m_ModelToVolumeMatrix.ele[0]  = 1.0f/m_VolSpacings[0];
	m_ModelToVolumeMatrix.ele[1]  = m_ModelToVolumeMatrix.ele[2] = m_ModelToVolumeMatrix.ele[3] = m_ModelToVolumeMatrix.ele[4] = 0.0f;
	m_ModelToVolumeMatrix.ele[5]  = 1.0f/m_VolSpacings[1];
	m_ModelToVolumeMatrix.ele[6]  = m_ModelToVolumeMatrix.ele[7] = m_ModelToVolumeMatrix.ele[8] = m_ModelToVolumeMatrix.ele[9] = 0.0f;
	m_ModelToVolumeMatrix.ele[10] = 1.0f/m_VolSpacings[2];
	m_ModelToVolumeMatrix.ele[11] = 0.0f;
	m_ModelToVolumeMatrix.ele[12] = 0.0f;
	m_ModelToVolumeMatrix.ele[13] = 0.0f;
	m_ModelToVolumeMatrix.ele[14] = 0.0f;
	m_ModelToVolumeMatrix.ele[15] = 1.0f;

	m_VolumeToModelMatrix=m_ModelToVolumeMatrix;
	m_VolumeToModelMatrix.Inverse();
	// view to camera
	Matrix4x4 ViewToPixelsMatrix;
	ViewToPixelsMatrix.ele[0]  =0.5f * (m_ViewPort[0]-1);
	ViewToPixelsMatrix.ele[1]  = ViewToPixelsMatrix.ele[2] = ViewToPixelsMatrix.ele[3] = ViewToPixelsMatrix.ele[4] = 0.0f;

	

	ViewToPixelsMatrix.ele[5]  =0.5f * (m_ViewPort[1]-1);
	ViewToPixelsMatrix.ele[6]  = ViewToPixelsMatrix.ele[7] = ViewToPixelsMatrix.ele[8] = ViewToPixelsMatrix.ele[9] = 0.0f;

	//ViewToPixelsMatrix.ele[3]=m_offx;
	//ViewToPixelsMatrix.ele[7]=m_offy;

	ViewToPixelsMatrix.ele[10] = 0.5f;
	//ViewToPixelsMatrix.ele[10] = 0.0f;
	ViewToPixelsMatrix.ele[11] = 0.0f;
	

	ViewToPixelsMatrix.ele[12] = ViewToPixelsMatrix.ele[0]; 
	ViewToPixelsMatrix.ele[13] = ViewToPixelsMatrix.ele[5];
	ViewToPixelsMatrix.ele[14] = ViewToPixelsMatrix.ele[10];
	ViewToPixelsMatrix.ele[15] = 1.0;
	

	m_PixelsToViewMatrix=ViewToPixelsMatrix;
	m_PixelsToViewMatrix.Inverse();

	// pixel to model
	m_ModelToPixelsMatrix = ViewToPixelsMatrix*m_ProjectionMatrix*modelViewMatrix;
	//m_ModelToPixelsMatrix = m_ProjectionMatrix*modelViewMatrix;

	m_PixelsToModelMatrix=m_ModelToPixelsMatrix;
	m_PixelsToModelMatrix.Inverse();


}




bool TypicalRayProducer::_calculateBounds()
{
	float minX, minY, maxX, maxY;
	Vector4 ModelPoint;   
	Vector4 PixelPoint;
	
	
	
	m_ImageOrigin[0]=0;//IMAGE_WITH*(m_ds-1.0f)/2.0f;
	m_ImageOrigin[1]=0;//IMAGE_HEIGHT*(m_ds-1.0f)/2.0f;
	m_ImageInUseSize[0]=IMAGE_WITH;
	m_ImageInUseSize[1]=IMAGE_HEIGHT;

	

	return true;
}

void TypicalRayProducer::UpdateMatrixAndBounds()
{
	_calculateMatrix();
	if (!_calculateBounds()) return;

	

}

void TypicalRayProducer::GetImageInUseSize(int imageSize[2])
{
	imageSize[0]=m_ImageInUseSize[0];
	imageSize[1]=m_ImageInUseSize[1];
}



void TypicalRayProducer::GetImageOrigin(int imageOrigin[2])
{
	imageOrigin[0]=m_ImageOrigin[0];
	imageOrigin[1]=m_ImageOrigin[1];
}
void TypicalRayProducer::RenderImage()  
{
	
	if (!m_ImageData) return;
	//if (!m_UC_data) return;
	
	glEnable(GL_BLEND);
	
	glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA);
	
	glWindowPos2i(m_ImageOrigin[0], m_ImageOrigin[1]);


	

}
