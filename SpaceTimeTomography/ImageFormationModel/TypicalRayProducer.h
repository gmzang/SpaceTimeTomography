// TypicalRayProducer.h

#ifndef __TypicalRayProducer_h_h
#define __TypicalRayProducer_h_h

#include "Matrix4x4.h"

class TypicalRayProducer
{
public:
	TypicalRayProducer(int volX,int volY,int volZ, int imgW, int imgH,double ds, const char * _bp);
	void Initial();
	virtual ~TypicalRayProducer();
	void SetBackProjectionMode(){m_BackProjection=true; } 
	
	bool IsBackProjectionMode(){return m_BackProjection;}
	
	void SetViewPort(int width,int height);  
	void SetVolSize(const int volDims[3],const float volSpacings[3]);
	void GetVolSize(int volDims[3], float volSpacings[3]);

	void SetViewMatrix(const Matrix4x4& ViewMatrix);
	void SetProjectionMatrix(const Matrix4x4& ProjectionMatrix);

	void UpdateMatrixAndBounds();

	
	const Matrix4x4& GetPixelsToModelMatrix() {return m_PixelsToModelMatrix;}
	const Matrix4x4& GetModelToVolumeMatrix() {return m_ModelToVolumeMatrix;}

	
	void GetImageInUseSize(int imageSize[2]);
	void GetImageOrigin(int imageOrigin[2]);
		float *GetTotalWeight() { return m_EachVoxel_TotalWeight;}
	
	float *GetImageData() { return m_ImageData; }
	float *GetWeightOfImage(){return m_WeightOfImage;}

	float *GetImageRaylength() { return m_imageraylength; }

	float *GetContributions() {return m_Voxel_Contributions;}
	float* GetCorrection() {return m_Correction;}
		
	void RenderImage();

public:

	void _calculateMatrix();

	bool _calculateBounds();

	int m_ViewPort[2];  
	int m_VolDims[3]; 
	float m_VolSpacings[3]; 
	float m_VolSize[3]; 

	
	Matrix4x4 m_ViewMatrix;
	Matrix4x4 m_ProjectionMatrix;

	
	Matrix4x4 m_PixelsToModelMatrix;
	Matrix4x4 m_ModelToPixelsMatrix;

	Matrix4x4 m_ModelToVolumeMatrix;
	Matrix4x4 m_PixelsToViewMatrix;



	Matrix4x4 m_VolumeToModelMatrix;

	Matrix4x4 m_model_view;
	float *m_Correction;
	
	Vector4 m_CameraPos;

	
	int m_ImageInUseSize[2];
	int m_ImageOrigin[2];

	
	float *m_ImageData;
	float *m_WeightOfImage;


	float *m_EachVoxel_TotalWeight;
	float *m_Voxel_Contributions;

	float *m_imageraylength;
	bool  m_BackProjection;
	
	double m_offx;
	double m_offy;
	

};

#endif