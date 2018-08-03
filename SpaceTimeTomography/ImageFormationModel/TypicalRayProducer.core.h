// TypicalRayProducer.core.h


#include "TypicalRayProducer.h"
//#include "ImageIO.h"
#include <iostream>

using namespace std;

//static inline Vector4 CastRay(const Vector4& voxel_zero,const Vector4& voxel_unit,const Vector4& CameraDirection,float tnear,float tfar,float density);

//static inline Vector4 CastingRay(const Vector4& voxel_zero,const Vector4& voxel_unit,const Vector4& CameraDirection,float tnear,float tfar);
static inline double  CastingRay(const Vector4& voxel_zero,const Vector4& voxel_unit,const Vector4& CameraDirection,float tnear,float tfar,float &SumWeight_Ray, int indexofImage);
//static inline Vector4 CastRay(const Vector4& voxel_zero,const Vector4& voxel_unit,const Vector4& CameraDirection,float tnear,float tfar,float &density)

static Matrix4x4 P2MMatrix;
static Matrix4x4 M2VMatrix;
static int Volume_Dims[3];

static float CroppingBounds[6];	
static Vector4 EyePos;

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

// calculate original for source, directions for ray, camera, and distance (near and far) from position of the measurement from detector
static inline void ProduceRay(float pixX,float pixY, Vector4& voxel_zero,Vector4& voxel_unit, Vector4& CameraDirection, float& tnear, float& tfar)
{	



	Vector4 model_pix=P2MMatrix*Vector4(pixX,pixY,0.0f,1.0f);
	model_pix*=1.0f/model_pix[3];

	CameraDirection=model_pix-EyePos;
	CameraDirection.Normalize();

	voxel_zero=M2VMatrix*EyePos;
	voxel_unit=M2VMatrix*CameraDirection;

	/// Find Start Point
	float t1,t2;
	//cout<<"CroppingBounds: "<<CroppingBounds[0]<<" "<<CroppingBounds[1]<<" "<<CroppingBounds[2]<<" "<<CroppingBounds[3]<<" "<<CroppingBounds[4]<<" "<<CroppingBounds[5]<<" "<<endl;
	t1=(CroppingBounds[0]-voxel_zero[0])/voxel_unit[0];
	t2=(CroppingBounds[1]-voxel_zero[0])/voxel_unit[0];

	float tminX=min(t1,t2);
	float tmaxX=max(t1,t2);

	t1=(CroppingBounds[2]-voxel_zero[1])/voxel_unit[1];
	t2=(CroppingBounds[3]-voxel_zero[1])/voxel_unit[1];

	float tminY=min(t1,t2);
	float tmaxY=max(t1,t2);

	t1=(CroppingBounds[4]-voxel_zero[2])/voxel_unit[2];
	t2=(CroppingBounds[5]-voxel_zero[2])/voxel_unit[2];

	float tminZ=min(t1,t2);
	float tmaxZ=max(t1,t2);

	tnear=max(max(tminX, tminY), tminZ);
	tnear=max(tnear,0.0f);
	tfar=min(min(tmaxX,tmaxY),tmaxZ);

//
//	//change about the world coordnation of image
//	Vector4 model_pix=P2MMatrix*Vector4(pixX,pixY,0.0f);
//	//Vector4 model_pix=P2MMatrix*Vector4(pixX,pixY,-345.924f);
//	model_pix*=1.0f/model_pix[3];
//
//	/*if(pixX==0&&pixY==0) cout<<"		0 0 model_pix"<<model_pix[0]-64<<" "<<model_pix[1]-64<<" "<<model_pix[2]-64<<" "<<model_pix[3]<<" "<<endl;
//	if(pixX==64&&pixY==64) cout<<"		64 64 model_pix"<<model_pix[0]-64<<" "<<model_pix[1]-64<<" "<<model_pix[2]-64<<" "<<model_pix[3]<<" "<<endl;
//	if(pixX==128&&pixY==128) cout<<"	128 128 model_pix"<<model_pix[0]-64<<" "<<model_pix[1]-64<<" "<<model_pix[2]-64<<" "<<model_pix[3]<<" "<<endl;
//*/
//
//
//	CameraDirection=model_pix-EyePos;
//	CameraDirection.Normalize();
//	/*if(pixX==0&&pixY==0)
//		cout<<"		0 0 Eyepos:"<<EyePos[0]-64<<"	"<<EyePos[1]-64<<" "<<EyePos[2]-64<<" "<<EyePos[3]<<endl;
//	if(pixX==64&&pixY==64)
//		cout<<"		64 64 Eyepos:"<<EyePos[0]-64<<"	"<<EyePos[1]-64<<" "<<EyePos[2]-64<<" "<<EyePos[3]<<endl;
//	if(pixX==128&&pixY==128)
//		cout<<"		128 128 Eyepos:"<<EyePos[0]-64<<"	"<<EyePos[1]-64<<" "<<EyePos[2]-64<<" "<<EyePos[3]<<endl;*/
//	voxel_zero=M2VMatrix*EyePos;
//	voxel_unit=M2VMatrix*CameraDirection;
//	//cout<<voxel_unit[0]<<"   "<<voxel_unit[1]<<"   "<<voxel_unit[2]<<"   "<<voxel_unit[3]<<endl;
//	float t1,t2;
////	cout<<"voxel_zero"<<voxel_zero[0]<<" "<<voxel_zero[1]<<" "<<voxel_zero[2]<<endl;
//	t1=-voxel_zero[0]/voxel_unit[0];
//	t2=(Volume_Dims[0]-1.0f-voxel_zero[0])/voxel_unit[0];
//
//	float tminX=min(t1,t2);
//	float tmaxX=max(t1,t2);
//
//	t1=-voxel_zero[1]/voxel_unit[1];
//	t2=(Volume_Dims[1]-1.0f-voxel_zero[1])/voxel_unit[1];
//	
//	float tminY=min(t1,t2);
//	float tmaxY=max(t1,t2);
//
//	t1=-voxel_zero[2]/voxel_unit[2];
//	t2=(Volume_Dims[2]-1.0f-voxel_zero[2])/voxel_unit[2];
//
//	float tminZ=min(t1,t2);
//	float tmaxZ=max(t1,t2);
//	//cout<<t1<<"   "<<t2<<endl;//"   "<<voxel_unit[2]<<"   "<<voxel_unit[3]<<endl;
//	tnear=max(max(tminX, tminY), max(tminX, tminZ));
//	tnear=max(tnear,0.0f);
//	tfar=min(min(tmaxX,tmaxY),min(tmaxX,tmaxZ));
//
//

}


void GenerateImage_MT(void* userData )
{
	TypicalRayProducer *RayProducer=(TypicalRayProducer *)userData;

	int imageSize[2];
	int imageOrigin[2];
	//int detectorSize[2];
//	RayProducer->GetDetectorSize(detectorSize);
	RayProducer->GetImageInUseSize(imageSize);
	RayProducer->GetImageOrigin(imageOrigin);
	float* image=RayProducer->GetImageData();//ZGM    get the data
	float *weight=RayProducer->GetWeightOfImage();
	float *imgLen=RayProducer->GetImageRaylength();
	if (!image) return;

	//int i,j;

	//int boundx1 = 30;
	//int boundx2 = 30;
	//int boundy1 = 200;
	//int boundy2 = 2;

	/*int boundx1 = 700;
	int boundx2 = 50;
	int boundy1 = 300;
	int boundy2 = 300;*/

	//int boundx1 = 60;
	//int boundx2 = 60;
	//int boundy1 = 20;
	//int boundy2 = 850;
	//for (j= boundy1;j<imageSize[1]- boundy2;j++)
	//{
	//	#pragma  omp parallel for
	//	for (i= boundx1;i<imageSize[0]- boundx2;i++)
	//int boundx1 = 60;
	//int boundx2 = 60;
	//int boundy1 = 20;
	//int boundy2 = 314;
	int boundx1 = 5;
	int boundx2 = 5;
	int boundy1 = 5;
	int boundy2 = 5;

// should be test!!! for openmp!!!
#pragma  omp parallel for
for (int j = 0;j<imageSize[1] ;j++)
	{
		
		for (int i = 0;i<imageSize[0] ;i++)
		{


			int m_imageindex = i + j*imageSize[0];
			if(i<boundx1||i>imageSize[0]- boundx2||j<boundy1|| j>imageSize[1]-boundy2)
			{
				RayProducer->m_WeightOfImage[m_imageindex] = (float)0.0f; //ray weight

																				   //float *pImage=image+(i-imageOrigin[0]+(j-imageOrigin[1])*imageSize[0]); 

				image[m_imageindex] = 0.0f;
			}


			

			Vector4 voxel_zero, voxel_unit, CameraDirection;
			float tnear,tfar;
		
			
		    //ProduceRay(i+0.5f,j+0.5f,voxel_zero,voxel_unit,CameraDirection,tnear,tfar);
			//ProduceRay(i+0.5,j+0.5,voxel_zero,voxel_unit,CameraDirection,tnear,tfar);
			ProduceRay(i,j,voxel_zero,voxel_unit,CameraDirection,tnear,tfar);

			//Vector4 value(0.0f,0.0f,0.0f,0.0f);
			double value=0.0f;
			float  SumWeight_Ray=0.0f;
			
			//cout<<imageSize[0]<<"	"<<imageSize[1]<<endl;
			//if(m_imageindex>=imageSize[1]*imageSize[0]||m_imageindex<0.0)
			//	cout<<m_imageindex<<endl;
			
			if (tnear<tfar) 
			{
					RayProducer->GetImageRaylength()[m_imageindex]=tfar-tnear;
			
					value= CastingRay(voxel_zero,voxel_unit,CameraDirection,tnear,tfar,SumWeight_Ray,m_imageindex); 
			}
		//	cout<<"SumWeight"<<SumWeight_Ray<<endl;
			RayProducer->m_WeightOfImage[m_imageindex]=(float)SumWeight_Ray; //ray weight
		
			//float *pImage=image+(i-imageOrigin[0]+(j-imageOrigin[1])*imageSize[0]); 


			image[m_imageindex] = value;
				//image[m_imageindex]=max(min(value,10000000.0f),-10000000.0f);
			//	if(image[m_imageindex]>1000000||image[m_imageindex]<-1000000) cout<<"Something wrong for the voxel valuse" <<image[m_imageindex]<<endl;
			

		
		}
	}
}



//  generating image
static inline void GenerateImage_GetCorrections(TypicalRayProducer *RayProducer)
{
	
	if (!RayProducer) return;
	RayProducer->UpdateMatrixAndBounds();
	P2MMatrix=RayProducer->GetPixelsToModelMatrix();
	M2VMatrix=RayProducer->GetModelToVolumeMatrix();
	float Volume_Spacings[3];
	RayProducer->GetVolSize(Volume_Dims,Volume_Spacings);


	CroppingBounds[0]=CroppingBounds[2]=CroppingBounds[4]=0.0f;
	CroppingBounds[1]=(float)(Volume_Dims[0]-1);
	CroppingBounds[3]=(float)(Volume_Dims[1]-1);
	CroppingBounds[5]=(float)(Volume_Dims[2]-1);
	//CroppingBounds[1] = (float)(Volume_Dims[0] + 4);
	//CroppingBounds[3] = (float)(Volume_Dims[1] + 4);
	//CroppingBounds[5] = (float)(Volume_Dims[2] + 4);
	//CroppingBounds[2] = 5;
	//CroppingBounds[3] = (float)(Volume_Dims[1] - 5);

	/*CroppingBounds[0] = CroppingBounds[2] = CroppingBounds[4] = 5.0f;
	CroppingBounds[1] = (float)(Volume_Dims[0] - 20);
	CroppingBounds[3] = (float)(Volume_Dims[1] - 20);
	CroppingBounds[5] = (float)(Volume_Dims[2] - 20);*/
	EyePos[0]=P2MMatrix[8]/P2MMatrix[11];
	EyePos[1]=P2MMatrix[9]/P2MMatrix[11];
	EyePos[2]=P2MMatrix[10]/P2MMatrix[11];
	EyePos[3]=1.0f;

	//LaunchSPMD(GenerateImage_MT,RayProducer);

	GenerateImage_MT(RayProducer);
	//TypicalRayProducer *RayProducer=(TypicalRayProducer *)userData;




}
