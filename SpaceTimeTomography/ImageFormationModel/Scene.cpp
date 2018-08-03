// Scene.cpp
#include "Scene.h"
#include "GL/glut.h"
#include <iostream>
using namespace std;
static const float pi=3.1415926536f;




int mmm=0;
Scene::Scene()
{


	//m_Fovy=10.0;
	//m_Fovy=20.0f;
	m_Near=1.0f;
	m_Far=5.0f;
	//m_Far=1.002f;
	m_Translate[0]=m_Translate[1]=m_Translate[2]=0.0f;
	//m_Translate[1]=-
	m_Scale=1.0f;

	//m_NeedCalculateCamera=true;
	//m_NeedCalculateMatrix=true;


	m_NeedCalculateCamera=false;
	m_NeedCalculateMatrix=false;

//	m_RadPerPixel=0.01f;
//	m_RadPerPixel=0.005f;// should pay attention!~ 
//	m_RadPerPixel=0.008646f;  //1 DU  0.008646
	//m_RadPerPixel=0.0087250f;   //0.0087250 0.0087266453f
	m_RadPerPixel=0.0087266453f;
	m_ScaleRatePerPixel=logf(1.002f);

	m_Rotating=false;
	m_Scaling=false;
	m_Moving=false;
	m_ds=1.0f;
}

Scene::~Scene()
{

}

void Scene::SetFovy(float fovy)
{
	//m_Fovy=fovy;//Fovy
	//m_Fovy=Fovy; //replace it simply
	m_NeedCalculateCamera=true;
}

void Scene::Resize(int width, int height)
{
	cout<<"W in scene before:" <<width<<"	H:"<<height<<endl;
	glViewport(0, 0, width, height);
	m_ViewPort[0]=width;
	m_ViewPort[1]=height;
	cout<<"W in scene:" <<width<<"	H:"<<height<<endl;
	m_NeedCalculateMatrix=true;

	//cout<<"PI:"<<pi<<endl;
	//cout<<"FOVY: "<<Fovy<<endl;
}
void Scene::_calculateCamera()
{
	//mmm++;
	//cout<<mmm<<endl;
	/// 在这里初始化
	float size[3];
	GetSize3D(size);
	
		m_ZTranslateBase=-1*m_sid;  //800

	// same_fan_fan_ is good  projections need left offset
		//m_Translate[0] = 0.0f; //-2.0061f and  0.84015f is samehalf which is good


							   // if  y offsetthe xtex file is "-" . then here is "+" .
		//m_Translate[1] = 0.0f;; // 20.9767737621844f 12.1218347737581f 
		//									 ////m_Translate[0]=0.0f;

		//									 ////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
		//									 // if  z offsetthe xtex file is "-" . then here is "-" . possible
		//									 // from 20-->0. from top to down. i.e. downside.
		//m_Translate[2] =0.0f; //2.9935

		//m_Translate[1] = -11.1644838684374f; // 20.9767737621844f 12.1218347737581f 
		//									 ////m_Translate[0]=0.0f;

		//									 ////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
		//									 // if  z offsetthe xtex file is "-" . then here is "-" . possible
		//									 // from 20-->0. from top to down. i.e. downside.
		//m_Translate[2] = -2.07001068685314f; //2.9935
		//m_Translate[1] = 6.38405943477408f; // 20.9767737621844f 12.1218347737581f 
		//									//									 ////m_Translate[0]=0.0f;

		//									//									 ////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
		//									//									 // if  z offsetthe xtex file is "-" . then here is "-" . possible
		//									//									 // from 20-->0. from top to down. i.e. downside.
		//m_Translate[2] = 2.62873035549521f; //2.9935

		// x offset xtect file is "-", then "-" here
	//m_Translate[0]= 0.0f; //-2.0061f and  0.84015f is samehalf which is good


	//// if  y offsetthe xtex file is "-" . then here is "+" .
	//m_Translate[1] = -3.93865108171589f; // 20.9767737621844f 12.1218347737581f 
	//									////m_Translate[0]=0.0f;

	//									////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
 //  // if  z offsetthe xtex file is "-" . then here is "-" . possible
	//// from 20-->0. from top to down. i.e. downside.
	//m_Translate[2] =  -5.00771351818163f; //2.9935


	//// if the xtex file is - . then here is + .
	//m_Translate[1] = 12.1218347737581f; //   dough  6.04102757758221f
	//									////m_Translate[0]=0.0f;
	//									////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
	//m_Translate[2] = -0.528369758972204f; //2.9935  dough -0.528369758972204f
	////m_Translate[1] = 0.0f; //12.1218347737581f 
	//									////m_Translate[0]=0.0f;

										////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
	//m_Translate[2] = 0.0f; //2.9935
	m_Near=m_sdd;
	//m_Far=(-m_ZTranslateBase+r)*5.0;
	m_Far=10000.0;

	m_NeedCalculateMatrix=true;
	

}

void Scene::SetSize3DModified()
{
	// x offset xtect file is "-", then "-" here
	//m_Translate[0] =0.0f; //-2.0061f and  0.84015f is samehalf which is good
						  // if  y offsetthe xtex file is "-" . then here is "+" .
	//m_Translate[1] = 6.38405943477408f; // 20.9767737621844f 12.1218347737581f 
	//									//									 ////m_Translate[0]=0.0f;

	//									//									 ////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
	//									//									 // if  z offsetthe xtex file is "-" . then here is "-" . possible
	//									//									 // from 20-->0. from top to down. i.e. downside.
	//m_Translate[2] = 2.62873035549521f; //2.9935
	////									// if  y offsetthe xtex file is "-" . then here is "+" .
	//m_Translate[1] = -3.93865108171589f; // 20.9767737621844f 12.1218347737581f 
	//									 ////m_Translate[0]=0.0f;

	//m_Translate[1] = -11.1644838684374f; // 20.9767737621844f 12.1218347737581f 
	//									 ////m_Translate[0]=0.0f;

	//									 ////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
	//									 // if  z offsetthe xtex file is "-" . then here is "-" . possible
	//									 // from 20-->0. from top to down. i.e. downside.
	//m_Translate[2] = -2.07001068685314f; //2.9935

	//									 ////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
	//									 // if  z offsetthe xtex file is "-" . then here is "-" . possible
	//									 // from 20-->0. from top to down. i.e. downside.
	//m_Translate[2] = -5.00771351818163f; //2.9935
	//m_Translate[0] = 0.0f;//m_Translate[1]=m_Translate[2]=0.0f;  //要不要平移
	//m_Translate[1] = 13.5578991445248f; //12.1218347737581f 
	//									////m_Translate[0]=0.0f;

	//									////m_Translate[1]=0.0f; //往下平移了一些。1.0744f  ,3.078
	//m_Translate[2] = 9.25380417800901f; //2.9935
	//m_Translate[1] = 0.0f;
	//m_Translate[2] = 0.0f;
	m_Scale=1.0f;
	m_Rotation.Identity();
	m_NeedCalculateCamera=true;
}

void Scene::_calculateMatrix()
{
	//Transpose
	//m_ViewMatrix=Matrix4x4::TranslateMatrix(m_Translate[0],m_Translate[1],m_ZTranslateBase+m_Translate[2])
	//	*Matrix4x4::ScaleMatrix(m_Scale,m_Scale,m_Scale)*m_Rotation.ToMatrix();

	m_ViewMatrix=Matrix4x4::TranslateMatrix(m_Translate[0],m_Translate[1],m_ZTranslateBase+m_Translate[2])
		*Matrix4x4::ScaleMatrix(m_Scale,m_Scale,m_Scale)*m_Rotation.ToMatrix();
	//m_ViewMatrix=m_Rotation.ToMatrix()*Matrix4x4::ScaleMatrix(m_Scale,m_Scale,m_Scale)*Matrix4x4::TranslateMatrix(m_Translate[0],m_Translate[1],m_ZTranslateBase+m_Translate[2]);



	
	float *projMat=m_ProjectionMatrix.ele;
	projMat[0]=m_Near/m_halfFrustumWidth;

	projMat[1]=projMat[2]=projMat[3]=0.0f;
	projMat[4]=0.0f;
	projMat[5]=m_Near/m_halfFrustumHeight;

	projMat[6]=projMat[7]=projMat[8]=projMat[9]=0.0f;
	projMat[10]=-(m_Far+m_Near)/(m_Far-m_Near);
	projMat[11]=-1.0f;
	projMat[12]=projMat[13]=0.0f;
	projMat[14]=-2.0f*m_Far*m_Near/(m_Far-m_Near);
	projMat[15]=0.0f;
}


void Scene::RenderContent(int iframe)
{
	//glClearColor(0.5f,0.5f,0.5f,1.0f);
	//glClearDepth(1.0f);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LEQUAL);
	//glBegin(GL_QUADS);

	//glEnd();

}


void Scene::GetSize3D(float size3D[3])
{
	size3D[0]=size3D[1]=size3D[2]=1.0f;
}

const Matrix4x4& Scene::GetViewMatrix()
{
	return m_ViewMatrix;
}
const Matrix4x4& Scene::GetProjectionMatrix()
{
	return m_ProjectionMatrix;
}

void Scene::_makeReference()
{
	m_OldTranslate[0]=m_Translate[0]; m_OldTranslate[1]=m_Translate[1]; m_OldTranslate[2]=m_Translate[2];
	m_OldScale=m_Scale;
	m_OldRotation=m_Rotation;
}




void Scene::GetTranslate(float translate[3])
{
	translate[0]=m_Translate[0];
	translate[1]=m_Translate[1];
	translate[2]=m_Translate[2];
}

void Scene::SetTranslate(const float translate[3])
{
	m_Translate[0]=translate[0];
	m_Translate[1]=translate[1];
	m_Translate[2]=translate[2];
	m_NeedCalculateMatrix=true;
}

void Scene::SetScale(float scale)
{
	m_Scale=scale;
	m_NeedCalculateMatrix=true;
}
void Scene::SetRotation(const Quaternion& rotation)
{
	m_Rotation=rotation;
	m_NeedCalculateMatrix=true;
}