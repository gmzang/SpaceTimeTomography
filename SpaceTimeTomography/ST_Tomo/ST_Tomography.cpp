//ST_Tomography.cpp

#include "ST_Tomography.h"
#include "RayCaster.h"
#include "TypicalRayProducer.h"
#include "VolumeData.h"
#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <string>
#include "GL/glew.h"
#include <malloc.h>
#include <fstream>
#include <time.h>
#include <math.h>
#include <sstream>
#include "tiffio.h"
#include "tiff.h"

#include"optical_flow.h"
#include"warping.h"

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <string>
#include<Eigen/Core>
#include <Eigen/Dense>
#include<Eigen/SVD>
#include <Eigen/Householder>
#include <time.h> 

#include "CoherentVolume.h"
#include<fstream>
#include <map>
#include"warping.h"
using namespace cimg_library;

namespace cl = cimg_library;

typedef cl::CImg<float>		image;
typedef CImgList<float>  imagelist;
imagelist Projslist;
imagelist Y_tv;//(dims[0], dims[1], dims[2], 3, 0.0);
imagelist Y_tp;// temporal prior
imagelist Y_bc;// bright consistency prior

imagelist X_bar;//
imagelist X_prev;//
image  Y1;//
image  Y2;//
image  Y3;//
image Y4;//
image tim;
image imROI;

// for flow
image uvw;
image I1;
image I2;
image tmpimage;


using namespace std;
int myrandom(int i) { return std::rand() % i; }
#if !defined(Circular_Angle)
# define Circular_Angle 360
#endif

#if !defined(IS_COUNTERCLOCK)
#  define IS_COUNTERCLOCK 1 
#define Voxelbased 1
#define Raybased  2
#define ATV   3
#endif

static int PRO_WIDTH;
static int PRO_HEIGHT;
static int iminY;
static int PRO_NO = 20; //init
static int ITERATIONS = 1;//
static int CP_ITERS = 1;//

double m_theta = 1.0;
static const float pi = 3.1415926536f;

float m_para;
float G_Max = 10.0f;
double Alpha = 0.1;
double m_alpha = Alpha;
#define USING_PRIOR ATV;
int  Prior = ATV;
int  Backprojector = Raybased;
int  Joint_iter = 20;


// params for volume density reconstruction in CP algorithm
float  L = 8.0;
float CP_tau = 1; //0.3

float CP_sigma = .9 / (L*CP_tau);

float	m_aabb[6];
// previous 0.03 VS 1.0
float TV_w = 0.1; // TV weight
float TP_w = 0.5;  // temporal prior weight  0.1  0.5  and 0.2
float BC_w = 0.5;  // temporal prior weight  0.1  0.5  and 0.2

float SART_lambda = CP_tau / 2;
#  define Factor sqrtf(2*SART_lambda) //change  //finish chang 1

// subvolume/region of interest
int cut[] = { 8,8,8,45,8,8,30,20 };

int boundx1 = 5;
int boundx2 = 5;
int boundy1 = 5;
int boundy2 = 5;


static inline void _Cal_Norm(float &u, float sigma)
{

	if (u > sigma) u = sigma;
	else if (u < -sigma) u = -sigma;

}
int dims[3];
int timestorender = 0;
int rot = 1;
//static int _index_= 0;
static int _ReadyToRotate_ = 1;
//static double Angle_Step= (double)Circular_Angle/(double)(PRO_NO-1.0);

static double Angle_Step = (double)Circular_Angle / (double)(PRO_NO); //test if identical to RTK results

static float center[3] = { dims[0] / 2.0f,dims[1] / 2.0f,dims[2] / 2.0f };

static int projsEachRound = 10;
static double degreeEachProj = 360.0 / projsEachRound;
float bak_densityaverage = 0.0f;
//bool SetOriginalPos = false;
int index_of_project = 0;
char* pOutputFile;
const char* outputfilename;
int len;
//string output = "t";  //RTK_Origin_0  SL128_RTK_KB
//string outputformat = ".mha";
int numofSART = 0;
//int NO_SART=1;
int *idxSART;
int noindex = 0;
float Views = 0;
double sampledistance = 0.0f;
int vh_x1 = 0;
int vh_x2 = 0;
int vh_y1 = 0;  
int vh_y2 = 0;
int vh_z1 = 0;
int vh_z2 = 0;
int No_PO = 0;
image warpI2;


Eigen::MatrixXf pY;
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef idx
#define idx(i,j,twidth)            (i + j*twidth)
#endif

#ifndef idx3D
#define idx3D(i,j,k,twidth,theight)            (i + j*twidth+ k*twidth*theight)
#endif

float inline cubicInterpolate(float p[4], float x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

}


float inline bicubicInterpolate(float p[4][4], float x, float y) {
	float arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);

}


ST_Tomography::ST_Tomography(double d1, double d2, double d3, double d4, double d5, double d6, const char *outputfile, int *imgwh, int projno, int iters, int sartiters,
	double alp, double sid, const char **projsfile, double *oxyz, double sdd, double ds, const char * _prior, const char * _bp, int *volxyz,
	int nframes, double *startdg, int nrounds)
{

	m_nframes = nframes;
	m_nrounds = nrounds;

	m_RayProducer = new TypicalRayProducer(volxyz[0], volxyz[1], volxyz[2], imgwh[0], imgwh[1], ds, _bp);
	Alpha = alp;
	m_alpha = Alpha;

	PRO_WIDTH = imgwh[0];
	PRO_HEIGHT = imgwh[1];
	iminY = imgwh[0];
	TV_w = d3;
	TP_w = d4;
	BC_w = d5;

	m_huberfactor = 1.0f / (1.0f + CP_sigma*d6 / d3);

	if (strcmp(_bp, "Raybased") == 0)
	{
		Backprojector = Raybased;
	}
	else
		Backprojector = Voxelbased;


	m_halfFrustumWidth = (imgwh[0] - 1)*ds*0.5f;
	m_halfFrustumHeight = (imgwh[1] - 1)*ds*0.5f;

	PRO_NO = projno;
	Angle_Step = (double)Circular_Angle / (double)(PRO_NO);
	Views = Angle_Step;

	ITERATIONS = sartiters;
	CP_ITERS = iters;

	m_sid = sid;
	m_sdd = sdd;
	m_ds = ds;
	m_ViewPort[0] = imgwh[0];
	m_ViewPort[1] = imgwh[1];
	m_Translate[0] = -oxyz[0];
	m_Translate[2] = +oxyz[1];
	m_Translate[1] = -oxyz[2];
	pOutput = outputfile;

	m_Data = NULL;
	m_BackgroundColor[0] = m_BackgroundColor[1] = m_BackgroundColor[2] = 0.0f;

	m_RayCaster = new RayCaster;

	pQuaternion = new Quaternion*[nframes];
	m_startDegree = new float*[nframes];
	pInput = new const char*[nframes];
	idxSART = new int[nframes];
	ifstream myfile("lds.txt");

	int idxrounds = 0;
	for (int i = 0;i < nframes;i++)
	{

		m_startDegree[i] = new float[m_nrounds];
		for (size_t idxr = 0; idxr < m_nrounds; idxr++)
		{
			myfile >> m_startDegree[i][idxr];
		}
		//m_startDegree[i] = startdg[i];
		pInput[i] = new char[PRO_NO];
		// need check this initilization
		pInput[i] = projsfile[i];
		pQuaternion[i] = new Quaternion[PRO_NO];
		idxSART[i] = 1;
	}
	myfile.close();
	cout << PRO_WIDTH << " " << PRO_HEIGHT << " " << volxyz[0] << " " << volxyz[1] << " " << volxyz[2] << endl;
	pY = Eigen::MatrixXf::Zero(PRO_NO, PRO_WIDTH*PRO_HEIGHT);
}

ST_Tomography::~ST_Tomography()
{
	delete m_RayProducer;
	delete m_RayCaster;
}

void ST_Tomography::SetData(VolumeData *data)
{
	m_Data = data;
	SetSize3DModified();

	if (!data) return;
	float spacings[3];
	m_Data->GetDimensions(dims);
	m_Data->GetSpacings(spacings);

	m_RayProducer->SetVolSize(dims, spacings);

	// change the sample rate
	m_RayCaster->SetSampleDistance(sqrt(spacings[0] * spacings[0] + spacings[1] * spacings[1] + spacings[2] * spacings[2]));
	sampledistance = m_RayCaster->GetSampleDistance();

	imROI = image(dims[0], dims[1], dims[2], 1.0);

	#pragma omp parallel for
	for (int k = 0; k < dims[2]; k++)
		for (int j = 0; j < dims[1]; j++)
			for (int i = 0; i < dims[0]; i++)
			{
				{
					{

						imROI(i, j, k) = 1.0f;
						if (((i < cut[6]) && (k < cut[7])) || ((i < cut[6]) && (k > (dims[2] - cut[7]))))
							imROI(i, j, k) = 0.0f;
						if (((i > (dims[0] - cut[6])) && (k < cut[7])) || ((i > (dims[0] - cut[6])) && (k > (dims[2] - cut[7]))))
							imROI(i, j, k) = 0.0f;

						if (((k < cut[6]) && (i < cut[7])) || ((k < cut[6]) && (i > (dims[0] - cut[7]))))
							imROI(i, j, k) = 0.0f;
						if (((k > (dims[2] - cut[6])) && (i < cut[7])) || ((k > (dims[2] - cut[6])) && (i > (dims[0] - cut[7]))))
							imROI(i, j, k) = 0.0f;

						if (i<cut[0] || i>dims[0] - cut[1] || j<cut[2] || j>dims[1] - 7 || k<cut[4] || k>dims[2] - cut[5])
							imROI(i, j, k) = 0.0f;

					}
				}
			}
	m_aabb[0] = 0;
	m_aabb[1] = dims[0];
	m_aabb[2] = 0;
	m_aabb[3] = dims[1];
	m_aabb[4] = 0;
	m_aabb[5] = dims[2];

	int psize = dims[0] * dims[1] * dims[2];

	Y1 = image(dims[0], dims[1], dims[2], 3, 0.0);
	Y2 = image(dims[0], dims[1], dims[2], 3, 0.0);
	Y3 = image(dims[0], dims[1], dims[2], 3, 0.0);
	Y4 = image(dims[0], dims[1], dims[2], 1, 0.0);

	Y_tv = imagelist(m_nframes + 1, dims[0], dims[1], dims[2], 3, 0.0);
	// remove 4
	Y_tp = imagelist(m_nframes + 1, dims[0], dims[1], dims[2], 1, 0.0);
	Y_bc = imagelist(m_nframes + 1, dims[0], dims[1], dims[2], 1, 0.0);

	X_bar = imagelist(m_nframes + 1, dims[0], dims[1], dims[2], 1, 0);
	X_prev = imagelist(m_nframes + 1, dims[0], dims[1], dims[2], 1, 0);

	float t_min = 1000.0f;
	float t_max = 0.0f;

	Projslist.insert(image(pInput[0]).get_split('z', m_nframes));
	t_max = Projslist[0].max_min(t_min);
	G_Max = t_max / 4.0f;

	InitViews();

}

void ST_Tomography::SetBackGroundColor(float r, float g, float b)
{
	m_BackgroundColor[0] = r; m_BackgroundColor[1] = g; m_BackgroundColor[2] = b;
}

void ST_Tomography::GetBackGroundColor(float& r, float& g, float& b)
{
	r = m_BackgroundColor[0]; g = m_BackgroundColor[1]; b = m_BackgroundColor[2];
}



void ST_Tomography::Init()
{
	glewInit(); // initial glew
}

void ST_Tomography::Resize(int width, int height)
{
	Scene::Resize(width, height);
	cout << "width: " << width << " height: " << height << endl;
	m_RayProducer->SetViewPort(width, height);
}

void  ST_Tomography::SetSampleDistance(float sd)
{
	m_RayCaster->SetSampleDistance(sd);
}
float ST_Tomography::GetSampleDistance()
{
	return m_RayCaster->GetSampleDistance();
}


void ST_Tomography::Render()
{

	cout << "==============================================" << endl;
	cout << "Joint-Reconstruction-Framework starts " << endl;
	cout << "==============================================" << endl;

	int dimension[3];
	m_Data->GetDimensions(dimension);
	int volumesize = dimension[0] * dimension[1] * dimension[2];

	int timeofPO = 1;

	image m_forward_warp;
	image m_backward_warp;

	for (size_t out_iter = 0; out_iter < Joint_iter; out_iter++)
	{
	

		cout << "==============================================" << endl;
		cout << "Flow-Reconstruction-subOptimization starts " << endl;
		cout << "==============================================" << endl;

			for (int iframe = 0; iframe < m_nframes - 1; iframe++)
			{
				char out_dir[1024];
				char basename[1024];
				char flow_field[1024];

				sprintf(basename, "FlowOutput");
				float dt = 1.0;

				I1 = 0.0;
				I2 = 0.0f;
				uvw = 0.0f;

				I1 = m_Data->m_Volumelist[iframe];
				I2 = m_Data->m_Volumelist[iframe + 1];
				int dim[] = { dimension[0], dimension[1], dimension[2] };

				uvw = m_Data->m_Flowlist[iframe];//.get_resize(166, 124, 148, 3, 5) for testing;
				const float aabb[] = { 0,dim[0], 0,dim[1], 0,dim[2] };

				float eta = 0.65;  //downsample
				float sigma = 0.5;  //blur 
				float precision = 0.1;
				int _norm = 4; //  huber
				int level = 3;

				float smoothness = 1.2;
				float huber = 1.0f;
				int iters = 3000;
				Y1 = 0.0f;
				Y2 = 0.0f;
				Y3 = 0.0f;
				Y4 = 0.0f;
				int nb = 3;
				uvw = optical_flow::L1Huber_OF_MultiScale3D(dim, I1, I2, uvw, Y1, Y2, Y3, Y4, precision,
					smoothness, iters, _norm, nb, eta, sigma, level, huber, cut);

				m_Data->m_Flowlist[iframe] = uvw;


				dt = 1.0f;
				tmpimage = warping::warp(dim, aabb, uvw, I2, -dt);

				tmpimage.save_tiff(str_format("%s/%03d.%03d.%03d.frame%02d.iter%02d.tiff", basename,
					dim[0], dim[1], dim[2], iframe, out_iter).c_str());

				uvw.save_analyze(str_format("%s/%03d.%03d.%03d.frame%02d.iter%02d.field.hdr", basename,
					dim[0], dim[1], dim[2], iframe, out_iter).c_str());
			}


		cout << "==============================================" << endl;
		cout << "Image-Reconstruction-subOptimization starts " << endl;
		cout << "==============================================" << endl;
		if (out_iter > 0)
		{
			for (No_PO = 0;No_PO < 2;No_PO++)
			{
				cout << "=============" << timeofPO++ << "  times in Proximal Operation iteration ==============" << endl;
				cout << "======= Huber threshold ==============" << m_huberfactor << endl;
				cout << "======= weight for each priors ==============" << m_huberfactor << endl;

				int kkk = 0;
				int inY = dimension[0];
				int inZ = dimension[0] * dimension[1];
				int inX = 1;
				string m_str;
				char temp[30];
				sprintf(temp, "%d", No_PO);
				m_str = temp;
				m_str += ".txt";
				for (int iframe = 0; iframe < m_nframes; iframe++)
				{
					if ((out_iter == 1) && (No_PO == 0))
					{
						X_bar[iframe] = m_Data->m_Volumelist[iframe];
						//if (iframe < m_nframes - 1)
						X_bar[iframe + 1] = m_Data->m_Volumelist[iframe + 1];

					}

					float dt = 1.0f;
					m_backward_warp = warping::warp(dims, m_aabb, m_Data->m_Flowlist[iframe], X_bar[iframe + 1], -dt);
					Y_tp[iframe] += CP_sigma*(X_bar[iframe + 1] - X_bar[iframe]);
					Y_bc[iframe] += CP_sigma*(m_backward_warp - X_bar[iframe]);

					//  the dual part

					#pragma omp parallel for
					for (int k = 1; k < dims[2] - 1; k++) {
						for (int j = 1; j < dims[1] - 1; j++) {
							for (int i = 1; i < dims[0] - 1; i++) {

								if (i<cut[0] || i>dims[0] - cut[1] || j<cut[2] || j>dims[1] - cut[3] || k<cut[4] || k>dims[2] - cut[5])
								{
									continue;
								}
								Y_tv[iframe](i, j, k, 0) += CP_sigma*(X_bar[iframe](i + 1, j, k) - X_bar[iframe](i, j, k))*m_huberfactor;
								Y_tv[iframe](i, j, k, 1) += CP_sigma*(X_bar[iframe](i, j + 1, k) - X_bar[iframe](i, j, k))*m_huberfactor;
								Y_tv[iframe](i, j, k, 2) += CP_sigma*(X_bar[iframe](i, j, k + 1) - X_bar[iframe](i, j, k))*m_huberfactor;


								Y_tv[iframe](i, j, k, 0) = max(-TV_w, min(TV_w, Y_tv[iframe](i, j, k, 0)));
								Y_tv[iframe](i, j, k, 1) = max(-TV_w, min(TV_w, Y_tv[iframe](i, j, k, 1)));
								Y_tv[iframe](i, j, k, 2) = max(-TV_w, min(TV_w, Y_tv[iframe](i, j, k, 2)));

								//for L2 norm in dual part
								Y_tp[iframe](i, j, k) = Y_tp[iframe](i, j, k)*TP_w / (TP_w + 2 * CP_sigma);

								Y_bc[iframe](i, j, k) = max(-BC_w, min(BC_w, Y_bc[iframe](i, j, k)));

							}
						}
					}

					X_prev[iframe] = m_Data->m_Volumelist[iframe];

					// adjoint operators
					// adjoint operator of backward warp and temporal prior
					if (iframe == 0)
					{
						m_Data->m_Volumelist[iframe] -= CP_tau*(-Y_tp[iframe]);
						m_Data->m_Volumelist[iframe] -= CP_tau*(-Y_bc[iframe]);

					}
					else
					{
						float dt = 1.0f;

						m_forward_warp = warping::warp(dims, m_aabb, m_Data->m_Flowlist[iframe - 1], Y_bc[iframe - 1], +dt);

						m_Data->m_Volumelist[iframe] -= CP_tau*(Y_tp[iframe - 1] - Y_tp[iframe]);
						m_Data->m_Volumelist[iframe] -= CP_tau*(m_forward_warp - Y_bc[iframe]);

					}


					// adjoint operator of TV

					#pragma omp parallel for
					for (int k = 1; k < dims[2] - 1; k++) {
						for (int j = 1; j < dims[1] - 1; j++) {
							for (int i = 1; i < dims[0] - 1; i++) {
								// for region of interest
								if (i<cut[0] || i>dims[0] - cut[1] || j<cut[2] || j>dims[1] - cut[3] || k<cut[4] || k>dims[2] - cut[5])
								{
									continue;
								}

								float y1 = -(Y_tv[iframe](i, j, k, 0) - Y_tv[iframe](i - 1, j, k, 0) + Y_tv[iframe](i, j, k, 1) -
									Y_tv[iframe](i, j - 1, k, 1) + Y_tv[iframe](i, j, k, 2) - Y_tv[iframe](i, j, k - 1, 2));

								m_Data->m_Volumelist[iframe](i, j, k) -= (CP_tau*y1);;

								m_Data->m_Volumelist[iframe](i, j, k) = max(min(G_Max, m_Data->m_Volumelist[iframe](i, j, k)), 0);
							}
						}
					}

					m_backward_warp = 0.0f;
					m_forward_warp = 0.0f;

					// SART as a solver for proximal operator of data term
					SART(iframe);

					#pragma omp parallel for
					for (int z = vh_z1;z < dimension[2] - vh_z2;z++)
					{
						for (int y = vh_y1;y < dimension[1] - vh_y2;y++)
						{
							for (int x = vh_x1;x < dimension[0] - vh_x2;x++)
							{
								if (x<cut[0] || x>dimension[0] - cut[1] || y<cut[2] || y>dimension[1] - cut[3] || z<cut[4] || z>dimension[2] - cut[5])
								{
									continue;
								}

								X_bar[iframe](x, y, z) = 2 * m_Data->m_Volumelist[iframe](x, y, z) - X_prev[iframe](x, y, z);

							}

						}

					}

					cout << "=========Image reconstruction=====-" << No_PO << "-CP iter " << iframe << "frame completed" << endl;
				}
			}  
		} 
	}// for outloop
	exit(0);
	return;

}

void ST_Tomography::SART(int iframes)
{
	pY.setZero();
	int kkk = 0;
	//////////////////////////////////////////////////////////////////////////
	// init volume data as U in  SART 
	//////////////////////////////////////////////////////////////////////////
	int dimension[3];
	//float spacings[3];
	m_Data->GetDimensions(dimension);
	int Yindex = dimension[0];
	int Zindex = dimension[0] * dimension[1];

	//cout<<"==============Begin the  SART algorithm=============="<<endl;
	time_t start, end;
	time(&start);

	for (int iter = 1;iter <= 1;iter++)
	{
		std::srand(unsigned(std::time(0)));
		std::vector<int> myvector;
		// set some values: 
		for (int i = 0; i < PRO_NO; ++i) myvector.push_back(i); // 1 2 3 4 5 6 7 8 9
		// using built-in random generator:
		std::random_shuffle(myvector.begin(), myvector.end());
		// using myrandom:
		std::random_shuffle(myvector.begin(), myvector.end(), myrandom);

		for (int i = 1;i <= PRO_NO;i++)
		{
			m_NeedCalculateMatrix = true;
			if (m_NeedCalculateMatrix)
			{
				this->m_Rotation = pQuaternion[iframes][myvector[i - 1]];
				cout << "i: " << i << " : " << myvector[i - 1] << "	 " << endl;;
				_calculateMatrix();
				m_NeedCalculateMatrix = false;
			}
			glMatrixMode(GL_PROJECTION);
			glLoadMatrixf(m_ProjectionMatrix);

			glMatrixMode(GL_MODELVIEW);
			glLoadMatrixf(m_ViewMatrix);

			////////////////////////////////////////////////////////////////////////////////////////////////////		
			////==================================ForwardProjection=======================================//////
			////////////////////////////////////////////////////////////////////////////////////////////////////
			RenderContent(iframes);
			////////////////////////////////////////////////////////////////////////////////////////////////////	
			////==================================Corrections==============================================//////
			////////////////////////////////////////////////////////////////////////////////////////////////////			
			#pragma omp parallel for
			for (int n = 0;n < PRO_HEIGHT;n++)
			{
				for (int m = 0;m < PRO_WIDTH;m++)
				{
					int k = idx(m, n, PRO_WIDTH);
					m_RayProducer->m_Correction[k] = (float)(Factor*Projslist[iframes](m, n, myvector[i - 1]) - Factor*m_RayProducer->GetImageData()[k] - pY(i - 1, k)) / (float)(Factor*m_RayProducer->GetImageRaylength()[k] + 1.0);

				}
			}
			//update pY for ST_Tomography algorithm
			#pragma omp parallel for
			for (int k = 0;k < PRO_WIDTH*PRO_HEIGHT;k++)
			{
				pY(i - 1, k) += m_alpha*m_RayProducer->m_Correction[k];
			}


			////backprojection

			if (Backprojector == Voxelbased)
			{

				#pragma omp parallel for
				for (int z = vh_z1;z < dimension[2] - vh_z2;z++)
				{
					for (int y = vh_y1;y < dimension[1] - vh_y2;y++)
					{
						for (int x = vh_x1;x < dimension[0] - vh_x2;x++)
						{
							int pos = idx3D(x, y, z, dimension[0], dimension[1]);
							Vector4 vec(x, y, z);
							Vector4 det_pix = m_RayProducer->m_model_view*m_RayProducer->m_VolumeToModelMatrix*vec;
							float _px = m_halfFrustumWidth - (det_pix[0] / det_pix[2])*m_sdd;
							float _py = m_halfFrustumHeight - (det_pix[1] / det_pix[2])*m_sdd;
							_px = _px / m_ds;
							_py = _py / m_ds;
							int px = (int)_px;
							int py = (int)_py;
							float dx = _px - px;
							float dy = _py - py;
							////bicubic interpolation
							if (px < 2 || px >= PRO_WIDTH - 2 || py < 2 || py >= PRO_HEIGHT - 2)
								continue;

							float *pC = m_RayProducer->m_Correction;
							pC += py*PRO_WIDTH + px;

							float pimg[4][4] = { pC[-1 - iminY],pC[-1],pC[-1 + iminY],pC[-1 + 2 * iminY],

								pC[-iminY],pC[0],pC[iminY],pC[2 * iminY],

								pC[1 - iminY],pC[1],pC[1 + iminY],pC[1 + 2 * iminY],

								pC[2 - iminY],pC[2],pC[2 + iminY],pC[2 + 2 * iminY]

							};

							float gray_value = bicubicInterpolate(pimg, dx, dy);
							m_Data->m_Volumelist[iframes](x, y, z) += m_alpha*gray_value;
							m_Data->m_Volumelist[iframes](x, y, z) = max(min(G_Max, m_Data->m_Volumelist[iframes](x, y, z)), 0);
						}
					}
				}
			}

			else   // for Raybased-backprojection
			{
				m_Data->SetProjectionMode(1);
				////////////////////////////////////////////////////////////////////////////////////////////////////
				////==================================BackProjection==============================================///
				////////////////////////////////////////////////////////////////////////////////////////////////////
				RenderContent(iframes);
				int dimens[3];
				m_Data->GetDimensions(dimens);

				int numberofMinus = 0;
				int nb_outliner = 0;

				#pragma omp parallel for
				for (int z = vh_z1;z < dimension[2] - vh_z2;z++)
				{
					for (int y = vh_y1;y < dimension[1] - vh_y2;y++)
					{
						for (int x = vh_x1;x < dimension[0] - vh_x2;x++)
						{
							int _index = idx3D(x, y, z, dimension[0], dimension[1]);
							float m_weight = m_RayProducer->m_EachVoxel_TotalWeight[_index];
							float m_contrb = (1.0*m_RayProducer->m_Voxel_Contributions[_index]) / m_weight;
							m_Data->m_Volumelist[iframes](x, y, z) += m_contrb;
							m_Data->m_Volumelist[iframes](x, y, z) = max(min(G_Max, m_Data->m_Volumelist[iframes](x, y, z)), 0);

						}

					}  //for 2
				} //for 3
				m_RayProducer->Initial();
				m_Data->SetProjectionMode(0);
			}

		}
	}



	int dimenst[3];
	m_Data->GetDimensions(dimenst);

	char Mybuf[10];
	sprintf(Mybuf, "%d", idxSART[iframes]);

	char Mybuf3[10];
	sprintf(Mybuf3, "%03d", iframes);

	string Myb = Mybuf;
	string Myb3 = Mybuf3;

	string _output = pOutput;

	string Mywritefile2 = _output + "-" + Myb3 + "-frame-" + Myb + "iter.tiff";
	cout << Mywritefile2 << endl;
	char* MyOutputfile2 = strdup(Mywritefile2.c_str());

	int dims2[3];
	float spacs2[3];

#pragma omp parallel for
	for (int z = 0;z < dimension[2] - 0;z++)
	{
		for (int y = 0;y < dimension[1] - 0;y++)
		{
			for (int x = 0;x < dimension[0] - 0;x++)
			{

				m_Data->m_Volumelist[iframes](x, y, z) *= imROI(x, y, z);
			}
		}
	}
	m_Data->m_Volumelist[iframes].save_tiff(MyOutputfile2);
	cout << "end---" << idxSART[iframes] << "th iteration for " << iframes << " frame" << endl;
	idxSART[iframes]++;
}


//settings for OPENGL
void ST_Tomography::RenderContent(int frame)
{
	glClearColor(m_BackgroundColor[0], m_BackgroundColor[1], m_BackgroundColor[2], 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	m_RayProducer->SetViewMatrix(this->GetViewMatrix());
	m_RayProducer->SetProjectionMatrix(this->GetProjectionMatrix());
	if (_renderCoreFuc(frame)) m_RayProducer->RenderImage(); //attention!

}


void ST_Tomography::GetSize3D(float size[3])
{
	int dims[3];
	float spacings[3];
	if (m_Data)
	{
		m_Data->GetDimensions(dims);
		m_Data->GetSpacings(spacings);
	}

	size[0] = dims[0] * spacings[0];
	size[1] = dims[1] * spacings[1];
	size[2] = dims[2] * spacings[2];
}

#include "SingleScalarVolume.core.h"
#include "RayCaster.core.h"
#include "TypicalRayProducer.core.h"
#include <iostream>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string>
using namespace std;


static VolumeData *vd;// 

bool ST_Tomography::_renderCoreFuc(int frame)
{

	if (!SingleScalarVolumePrepare(m_Data, m_RayProducer, frame)) return false;
	if (!RayCasterPrepare(m_RayCaster, dims)) return false;
	GenerateImage_GetCorrections(m_RayProducer);
	return true;

}
//



bool ST_Tomography::RotateVolume(int i)
{
	this->m_NeedCalculateMatrix = true;
	this->m_NeedCalculateCamera = true;
	this->_makeReference();
	double _Index = m_RadPerPixel;
	double angle = degreeEachProj*3.1415926 / 180;
	Vector4 v;
	v[0] = 0;
	v[1] = 1;
	v[2] = v[3] = 0.0f;
	if (!IS_COUNTERCLOCK)
		v[1] *= -1.0f;
	double length = v.Length();
	v *= 1.0f / length;
	this->m_Rotation = Quaternion(v[0], v[1], v[2], length*angle)*this->m_OldRotation;
	pQuaternion[i][index_of_project] = this->m_Rotation;
	index_of_project++;
	return true;
}


void ST_Tomography::InitViews()
{
	for (int i = 0; i < m_nframes; i++)
	{
		index_of_project = 0;
		for (size_t j = 0; j < m_nrounds; j++)
		{
			this->m_Rotation.Identity();
			m_OldRotation = m_Rotation;
			this->m_NeedCalculateMatrix = true;
			this->m_NeedCalculateCamera = true;
			this->_makeReference();
			double angle = this->m_startDegree[i][j] * 3.1415926 / 180;
			Vector4 v;
			v[0] = 0;
			v[1] = 1;
			v[2] = v[3] = 0.0f;
			if (!IS_COUNTERCLOCK)
				v[1] *= -1.0f;

			double length = v.Length();
			v *= 1.0f / length;
			//length*
			this->m_Rotation = Quaternion(v[0], v[1], v[2], length*angle)*this->m_OldRotation;
			pQuaternion[i][index_of_project] = this->m_Rotation;
			index_of_project++;
			for (int index = 0;index < projsEachRound - 1;index++)
			{
				if (m_NeedCalculateCamera)
				{
					_calculateCamera();
					m_NeedCalculateCamera = false;
				}
				m_NeedCalculateMatrix = true;
				if (m_NeedCalculateMatrix)
				{
					_calculateMatrix();
					m_NeedCalculateMatrix = false;
				}
				glMatrixMode(GL_PROJECTION);
				glLoadMatrixf(m_ProjectionMatrix);
				glMatrixMode(GL_MODELVIEW);
				glLoadMatrixf(m_ViewMatrix);
				RotateVolume(i);
			}
		}
	}
}



