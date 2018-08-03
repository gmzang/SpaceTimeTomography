//ST_Tomography.h

#ifndef __ST_Tomography_h
#define __ST_Tomography_h


#define cimg_use_tiff
#define cimg_use_tif
#include "CImg.h"
#include "Scene.h"
#include <string>

class VolumeData;
class RayCaster;
class TypicalRayProducer;
class CoherentVolume;
//class GridPartition;



class ST_Tomography : public Scene
{
public:

	

	ST_Tomography(double d1, double d2, double d3, double d4, double d5, double d6, const char *outputfile, int *imgwh, int projNo, int iters, int sartiters,
		double alp, double sid, const char **projsfile, double *oxyz,  double sdd, double ds, const char * _prior, 
		const char * _bp, int *volxyz,int nframes,double *startdg, int nrounds);
	virtual ~ST_Tomography();
		void InitViews();
	bool RotateVolume(int i );
	// data access 
	void SetData(VolumeData *data);
	void SART(int iframes);
	
	void SetBackGroundColor(float r,float g, float b);
	void GetBackGroundColor(float& r,float& g,float& b);

	void SetIsovalue(float isovalue);
	float GetIsovalue();

	void SetSurfaceColor(float r,float g,float b);
	void GetSurfaceColor(float &r,float &g,float &b);

	void SetLightDirection(float x,float y,float z);
	void SetLightIntensity(float intensity);
	void SetLightColor(float r,float g,float b);

	void SetAmbient(float value);
	float GetAmbient();

	void SetDiffuse(float value);
	float GetDiffuse();

	void SetSpecular(float value);
	float GetSpecular();

	void SetSpecularPower(float value);
	float GetSpecularPower();
	
	void  SetSampleDistance(float sd);
	float GetSampleDistance();

	
	virtual void Init();
	virtual void Resize(int width, int height);
	//void addMhaHeader(const int volDims[3],const float volSpacings[3]);
protected:
	
	virtual void RenderContent(int iframe);
	virtual void GetSize3D(float size[3]);
	virtual void Render();
	
private:
	// main loop
	bool _renderCoreFuc(int iframe);
	//bool _renderCoreFuc();

	// data
	VolumeData* m_Data;
	

	// ray caster
	RayCaster* m_RayCaster;
	// ray generator
	TypicalRayProducer* m_RayProducer;
	// main direction
	int m_mainDir;

	float m_BackgroundColor[3];
	//float *VoxelValue;
	float **pImg;

	float *pProImg;
	 //assign value to correction
	//Quaternion *pQuaternion;
	Quaternion **pQuaternion;
	int m_nframes;
	int m_nrounds;


	double m_s;
	double m_t;
	double m_l;
	const char *pOutput;
	const char **pInput;

	double m_mu;
	float **m_startDegree;
	float m_huberfactor;

};

#endif