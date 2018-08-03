// Scene.h
#ifndef __Scene_h
#define __Scene_h

#include "AbstractScene.h"
#include "Quaternion.h"


class Scene : public AbstractScene
{
public:
	Scene();
	virtual ~Scene();
	// set fovy angle
	void SetFovy(float fovy);

	// set params for roate amd zoom (debugging for rendering)
	void SetRotateRadPerPixel(float radPerPixel) { m_RadPerPixel=radPerPixel;}
	float GetRotateRadPerPixel() { return m_RadPerPixel; }

	void SetScaleRatePerPixel(float scaleRatePerPixel) { m_ScaleRatePerPixel=scaleRatePerPixel;}
	float GetScaleRatePerPixel() { return m_ScaleRatePerPixel; }

	// current parames and statues
	void GetTranslate(float translate[3]);//maybe useful!
	float GetScale() { return m_Scale;}
	const Quaternion& GetRotation() { return m_Rotation; }

	void SetTranslate(const float translate[3]);
	void SetScale(float scale);
	void SetRotation(const Quaternion& rotation);

	// current matrix
	const Matrix4x4& GetViewMatrix();
	const Matrix4x4& GetProjectionMatrix();

	// z distance, the x-ray source to isocenter of volume
	float GetZTranslateBase() { return m_ZTranslateBase;}

	// size of view port
	void GetViewPort(int* viewport) { viewport[0]=m_ViewPort[0]; viewport[1]=m_ViewPort[1]; }





	virtual void Resize(int width, int height);

	

protected:
	// render context
	virtual void RenderContent(int iframe);
	// box (ROI) for object
	virtual void GetSize3D(float size3D[3]);
	void SetSize3DModified();
	int m_OldMouseX,m_OldMouseY;


public:
	// initialization
	void _calculateCamera();
	//  update the matrix
	void _calculateMatrix();
	// reference point
	void _makeReference();


public:
	// size of view port
	int m_ViewPort[2];

	
	//float m_Fovy;

	// nearest plane and farest plane
	float m_Near, m_Far;

	// z offset
	float m_ZTranslateBase;

	// translate, scale, and rotation 
	float m_Translate[3];
	float m_Scale;
	Quaternion m_Rotation;

	// tag for whether update is needed?
	bool m_NeedCalculateCamera; // initialization again
	bool m_NeedCalculateMatrix; // recomputing matrix 

	// matrix
	Matrix4x4 m_ViewMatrix;
	Matrix4x4 m_ProjectionMatrix;

	////// params for interaction (in rendering)

	
	float m_RadPerPixel;
	float m_ScaleRatePerPixel;
	
	
	float m_OldTranslate[3];

	float m_OldScale;
	Quaternion m_OldRotation;  

	
	bool m_Rotating;
	bool m_Scaling;
	bool m_Moving;
	// distance
	float m_sid;
	float m_sdd;
	double m_halfFrustumHeight;
	double m_halfFrustumWidth;
	double m_ds;
};

#endif