#include "GLUTSceneContainer.h"
#include "AbstractScene.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif

//static int detectorX ;
//static int detectorY ;
static int W_WIDTH;
static int W_HEIGHT;
static double W_ds;
#include <iostream>
using namespace  std;
static bool  IsNeedR=true;
GLUTSceneContainer* GLUTSceneContainer::GetMainWindow(int *wh, double ds)
{
	W_ds=ds;
	W_WIDTH=wh[0];
	W_HEIGHT=wh[1];
	//detectorX=w*W_ds;
	//detectorY=h*W_ds;
	if (!s_pMainWindow)
		s_pMainWindow=new GLUTSceneContainer;
	return s_pMainWindow;
}

void GLUTSceneContainer::SetSize(int width,int height)
{
	m_Size[0]=width;
	m_Size[1]=height;
//


}

void GLUTSceneContainer::GetSize(int &width,int &height) const
{
	width=m_Size[0];
	height=m_Size[1];

	
}

void GLUTSceneContainer::SetTitle(const char* title)
{
	strcpy(m_Title,title);
}

const char* GLUTSceneContainer::GetTitle() const
{
	return m_Title;
}

void GLUTSceneContainer::SetScene(AbstractScene* scene)
{
	m_Scene=scene;
}

AbstractScene* GLUTSceneContainer::GetScene() const
{
	return m_Scene;
}

void GLUTSceneContainer::Create(int posX,int posY)
{
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
	glutInitWindowPosition(posX,posY);
	glutInitWindowSize(W_WIDTH,W_HEIGHT);
	m_WindowID=glutCreateWindow(m_Title);
	if (m_FullScreen) glutFullScreen();
	glutDisplayFunc(s_display);
//	glutKeyboardFunc(s_keyboard);
	glutSpecialFunc(s_speckey);
	glutReshapeFunc(s_reshape);
//	glutMouseFunc(s_mouse);
	glutMotionFunc(s_motion);
	
	if (m_Scene) m_Scene->Init();
	glutMainLoop();
}

void GLUTSceneContainer::s_display()
{
	if(IsNeedR)
	{
	if (s_pMainWindow->m_Scene) s_pMainWindow->m_Scene->Render();
	glutSwapBuffers();
	IsNeedR=false;
	}
}



void GLUTSceneContainer::s_speckey(GLint key, GLint x, GLint y)
{
	bool redis=false;
	if (s_pMainWindow->m_Scene) s_pMainWindow->m_Scene->Speckey(key,x,y,redis);
	if (redis) glutPostRedisplay(); 
}

void GLUTSceneContainer::s_reshape(int width, int height)
{

	if (s_pMainWindow->m_Scene) s_pMainWindow->m_Scene->Resize(W_WIDTH,W_HEIGHT); 
	
	glutPostRedisplay();
}



void GLUTSceneContainer::s_motion(int x, int y)
{
	bool redis=false;
	if (s_pMainWindow->m_Scene) s_pMainWindow->m_Scene->MouseMove(x,y,redis);
	if (redis) glutPostRedisplay();
}

GLUTSceneContainer::GLUTSceneContainer()
{
	
	m_Size[0]=W_WIDTH;
	m_Size[1]=W_HEIGHT;
	m_Title[0]=0;
	m_FullScreen=false;
	m_WindowID=0;
	m_Scene=0;
}

GLUTSceneContainer::~GLUTSceneContainer()
{
}

void GLUTSceneContainer::Exit(int i)
{
	if (s_pMainWindow) delete s_pMainWindow;
	exit(i);
}

GLUTSceneContainer* GLUTSceneContainer::s_pMainWindow=0;