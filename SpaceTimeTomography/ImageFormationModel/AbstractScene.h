#ifndef Abstract_Scene_h
#define Abstract_Scene_h

class AbstractScene
{
public:
	AbstractScene(){}
	virtual ~AbstractScene(){}

	// Rendering, Resize, recycling
	virtual void Init(){};
	virtual void Render()=0;
	virtual void Resize(int width, int height)=0;
	virtual void Destroy(){};
	// CG standard functions, for update and display the windows for debugging
	virtual void Keyboard ( unsigned char key, int x, int y, bool &redisplay) {}
	virtual void Speckey(int key, int x, int y, bool &redisplay) {}

	virtual void LButtonDown(int x, int y, bool &redisplay){}
	virtual void LButtonUp(int x, int y, bool &redisplay){}

	virtual void MButtonDown(int x, int y, bool &redisplay){}
	virtual void MButtonUp(int x, int y, bool &redisplay){}

	virtual void MButtonDown(int x, int y){}
	virtual void MButtonUp(int x, int y){}


	virtual void RButtonDown(int x, int y, bool &redisplay){}
	virtual void RButtonUp(int x, int y, bool &redisplay){}

	virtual void Wheel(int value, bool &redisplay){}

	virtual void MouseMove(int x, int y, bool &redisplay){}

};

#endif
