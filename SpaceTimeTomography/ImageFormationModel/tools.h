


int index3DtoLinear( int *sizeMat, int i, int j, int k);
float myAbs(float x);
float myMin(float a, float b);
float myMax(float a, float b);
float myPow2(float x);

float dxp3( float *data,  int *sizeMat, int i, int j, int k);
float dyp3( float *data,  int *sizeMat, int i, int j, int k);
float dzp3( float *data,  int *sizeMat, int i, int j, int k);
float dtp3( float *data,  int *sizeMat, int i, int j, int k);
float dxm3( float *data,  int *sizeMat, int i, int j, int k);
float dym3( float *data,  int *sizeMat, int i, int j, int k);
float dzm3( float *data,  int *sizeMat, int i, int j, int k);
float dtm3( float *data,  int *sizeMat, int i, int j, int k);

float dxc3( float *data,  float *u,  int *sizeMat, int i, int j, int k);
float dyc3( float *data,  float *u,  int *sizeMat, int i, int j, int k);
float dxcT3( float *data,  float *u,  int *sizeMat, int i, int j, int k);
float dycT3( float *data,  float *u,  int *sizeMat, int i, int j, int k);


int linearTo3Di( int *sizeMat, int index);
int linearTo3Dj( int *sizeMat, int index);
int linearTo3Dk( int *sizeMat, int index);