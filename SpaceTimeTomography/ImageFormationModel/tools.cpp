
#include "tools.h"


int index3DtoLinear( int *sizeMat, int i, int j,int k)
{
	return (int)(i + j*sizeMat[0] + k*sizeMat[0]*sizeMat[1]);
}

int linearTo3Di( int *sizeMat, int index)
{
	return index%sizeMat[0];
}

int linearTo3Dj( int *sizeMat, int index)
{
	return (index-linearTo3Di(sizeMat,index))/sizeMat[0] % sizeMat[1];
}

int linearTo3Dk( int *sizeMat, int index)
{
	return ((index-linearTo3Di(sizeMat,index))/sizeMat[0] - linearTo3Dj(sizeMat,index)) / sizeMat[1] % sizeMat[2];
}

float myAbs(float x)
{
	return x > 0 ? x : -x;
}

float myMin(float a, float b)
{
	return a > b ? b : a;
}

float myMax(float a, float b)
{
	return a < b ? b : a;
}

float myPow2(float x)
{
	return x * x;
}

float dxc3( float *data,  float *u,  int *sizeMat, int i, int j, int k)
{
	int indexIP = index3DtoLinear(sizeMat, i + 1, j,k);
	int indexIM = index3DtoLinear(sizeMat, i - 1, j,k);
	if (i > 0 && i < sizeMat[0] - 1 && k < sizeMat[2] - 1)
	{
		return 0.5f*(data[indexIP] - data[indexIM]);
	}
	else
	{
		return 0.0f;
	}
}

float dxcT3( float *data,  float *u,  int *sizeMat, int i, int j, int k )
{
	//Attention! Transposed operator contains a different weight for u
	int indexIP = index3DtoLinear(sizeMat, i + 1, j, k);
	int indexIM = index3DtoLinear(sizeMat, i - 1, j, k);
	int indexIJ = index3DtoLinear(sizeMat, i, j, k);

	if (i > 1 && i < sizeMat[0] - 2 && k < sizeMat[2] - 1)
	{
		return 0.5f*(data[indexIP] * u[indexIP] - data[indexIM] * u[indexIM]);
	}
	else if (i <= 1 && k < sizeMat[2] - 1)
	{
		return 0.5f* data[indexIP] * u[indexIP];
	}
	else if (k < sizeMat[2] - 1)
	{
		return -0.5f* data[indexIM] * u[indexIM];
	}
	else
	{
		return 0.0f;
	}
}

float dyc3( float *data,  float *u,  int *sizeMat, int i, int j, int k)
{
	int indexJP = index3DtoLinear(sizeMat, i, j + 1, k);
	int indexJM = index3DtoLinear(sizeMat, i, j - 1, k);
	if (j > 0 && j < sizeMat[1] - 1 && k < sizeMat[2] - 1)
	{
		return 0.5f*(data[indexJP] - data[indexJM]);
	}
	else
	{
		return 0.0f;
	}
}

float dycT3( float *data,  float *u,  int *sizeMat, int i, int j, int k)
{
	//Attention! Transposed operator contains a different weight for u
	int indexIP = index3DtoLinear(sizeMat, i, j + 1, k);
	int indexIM = index3DtoLinear(sizeMat, i, j - 1, k);
	int indexIJ = index3DtoLinear(sizeMat, i, j, k);

	if (j > 1 && j < sizeMat[1] - 2 && k < sizeMat[2] - 1)
	{
		return 0.5f*(data[indexIP] * u[indexIP] - data[indexIM] * u[indexIM]);
	}
	else if (j <= 1 && k < sizeMat[2] - 1)
	{
		return 0.5f* data[indexIP] * u[indexIP];
	}
	else if (k < sizeMat[2] - 1)
	{
		return -0.5f* data[indexIM] * u[indexIM];
	}
	else
	{
		return 0.0f;
	}
}




float dxp3( float *data,  int *sizeMat, int i, int j, int k)
{
	return i < sizeMat[0] - 1 ? data[index3DtoLinear(sizeMat, i + 1, j, k)] - data[index3DtoLinear(sizeMat, i, j, k)] : 0.0f;
}

float dyp3( float *data,  int *sizeMat, int i, int j, int k)
{
	return j < sizeMat[1] - 1 ? data[index3DtoLinear(sizeMat, i, j + 1, k)] - data[index3DtoLinear(sizeMat, i, j, k)] : 0.0f;
}
float dzp3( float *data,  int *sizeMat, int i, int j, int k)
{
	return k < sizeMat[2] - 1 ? data[index3DtoLinear(sizeMat, i, j , k+1)] - data[index3DtoLinear(sizeMat, i, j, k)] : 0.0f;
}
//forward differenc  
float dtp3( float *data,  int *sizeMat, int i, int j, int k)
{
	return k < sizeMat[2] - 1 ? data[index3DtoLinear(sizeMat, i, j, k + 1)] - data[index3DtoLinear(sizeMat, i, j, k)] : 0.0f;
}

// backward difference
float dxm3( float *data,  int *sizeMat, int i, int j, int k)
{
	if (i > 0 && i < sizeMat[0] - 1)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)] - data[index3DtoLinear(sizeMat, i - 1, j, k)];
	}
	else if (i == 0)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)];
	}
	else
	{
		return -data[index3DtoLinear(sizeMat, i - 1, j, k)];
	}
}

float dym3( float *data,  int *sizeMat, int i, int j, int k)
{
	if (j > 0 && j < sizeMat[1] - 1)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)] - data[index3DtoLinear(sizeMat, i, j - 1, k)];
	}
	else if (j == 0)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)];
	}
	else
	{
		return -data[index3DtoLinear(sizeMat, i, j - 1, k)];
	}
}

float dzm3( float *data,  int *sizeMat, int i, int j, int k)
{
	if (k > 0 && k < sizeMat[2] - 1)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)] - data[index3DtoLinear(sizeMat, i, j, k-1)];
	}
	else if (k == 0)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)];
	}
	else
	{
		return -data[index3DtoLinear(sizeMat, i, j, k-1)];
	}
}


float dtm3( float *data,  int *sizeMat, int i, int j, int k)
{
	if (k > 0 && k < sizeMat[2] - 1)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)] - data[index3DtoLinear(sizeMat, i, j, k - 1)];
	}
	else if (k == 0)
	{
		return data[index3DtoLinear(sizeMat, i, j, k)];
	}
	else
	{
		return -data[index3DtoLinear(sizeMat, i, j, k - 1)];
	}
}
