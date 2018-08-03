#ifndef OPTICAL_FLOW_H
#define OPTICAL_FLOW_H

#include<thread>

#define cimg_use_tiff
#define cimg_use_tif
#include "CImg.h"
#include"linear_solver.h"
namespace cl = cimg_library;
namespace ls = linear_solver;
#include<vec3.h>
#include<options.h>

#include"cg.h"
#include"warping.h"
#include "tools.h"

#include<scope.h>
utilities_scope_defines
int interpo = 5;
#include<string_utils.h>
#define DIMENSION_SIZE 6
int ROI[] = { 1,1,1,1,1,1 };

namespace optical_flow {

	class cell_indexer {
	private:
		int m_dim[3];
	public:
		cell_indexer(const int *dim) {
			m_dim[0] = dim[0];
			m_dim[1] = dim[1];
			m_dim[2] = dim[2];
		}

		inline bool is_boundary(int i, int j, int k) const {
			return i < 0 || i >= m_dim[0] || j < 0 || j >= m_dim[1] || k < 0 || k >= m_dim[2];
		}

		inline int num_cells() const {
			return m_dim[0] * m_dim[1] * m_dim[2];
		}

		inline int operator()(int i, int j, int k) const {
			i = std::max(0, std::min(m_dim[0] - 1, i));
			j = std::max(0, std::min(m_dim[1] - 1, j));
			k = std::max(0, std::min(m_dim[2] - 1, k));
			return i + m_dim[0] * (j + k*m_dim[1]);
		}
	};

	class cell_var_indexer : public cell_indexer {
	private:
		int m_nvars;
	public:
		cell_var_indexer(const int *dim, const int nvars) : cell_indexer(dim), m_nvars(nvars) { }
		inline int operator()(int i, int j, int k, int c) const {
			return cell_indexer::operator()(i, j, k)*m_nvars + c;
		}
		inline int num_vars() {
			return m_nvars*num_cells();
		}
	};


	template< typename real, typename image >
	image blur_and_downsample(const image &I, const int nx, const int ny, const int nz, const real sigma) {
		return I.get_blur(sigma).get_resize(nx, ny, nz, I.spectrum(), 5);
	}

	
	template< typename real, typename image >
	image blur_and_downsample(const image &I, const real eta, const real sigma) {
		return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth(), sigma);
	}
	template< typename real, typename image >
	image Y_blur_and_downsample(const image &I, const real eta, const real sigma) {
		return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth(), sigma);
	}

	template< typename real, typename image >
	void L1Huber_OF(int *dim, const image & image1, const image & image2, image & _inputV, image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4,
		float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps, real *_aabb, float huber)
	{

		const float tol = _tol;
		const float lambda = _lambda;

		const int maxIterations = _maxIterations;


		int typeNorm = _norm;


		float stepsize[4] = { 1.0f, 1.0f, 1.0f ,1.0f };
		float stepsizeD[4] = { 1.0f / stepsize[0],
			1.0f / stepsize[1],
			1.0f / stepsize[2],
			1.0f / stepsize[3] };

		int numberOfWarps = _numberOfWarps;


		float huberEpsilon = huber;



		int gradientConstancy = 0;
	
		const int nPx = (int)(dim[0] * dim[1] * dim[2]);
		cout << "nPX" << nPx << endl;
		

		float* v1 = new float[nPx];
		float* v2 = new float[nPx];
		float* v3 = new float[nPx];

		float* v1Old = new float[nPx];
		float* v2Old = new float[nPx];
		float* v3Old = new float[nPx];

		float* image1f = new float[nPx];
		float* image2f = new float[nPx];

		float* ux = new float[nPx];
		float* uy = new float[nPx];
		float* uz = new float[nPx];
		float* ut = new float[nPx];
		

		float* y11 = new float[nPx];
		float* y12 = new float[nPx];
		float* y13 = new float[nPx];

		float* y21 = new float[nPx];
		float* y22 = new float[nPx];
		float* y23 = new float[nPx];

		float* y31 = new float[nPx];
		float* y32 = new float[nPx];
		float* y33 = new float[nPx];
		float* y4 = new float[nPx];
		

		float* Kty1 = new float[nPx];
		float* Kty2 = new float[nPx];
		float* Kty3 = new float[nPx];

		float* Kty1Old = new float[nPx];
		float* Kty2Old = new float[nPx];
		float* Kty3Old = new float[nPx];

		float* Kx11 = new float[nPx];
		float* Kx12 = new float[nPx];
		float* Kx13 = new float[nPx];

		float* Kx21 = new float[nPx];
		float* Kx22 = new float[nPx];
		float* Kx23 = new float[nPx];

		float* Kx31 = new float[nPx];
		float* Kx32 = new float[nPx];
		float* Kx33 = new float[nPx];
		float* Kx4 = new float[nPx];
		

		float sigma1 = myMin(stepsize[3] / 3.0f, myMin(stepsize[1] / 3.0f, stepsize[2] / 3.0f));
		float* sigma2 = new float[nPx];

	

		float* tau1 = new float[nPx];
		float* tau2 = new float[nPx];
		float* tau3 = new float[nPx];

		int * tableI = new int[nPx];
		int * tableJ = new int[nPx];
		int * tableK = new int[nPx];

		//Huber Factor

		const float huberFactor = 1.0f / (1.0f + sigma1* huberEpsilon / lambda);

		//residuals
		float p = 0.0f;
		float d = 0.0f;

		#pragma omp parallel for
		for (int k = 0; k < dim[2]; ++k)
		{
			for (int j = 0; j < dim[1]; ++j)
			{
				for (int i = 0; i < dim[0]; ++i)
				{
					//int tmpIndex = index2DtoLinear(dim, i, j);
					int tmpIndex = index3DtoLinear(dim, i, j, k);

					tableI[tmpIndex] = i;
					tableJ[tmpIndex] = j;
					tableK[tmpIndex] = k;

				
					v1[tmpIndex] = _inputV(i, j, k, 0);
					v2[tmpIndex] = _inputV(i, j, k, 1);
					v3[tmpIndex] = _inputV(i, j, k, 2);
		
					Kty1[tmpIndex] = 0.0f;
					Kty2[tmpIndex] = 0.0f;
					Kty3[tmpIndex] = 0.0f;

					
					y11[tmpIndex] = _inputY1(i, j, k, 0);
					y12[tmpIndex] = _inputY1(i, j, k, 1);
					y13[tmpIndex] = _inputY1(i, j, k, 2);

					y21[tmpIndex] = _inputY2(i, j, k, 0);
					y22[tmpIndex] = _inputY2(i, j, k, 1);
					y23[tmpIndex] = _inputY2(i, j, k, 2);

					y31[tmpIndex] = _inputY3(i, j, k, 0);
					y32[tmpIndex] = _inputY3(i, j, k, 1);
					y33[tmpIndex] = _inputY3(i, j, k, 2);

					y4[tmpIndex] = _inputY4(i, j, k);

				

					Kx11[tmpIndex] = 0.0f;
					Kx12[tmpIndex] = 0.0f;
					Kx13[tmpIndex] = 0.0f;

					Kx21[tmpIndex] = 0.0f;
					Kx22[tmpIndex] = 0.0f;
					Kx23[tmpIndex] = 0.0f;

					Kx31[tmpIndex] = 0.0f;
					Kx32[tmpIndex] = 0.0f;
					Kx33[tmpIndex] = 0.0f;

					Kx4[tmpIndex] = 0.0f;

					
				}
			}
		}
		cout << "Finish the initialization, open clock " << endl;
		clock_t begin = clock();
		//do k warpings
		for (int idxwarp = 0; idxwarp < numberOfWarps; ++idxwarp)
		{


		
			real dt = 1;
			image _warpI2 = warping::warp(dim, _aabb, _inputV, image2, -dt);
			imagelist warpI2_xyz = _warpI2.get_gradient("xyz", 0);
			

			#pragma omp parallel for
			for (int k = ROI[4]; k < dim[2] - ROI[5]; ++k)
			{
				for (int j = ROI[2]; j < dim[1] - ROI[3]; ++j)
				{
					for (int i = ROI[0]; i < dim[0] - ROI[1]; ++i)
					{
						int tmpIndex = index3DtoLinear(dim, i, j, k);
						
						ux[tmpIndex] = warpI2_xyz[0](i, j, k);
						uy[tmpIndex] = warpI2_xyz[1](i, j, k);
						uz[tmpIndex] = warpI2_xyz[2](i, j, k);
						
						ut[tmpIndex] = _warpI2(i, j, k) - image1(i, j, k) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex] - uz[tmpIndex] * v3[tmpIndex];

					}
				}
			}


			#pragma omp parallel for
			for (int i = 0; i < nPx; ++i)
			{
			
				tau1[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(ux[i]);
				tau2[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(uy[i]);
				tau3[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(uz[i]);
				
				sigma2[i] = std::abs(ux[i]) + std::abs(uy[i]) + std::abs(uz[i]);
				
				tau1[i] = 1.0f / tau1[i];
				tau2[i] = 1.0f / tau2[i];
				tau3[i] = 1.0f / tau3[i];
				sigma2[i] = 1.0f / sigma2[i];

			
			}


			int iterations = 0;
			float err = 1.0f;

			while (err > tol && iterations <= maxIterations)
			{
				++iterations;

				if (iterations % 50 == 0)
				{
					p = 0.0f;
					d = 0.0f;
				}

			
				#pragma omp parallel for
				for (int k = ROI[4]; k < dim[2] - ROI[5]; ++k)
				{
					for (int j = ROI[2]; j < dim[1] - ROI[3]; ++j)
					{
						for (int i = ROI[0]; i < dim[0] - ROI[1]; ++i)
						{
							int tmpIndex = index3DtoLinear(dim, i, j, k);

							
							if (iterations % 50 == 0)
							{
								v1Old[tmpIndex] = v1[tmpIndex];
								v2Old[tmpIndex] = v2[tmpIndex];
								v3Old[tmpIndex] = v3[tmpIndex];

								Kty1Old[tmpIndex] = Kty1[tmpIndex];
								Kty2Old[tmpIndex] = Kty2[tmpIndex];
								Kty3Old[tmpIndex] = Kty3[tmpIndex];
							}

							//transpose equals -div  
							Kty1[tmpIndex] = -stepsizeD[1] * dxm3(y11, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[2] * dym3(y12, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[3] * dzm3(y13, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
								ux[tmpIndex] * y4[tmpIndex];

							Kty2[tmpIndex] = -stepsizeD[1] * dxm3(y21, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[2] * dym3(y22, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[3] * dzm3(y23, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
								uy[tmpIndex] * y4[tmpIndex];

							Kty3[tmpIndex] = -stepsizeD[1] * dxm3(y31, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[2] * dym3(y32, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[3] * dzm3(y33, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
								uz[tmpIndex] * y4[tmpIndex];


							v1[tmpIndex] -= tau1[tmpIndex] * Kty1[tmpIndex];
							v2[tmpIndex] -= tau2[tmpIndex] * Kty2[tmpIndex];
							v3[tmpIndex] -= tau3[tmpIndex] * Kty3[tmpIndex];
							_inputV(i, j, k, 0) = v1[tmpIndex]; 
							_inputV(i, j, k, 1) = v2[tmpIndex];
							_inputV(i, j, k, 2) = v3[tmpIndex];
															 


							if (iterations % 50 == 0)
							{
								//residuals
								p += std::abs((v1Old[tmpIndex] - v1[tmpIndex]) / tau1[tmpIndex] - Kty1Old[tmpIndex] + Kty1[tmpIndex])
									+ std::abs((v2Old[tmpIndex] - v2[tmpIndex]) / tau2[tmpIndex] - Kty2Old[tmpIndex] + Kty2[tmpIndex])
									+ std::abs((v3Old[tmpIndex] - v3[tmpIndex]) / tau3[tmpIndex] - Kty3Old[tmpIndex] + Kty3[tmpIndex]);

							}
						}
					}
				}
				

				#pragma omp parallel for reduction(+:d)
				for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
				{
					float y11Tilde, y12Tilde, y13Tilde, y21Tilde, y22Tilde, y23Tilde, y31Tilde, y32Tilde, y33Tilde;
					float y11Old, y12Old, y13Old, y21Old, y22Old, y23Old, y31Old, y32Old, y33Old, y4Old;
					float Kx11Old, Kx12Old, Kx13Old, Kx21Old, Kx22Old, Kx23Old, Kx31Old, Kx32Old, Kx33Old, Kx4Old;
					if (iterations % 50 == 0)
					{


						y11Old = y11[tmpIndex];
						y12Old = y12[tmpIndex];
						y13Old = y13[tmpIndex];
						y21Old = y21[tmpIndex];
						y22Old = y22[tmpIndex];
						y23Old = y23[tmpIndex];
						y31Old = y31[tmpIndex];
						y32Old = y32[tmpIndex];
						y33Old = y33[tmpIndex];

						y4Old = y4[tmpIndex];
					}

					Kx11Old = Kx11[tmpIndex];
					Kx12Old = Kx12[tmpIndex];
					Kx13Old = Kx13[tmpIndex];
					Kx21Old = Kx21[tmpIndex];
					Kx22Old = Kx22[tmpIndex];
					Kx23Old = Kx23[tmpIndex];
					Kx31Old = Kx31[tmpIndex];
					Kx32Old = Kx32[tmpIndex];
					Kx33Old = Kx33[tmpIndex];

					Kx4Old = Kx4[tmpIndex];

					
					Kx11[tmpIndex] = stepsizeD[1] * dxp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx12[tmpIndex] = stepsizeD[2] * dyp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx13[tmpIndex] = stepsizeD[3] * dzp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);

					Kx21[tmpIndex] = stepsizeD[1] * dxp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx22[tmpIndex] = stepsizeD[2] * dyp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx23[tmpIndex] = stepsizeD[3] * dzp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);

					Kx31[tmpIndex] = stepsizeD[1] * dxp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx32[tmpIndex] = stepsizeD[2] * dyp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx33[tmpIndex] = stepsizeD[3] * dzp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);

					Kx4[tmpIndex] = ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex] + uz[tmpIndex] * v3[tmpIndex];  
																															

					if (typeNorm == 4) // Huber
					{
						y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old)) * huberFactor;
						y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old)) * huberFactor;
						y13Tilde = (y13[tmpIndex] + sigma1*(Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old)) * huberFactor;

						y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old)) * huberFactor;
						y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old)) * huberFactor;
						y23Tilde = (y23[tmpIndex] + sigma1*(Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old)) * huberFactor;

						y31Tilde = (y31[tmpIndex] + sigma1*(Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old)) * huberFactor;
						y32Tilde = (y32[tmpIndex] + sigma1*(Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old)) * huberFactor;
						y33Tilde = (y33[tmpIndex] + sigma1*(Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old)) * huberFactor;
					}
					else
					{
						

						y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old));
						y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old));
						y13Tilde = (y13[tmpIndex] + sigma1*(Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old));

						y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old));
						y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old));
						y23Tilde = (y23[tmpIndex] + sigma1*(Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old));

						y31Tilde = (y31[tmpIndex] + sigma1*(Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old));
						y32Tilde = (y32[tmpIndex] + sigma1*(Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old));
						y33Tilde = (y33[tmpIndex] + sigma1*(Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old));
					}

					float divisor1 = std::max(1.0f, std::sqrt(y11Tilde*y11Tilde + y12Tilde*y12Tilde + y13Tilde*y13Tilde) / lambda);
					float divisor2 = std::max(1.0f, std::sqrt(y21Tilde*y21Tilde + y22Tilde*y22Tilde + y23Tilde*y23Tilde) / lambda);
					float divisor3 = std::max(1.0f, std::sqrt(y31Tilde*y31Tilde + y32Tilde*y32Tilde + y33Tilde*y33Tilde) / lambda);


					y11[tmpIndex] = y11Tilde / divisor1;
					y12[tmpIndex] = y12Tilde / divisor1;
					y13[tmpIndex] = y13Tilde / divisor1;

					y21[tmpIndex] = y21Tilde / divisor2;
					y22[tmpIndex] = y22Tilde / divisor2;
					y23[tmpIndex] = y23Tilde / divisor2;

					y31[tmpIndex] = y31Tilde / divisor3;
					y32[tmpIndex] = y32Tilde / divisor3;
					y33[tmpIndex] = y33Tilde / divisor3;



					y4[tmpIndex] = std::max(-1.0f, std::min(1.0f, y4[tmpIndex] + sigma2[tmpIndex] * (Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old + ut[tmpIndex])));

					if (iterations % 50 == 0)
					{

					
						d += std::abs((y11Old - y11[tmpIndex]) / sigma1 - Kx11Old + Kx11[tmpIndex]) +
							std::abs((y12Old - y12[tmpIndex]) / sigma1 - Kx12Old + Kx12[tmpIndex]) +
							std::abs((y13Old - y13[tmpIndex]) / sigma1 - Kx13Old + Kx13[tmpIndex]) +

							std::abs((y21Old - y21[tmpIndex]) / sigma1 - Kx21Old + Kx21[tmpIndex]) +
							std::abs((y22Old - y22[tmpIndex]) / sigma1 - Kx22Old + Kx22[tmpIndex]) +
							std::abs((y23Old - y23[tmpIndex]) / sigma1 - Kx23Old + Kx23[tmpIndex]) +

							std::abs((y31Old - y31[tmpIndex]) / sigma1 - Kx31Old + Kx31[tmpIndex]) +
							std::abs((y32Old - y32[tmpIndex]) / sigma1 - Kx32Old + Kx32[tmpIndex]) +
							std::abs((y33Old - y33[tmpIndex]) / sigma1 - Kx33Old + Kx33[tmpIndex]) +

							std::abs((y4Old - y4[tmpIndex]) / sigma2[tmpIndex] - Kx4Old + Kx4[tmpIndex]);

				
					}
				}

				if (iterations % 50 == 0)
				{
					err = (d*d + p*p) / (float)nPx;
				}

				if (iterations % 1000 == 0)
				{

					cout << "Iteration: " << iterations << " Residual " << err << endl;

					//mexPrintf("Iteration %d,Residual %e\n", iterations, err);
					//mexEvalString("drawnow;");
				}
			}
		}


#pragma omp parallel for
		for (int k = ROI[4]; k < dim[2] - ROI[5]; ++k)
		{
			for (int j = ROI[2]; j < dim[1] - ROI[3]; ++j)
			{
				for (int i = ROI[0]; i < dim[0] - ROI[1]; ++i)
				{
					int tmpIndex = index3DtoLinear(dim, i, j, k);
				

					_inputY1(i, j, k, 0) = (float)y11[tmpIndex];
					_inputY1(i, j, k, 1) = (float)y12[tmpIndex];
					_inputY1(i, j, k, 2) = (float)y13[tmpIndex];

					_inputY2(i, j, k, 0) = (float)y21[tmpIndex];
					_inputY2(i, j, k, 1) = (float)y22[tmpIndex];
					_inputY2(i, j, k, 2) = (float)y23[tmpIndex];

					_inputY3(i, j, k, 0) = (float)y31[tmpIndex];
					_inputY3(i, j, k, 1) = (float)y32[tmpIndex];
					_inputY3(i, j, k, 2) = (float)y33[tmpIndex];

					_inputY4(i, j, k, 0) = (float)y4[tmpIndex];

					
					_inputV(i, j, k, 0) = (float)v1[tmpIndex];
					_inputV(i, j, k, 1) = (float)v2[tmpIndex];
					_inputV(i, j, k, 2) = (float)v3[tmpIndex];

				}

			}
		}
		cout << " in the L1TVOpticalFlowNonlinear  5 " << endl;

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cout << "=========Total time is :" << elapsed_secs << " s ======" << endl;
		int num1 = _inputV.width();
		int num2 = _inputV.height();
		int num3 = _inputV.depth();


		char out_dir[1024];
		char basename[1024];
		char flow_field[1024];
		//sprintf(out_dir, ".");
		sprintf(basename, "outputf1f2Roi");


		delete[] tableI;
		delete[] tableJ;
		delete[] tableK;

		delete[] image1f;
		delete[] image2f;

		delete[] sigma2;

		delete[] tau1;
		delete[] tau2;
		delete[] tau3;

		delete[] v1;
		delete[] v2;
		delete[] v3;
		delete[] v1Old;
		delete[] v2Old;
		delete[] v3Old;

		delete[] ux;
		delete[] uy;
		delete[] uz;
		delete[] ut;


		delete[] y11;
		delete[] y12;
		delete[] y13;
		delete[] y21;
		delete[] y22;
		delete[] y23;
		delete[] y31;
		delete[] y32;
		delete[] y33;
		delete[] y4;


		delete[] Kty1;
		delete[] Kty2;
		delete[] Kty3;

		delete[] Kty1Old;
		delete[] Kty2Old;
		delete[] Kty3Old;

		delete[] Kx11;
		delete[] Kx12;
		delete[] Kx13;
		delete[] Kx21;
		delete[] Kx22;
		delete[] Kx23;
		delete[] Kx31;
		delete[] Kx32;
		delete[] Kx33;
		delete[] Kx4;

	}






	template< typename real, typename image >
	image L1Huber_OF_MultiScale3D(int *dim, const image &image1, const image &image2, image &_inputV,
		image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4, real _tol, real _lambda, int _maxIterations, int _norm, int _numberOfWarps,
		real eta, real sigma, int scales, real _huber, int *cut)
	{


		image I1_blur = image1.get_blur(sigma);
		image I2_blur = image2.get_blur(sigma);


		int size = image1.width()*image1.height()*image1.depth();

		imagelist Im1s;
		imagelist Im2s;
		imagelist inputVs;

		imagelist inputY1s;
		imagelist inputY2s;
		imagelist inputY3s;
		imagelist inputY4s;

	


		//imagelist dims;
		int twidth = I1_blur.width();
		int theight = I1_blur.height();
		int tdepth = I1_blur.depth();
		float t_eta = 1.0f;

		int tdims[] = { twidth ,theight, tdepth };

		// reture at last second scale;
		for (int iscale = 0; iscale < scales; iscale++)
		{



			Im1s.insert(I1_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
			Im2s.insert(I2_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
			inputVs.insert(_inputV.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY1s.insert(_inputY1.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY2s.insert(_inputY2.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY3s.insert(_inputY3.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));

			t_eta = t_eta*eta;


		}
		Im1s.reverse();
		Im2s.reverse();
		inputVs.reverse();

		inputY1s.reverse();
		inputY2s.reverse();
		inputY3s.reverse();
		inputY4s.reverse();
	
		
		
	
	for (int iscale = 0; iscale < inputVs.size(); iscale++)
		{
			int levels = inputVs.size() - iscale;
			cout << "At levels " << levels << endl;
			int _dim[] = { Im1s[iscale].width(),Im1s[iscale].height(),Im1s[iscale].depth() };
			cout << "_dim[]  :" << _dim[0] << " " << _dim[1] << " " << _dim[2] << endl;
	
			real _aabb[] = { 0,tdims[0], 0,tdims[1], 0,tdims[2] };

			cout << "tdims[]  :" << tdims[0] << " " << tdims[1] << " " << tdims[2] << endl;

			
			ROI[0] =(float)cut[0]*((float)_dim[0] / (float)tdims[0]);
			ROI[1] = (float)cut[1]*((float)_dim[0] / (float)tdims[0]);
			ROI[2] = (float)cut[2] * ((float)_dim[1] / (float)tdims[1]);
			ROI[3] = (float)cut[3] * ((float)_dim[1] / (float)tdims[1]);
			ROI[4] = (float)cut[4] * ((float)_dim[2] / (float)tdims[2]);
			ROI[5] = (float)cut[5] * ((float)_dim[2] / (float)tdims[2]);
			cout << "ROI: " << ROI[0] << " " << ROI[1] << " " << ROI[2] << " " << ROI[3] << " " << ROI[4] << " " << ROI[5] << " " << endl;

			L1Huber_OF(_dim, Im1s[iscale], Im2s[iscale], inputVs[iscale],
				inputY1s[iscale], inputY2s[iscale], inputY3s[iscale], inputY4s[iscale], _tol, _lambda, _maxIterations, _norm, _numberOfWarps, _aabb, _huber);


			if (iscale < inputVs.size() - 1)
			{

				int _dims[] = { inputVs[iscale + 1].width(), inputVs[iscale + 1].height(), inputVs[iscale + 1].depth() };
				//			int _olddims[] = { inputVs[iscale].width(), inputVs[iscale].height(), inputVs[iscale].depth() };
		
				inputVs[iscale + 1] = inputVs[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);

				inputY1s[iscale + 1] = inputY1s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
				inputY2s[iscale + 1] = inputY2s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
				inputY3s[iscale + 1] = inputY3s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
				inputY4s[iscale + 1] = inputY4s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 1, interpo);
			}

		}

		return inputVs[inputVs.size()-1];
	}

}  // flow namespace


#endif



