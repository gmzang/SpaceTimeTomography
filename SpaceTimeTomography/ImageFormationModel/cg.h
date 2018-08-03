#ifndef CG_H
#define CG_H

#include<cmath>
#include<iostream>

#include<image_blas.h>

namespace image_cg {
	
	template< typename real, typename image, typename linear_op >
	void cg( linear_op &A, image &x, const image &b, const real rel_tol=1e-6, const int max_iters=500, const real abs_tol=1e-9 ){
		
		image Ap, p, r;
		real rr, rr_old, rr_init, alpha, norm;
		
		
		// evaluate the residual
		A(x,Ap);
		image_blas::add( image_blas::scaled(Ap,-1.0), b, r );
		image_blas::copy( r, p );
		rr = rr_old = rr_init = image_blas::dot<real>( r, r );
		if( (norm=sqrt(rr)) < abs_tol ){
			std::cout << "cg finished with residual: " << norm << std::endl;
			return;
		}

		for( int iter=0; iter<max_iters; iter++ ){
			// P = Ap
			A(p,Ap);
			
			// alpha = dot(r_old,r_old)/dot(p,Ap)
			alpha = rr_old/image_blas::dot<real>(p,Ap);
			
			// x += alpha*p
			image_blas::add( image_blas::scaled(p, alpha), x );
			
			// r -= alpha*Ap
			image_blas::add( image_blas::scaled(Ap,-alpha), r );
			
			// compute dot(r,r)
			rr = image_blas::dot<real>( r, r );
			norm = sqrt(rr);
			if( sqrt(rr/rr_init) < rel_tol ){
				std::cout << "cg_finished with residual: " << sqrt(rr/rr_init) << std::endl;
				return;
			} else {
				if( iter == 0 || iter % 10 == 0 )
					std::cout << "  iteration: " << iter << ": rel. resid.: " << sqrt(rr/rr_init) << std::endl;
			}
			// p = r + p*dot(r,r)/dot(r_old,r_old)
			image_blas::add( image_blas::scaled(p,rr/rr_old), r, p );
			rr_old = rr;
		}
	}

};

#endif