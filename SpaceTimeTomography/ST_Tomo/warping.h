#ifndef WARPING_H
#define WARPING_H

#include<thread>

#include"vec3.h"

namespace warping {
	
	template< typename real, typename real3 >
	inline real3 point_grid_to_world( const int *dim, const real *aabb, const real3 &grid ){
		return real3(aabb[0] + grid[0]*(aabb[1]-aabb[0])/real(dim[0]),
					 aabb[2] + grid[1]*(aabb[3]-aabb[2])/real(dim[1]),
					 aabb[4] + grid[2]*(aabb[5]-aabb[4])/real(dim[2]) );
	}
	
	template< typename real, typename real3 >
	inline real3 world_to_point_grid( const int *dim, const real *aabb, const real3 &world ){
		return real3(real(dim[0])*(world[0]-aabb[0])/(aabb[1]-aabb[0]),
					 real(dim[1])*(world[1]-aabb[2])/(aabb[3]-aabb[2]),
					 real(dim[2])*(world[2]-aabb[4])/(aabb[5]-aabb[4]) );
	}
	
	template< typename real, typename real3 >
	inline real3 cell_grid_to_world( const int *dim, const real *aabb, const real3 &grid ){
		return real3(aabb[0] + (grid[0]+0.5)*(aabb[1]-aabb[0])/real(dim[0]),
					 aabb[2] + (grid[1]+0.5)*(aabb[3]-aabb[2])/real(dim[1]),
					 aabb[4] + (grid[2]+0.5)*(aabb[5]-aabb[4])/real(dim[2]) );
	}
	
	template< typename real, typename real3 >
	inline real3 world_to_cell_grid( const int *dim, const real *aabb, const real3 &world ){
		return real3(real(dim[0])*(world[0]-aabb[0])/(aabb[1]-aabb[0])-0.5,
					 real(dim[1])*(world[1]-aabb[2])/(aabb[3]-aabb[2])-0.5,
					 real(dim[2])*(world[2]-aabb[4])/(aabb[5]-aabb[4])-0.5 );
	}
	
	inline void check( const bool &val, const char *msg ){
		if( !val ){
			std::cout << "ERROR: " << msg << std::endl;
		}
	}
	
	template< typename real, typename real3 >
	inline bool point_in_bounding_box( const real* aabb, const real3 &p ){
		return p[0] >= aabb[0] && p[0] <= aabb[1] && p[1] >= aabb[2] && p[1] <= aabb[3] && p[2] >= aabb[4] && p[2] <= aabb[5];
	}
	
	template< typename real, typename real3, typename image >
	inline real sample_image( const int *dim, const real *aabb, const image &img, const real3 &p, const int c ){
		check( dim[0] == img.width() && dim[1] == img.height() && dim[2] == img.depth(), "grid dimensions sanity check" );
		check( aabb[1] > aabb[0] && aabb[3] > aabb[2] && aabb[5] > aabb[4], "bounding box sanity check" );
		
		real3 g = world_to_cell_grid( dim, aabb, p );
		
		check( dim[0] == img.width() && dim[1] == img.height() && dim[2] == img.depth(), "grid dimensions sanity check" );
		check( aabb[1] > aabb[0] && aabb[3] > aabb[2] && aabb[5] > aabb[4], "bounding box sanity check" );

		
		g[0] = std::max( real(1), std::min( real(dim[0]-2), g[0] ) );
		g[1] = std::max( real(1), std::min( real(dim[1]-2), g[1] ) );
		g[2] = std::max( real(1), std::min( real(dim[2]-2), g[2] ) );
		int ix=g[0], iy=g[1], iz=g[2];
		
		check( ix >= 0 && ix < img.width() && iy >= 0 && iy < img.height() && iz >= 0 && iz < img.depth(), "point clamped correctly" );
		//check(ix >= 0 && ix < img.width() && iy >= 0 && iy < img.height() && iz >= 0 && iz <= img.depth(), "point clamped correctly");

		real minim = img(ix,iy,iz,c);
		real maxim = minim;
		minim = std::min( minim, img(ix+1,iy+0,iz+0,c) );
		maxim = std::max( maxim, img(ix+1,iy+0,iz+0,c) );
		minim = std::min( minim, img(ix+1,iy+1,iz+0,c) );
		maxim = std::max( maxim, img(ix+1,iy+1,iz+0,c) );
		minim = std::min( minim, img(ix+0,iy+1,iz+0,c) );
		maxim = std::max( maxim, img(ix+0,iy+1,iz+0,c) );
		minim = std::min( minim, img(ix+0,iy+0,iz+1,c) );
		maxim = std::max( maxim, img(ix+0,iy+0,iz+1,c) );
		minim = std::min( minim, img(ix+1,iy+0,iz+1,c) );
		maxim = std::max( maxim, img(ix+1,iy+0,iz+1,c) );
		minim = std::min( minim, img(ix+1,iy+1,iz+1,c) );
		maxim = std::max( maxim, img(ix+1,iy+1,iz+1,c) );
		minim = std::min( minim, img(ix+0,iy+1,iz+1,c) );
		maxim = std::max( maxim, img(ix+0,iy+1,iz+1,c) );
		
		//return std::max( minim, std::min( maxim, img.cubic_atXYZ( g[0], g[1], g[2], c ) ) );
		return std::max(minim, std::min(maxim, img.cubic_atXYZ(g[0], g[1], g[2], c)));

	}
	
	template< typename real, typename real3, typename image >
	inline real3 sample_velocity( const int *dim, const real *aabb, const image &vel, const real3 &p ){
		real3 v;
		v[0] = sample_image( dim, aabb, vel, p, 0 );
		v[1] = sample_image( dim, aabb, vel, p, 1 );
		v[2] = sample_image( dim, aabb, vel, p, 2 );
		return v;
	}
	
	
		template< typename real, typename real3, typename image >
		inline real3 trace(const int *dim, const real *aabb, const image &vel, const real3 &p, const real dt, const int steps = 10) {
	/*	real3 pf = p;
		for( int i=0; i<steps; i++ ){
			pf += dt*sample_velocity(dim,aabb,vel,pf)/real(steps);
		}
		real3 pi = pf;
		for( int i=0; i<steps; i++ ){
			pi -= dt*sample_velocity(dim,aabb,vel,pi)/real(steps);
		}
		return pf + 0.5*(p-pi);*/
		real3 p2 = p + dt*sample_velocity( dim, aabb, vel, p );
		real3 p0 = p2 - dt*sample_velocity( dim, aabb, vel, p2 );
		return p2 + 0.5*(p-p0);
	}
	
	template< typename real, typename image >
	image warp( const int *dim, const real *aabb, const image &vel, const image &src, real dt ){
		int nc = src.spectrum();
		image out( dim[0], dim[1], dim[2], nc, 0.0 );
		
		
		auto warp_job = [nc,dim,aabb,&vel,&src,dt,&out]( int job_id, int n_jobs ) -> void {
			for( int k=job_id; k<dim[2]; k+=n_jobs ){
				for( int j=0; j<dim[1]; j++ ){
					for( int i=0; i<dim[0]; i++ ){
						geom::vec3<real> p = trace( dim, aabb, vel, cell_grid_to_world( dim, aabb, geom::vec3<real>(i,j,k) ), -dt );
						for( int c=0; c<nc; c++ ){
							out(i,j,k,c) = sample_image( dim, aabb, src, p, c );
						}
					}
				}
			}
		};
		const int njobs=4;
		std::thread job[njobs-1];
		for( int i=0; i<njobs-1; i++ ){
			job[i] = std::thread( std::bind( warp_job, i, njobs ) );
		}
		warp_job( njobs-1, njobs );
		for( int i=0; i<njobs-1; i++ ){
			job[i].join();
		}
		
		return out;
	}
	
};


#endif