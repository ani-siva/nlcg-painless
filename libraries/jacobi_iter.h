#include<stdio.h>
#include<stdbool.h>

#ifndef _JACOBI_ITER_H
	#define _JACOBI_ITER_H
		
		struct jacobi_iter{
			int row,col;
			const char* matfile;
			const char* arrfile;
			const char* initfile;
			double errlimit;
			double** A;
			double* b;
			double* x;

			void(*calculate)(struct jacobi_iter *this);
			bool(*check_err)(double* newx,double* x,int row,double errlimit);
		};

		extern const struct jacobi_iterClass{
			struct jacobi_iter(*new)(const char* matfile,const char* arrfile,const char* initfile,int row,int col,double errlimit);
		}jacobi_iter;

#endif
