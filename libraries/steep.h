#include<stdio.h>

#ifndef _STEEP_H
	#define _STEEP_H

		struct steep {
			int row,col;
			const char* matfile;
			const char* arrfile;
			const char* initfile;
			double errlimit;
			double** A;
			double* b;
			double* x;
			
			void(*calculate)(struct steep *this);
			int(*check_error)(double* xi,double* x,double errlimit,int row);
			double*(*residual)(double** A,double *b,double* x,int row,int col);
			double(*alpha)(double** A,double* res,int row, int col);
			double*(*xupdate)(double* x0,double alpha,double* res,int row);
				

		};
		
		extern const struct steepClass {
			struct steep (*new)(const char* matfile, const char* arrfile, const char* initfile, int row, int col, double errlimit);
		}steep;

#endif
			
