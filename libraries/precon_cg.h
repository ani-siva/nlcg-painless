#include<stdio.h>
#ifndef _PRECON_CG_H
	#define _PRECON_CG_H
		struct precon_cg{
			int row,col;
			const char* matfile;
			const char* arrfile;
			const char* initfile;
			const char* preconfile;

			double** A;
			double* b;
			double** M;
			double** MInvA;
			double* x;
			double** u;
			double errlimit;
			void(*calculate)(struct precon_cg *this);

		};

		extern const struct precon_cgClass{
			struct precon_cg(*new)(const char* matfile,const char* arrfile,const char* initfile,const char* preconfile,int row,int col,double errlimit);
		}precon_cg;
#endif
