#include<stdio.h>

#ifndef _GRAM_SCHMIDT_CG_H
	#define _GRAM_SCHMIDT_CG_H

		struct gram_schmidt_cg{
				
			int row,col;
			const char* matfile;
			const char* arrfile;
			const char* initfile;
			const char* ufile;
			double** A;
			double* b;
			double* x;
			double** u;

			void(*calculate)(struct gram_schmidt_cg *this);
			double** (*gram_schmidt_orthogonal)(double** mat,double** u,int row,int col);
		};

		extern const struct gram_schmidt_cgClass{
			struct gram_schmidt_cg(*new)(const char* matfile,const char* arrfile,const char* initfile,const char* ufile,int row,int col);
		}gram_schmidt_cg;

#endif
			
