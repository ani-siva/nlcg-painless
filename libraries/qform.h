#include<stdio.h>

#ifndef _QFORM_H
	#define _QFORM_H
		
		struct qform {
			int row,col;
			int nx,ny;
			const char* matfile;
			const char* arrfile;
			double** A;
			double* b;
			double c;
			double ds;
			void (*calculate)(struct qform *this);
			void (*_vect_generate)(double* values,int width,int cur_col,double max,double ds,double** A,double* b,double c,FILE *fp,FILE *fp2);	
		};
		extern const struct qformClass {
			struct qform (*new)(const char* matfile,const char* arrfile,double c,int row,int col,int nx,int ny,double ds);
		}qform;

#endif


