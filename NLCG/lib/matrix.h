#include<stddef.h>
#include<stdio.h>

#ifndef _MATRIX_H
	#define _MATRIX_H

	double** init_mat(int r,int c);
	double** read_mat(size_t r,size_t c,const char* filename);
	void print_mat(double** A,int r,int c);
	double** lu_inverse_mat(double** A, int r, int c);
	double** mat_prod(double** mat1, double** mat2, int r,int c);
	double** gram_schmidt_orthog(double** mat, int r,int c);
	void write_mat(double** A,int r,int c,FILE *fp);	

	void print_vect(double* vect,int r);
	double vect_prod(double* vect1,double* vect2,int r);
	double* init_vect(int r);
	double* read_vect(size_t r,const char* filename);
	void write_vect(double* vect,int row,FILE *fp);

	double* mat_vect(double** mat,double* vect,int r,int c,char choice);
#endif
