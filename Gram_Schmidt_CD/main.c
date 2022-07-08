#include<stdio.h>
#include"matrix.h"
#include<math.h>
#include"gram_schmidt_cg.h"
#include"qform.h"
#include"jacobi_iter.h"

int main()
{
	const char* matfile = "mat.input";
	const char* arrfile = "arr.input";
	const char* initfile = "init.input";
	const char* ufile = "u.input";

	int row = 2 ,  col = 2;
	double errlimit = 1e-7;

	double c = 0;
	double ds = 0.1;
	int nx=5,ny=5;

	struct gram_schmidt_cg g = gram_schmidt_cg.new(matfile,arrfile,initfile,ufile,row,col);
	struct jacobi_iter j =  jacobi_iter.new(matfile,arrfile,initfile,row,col,errlimit);	
	struct qform q = qform.new(matfile,arrfile,c,row,col,nx,ny,ds);

	g.calculate(&g);	
	j.calculate(&j);	
	q.calculate(&q);	
}

	
	
		

