#include<stdio.h>
#include"qform.h"
#include"cg_method.h"
#include"gram_schmidt_cg.h"
int main()
{

	const char* matfile="mat.input";
	const char* arrfile="arr.input";
	const char* initfile="init.input";
	const char* ufile="u.input";
	double errlimit = 1e-7;
	int row=2,col=2;
	double const_c = 0;
	double ds = 0.1;
	int nx = 5,ny = 5;

	
	struct cg_method c = cg_method.new(matfile,arrfile,initfile,row,col,errlimit);
	struct gram_schmidt_cg g = gram_schmidt_cg.new(matfile,arrfile,initfile,ufile,row,col);
	struct qform q = qform.new(matfile,arrfile,const_c,row,col,nx,ny,ds);

	q.calculate(&q);
	g.calculate(&g);
	c.calculate(&c);
	
}
