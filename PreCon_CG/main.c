#include<stdio.h>
#include"precon_cg.h"

int main()
{
	const char* matfile="mat.input";
	const char* arrfile="arr.input";
	const char* initfile="init.input";
	const char* preconfile="precon.input";
	int row = 2, col = 2;
	double errlimit = 1e-6;
	struct precon_cg p = precon_cg.new(matfile,arrfile,initfile,preconfile,row,col,errlimit);
	p.calculate(&p);
}
