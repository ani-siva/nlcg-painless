#include"nlcg.h"
#include"matrix.h"
#include<stdio.h>
#include<string.h>
#include<math.h>

int main()
{
	const char* matfile="mat.input";
	const char* arrfile="arr.input";
	const char* initfile="init.input";
	const char* preconfile="precon.input";
	int row=2,col=2;
	int nx=5,ny=5;
	double c=0;
	double ds=1e-2;
	double errlimit=1e-2;
	
	struct nlcg n = nlcg.new(matfile,arrfile,initfile,preconfile,row,col,errlimit,ds);	
	n.calculate(&n);

}
