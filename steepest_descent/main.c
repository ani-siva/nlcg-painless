#include<stdio.h>
#include<math.h>
#include"matrix.h"
#include"steep.h"
#include"qform.h"

int main(void)
{

	const char* matfile = "mat.input";
	const char* arrfile = "arr.input";
	const char* initfile = "init.input";
	double errlimit = 1e-7;
	int nx = 5, ny = 5;
	int row = 2, col = 2;
	double ds = 0.1;
	double c = 0.0;
	struct steep s = steep.new(matfile,arrfile,initfile,row,col,errlimit);
	s.calculate(&s);

	struct qform q = qform.new(matfile,arrfile,c,row,col,nx,ny,ds);
	q.calculate(&q);	
	return 0;
}

