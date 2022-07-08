#include<stdio.h>
#include"jacobi_iter.h"
#include"qform.h"
#include"steep.h"

int main(void)
{
	const char* matfile = "mat.input";
	const char* arrfile = "arr.input";
	const char* initfile = "init.input";
	int row = 2 ; int col = 2;
	double errlimit = 1e-7;

	
	double c = 0;
	int nx = 5 , ny = 5;
	double ds = 0.1;

	struct jacobi_iter j = jacobi_iter.new(matfile,arrfile,initfile,row,col,errlimit);
	struct steep s = steep.new(matfile,arrfile,initfile,row,col,errlimit);
	struct qform q = qform.new(matfile,arrfile,c,row,col,nx,ny,ds); 	
	j.calculate(&j);	
	s.calculate(&s);
	q.calculate(&q);	

}
