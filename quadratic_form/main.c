#include "qform.h"
#include<stdio.h>

//Program to solve quadratic function f = x^t A x - bx + c at all points of a given space and fine the minimum value

int main()
{
	const char* matfile = "mat.txt"; //The value of transformation matrix A	
	const char* arrfile = "arr.txt";//The value of array b
	double c = 0; //Constant value c
	int nx = 5 , ny = 5; //Computational Box Size (Make sure they are same)
	int row = 2, col = 2; //Space Dimensions (Make sure they are same)
	double ds = 0.1; //Space step

	struct qform q = qform.new(matfile,arrfile,c,row,col,nx,ny,ds); //Initializing the input values
	q.calculate(&q); //Calculating and writing the f-values in file "qfoutput.txt"	
}
