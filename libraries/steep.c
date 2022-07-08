#include"steep.h"
#include"matrix.h"
#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>

double* xupdate(double* x0,double alpha,double* res,int row)
{
	double *x;
	x = (double *)malloc(sizeof(double)*row);
	for(int i=0;i<row;i++)
	{
		x[i] = x0[i] + (alpha * res[i]);
	}
	return x;
}

double alpha(double* res,double** A,int row,int col)
{
	double alpha;
	double* ra = mat_vect(A,res,row,col,'r');
	alpha = (vect_prod(res,res,row))/(vect_prod(ra,res,row));
	return alpha;
}
double* residual(double** A,double* b,double* x,int row,int col)
{
	double* res;
	res = (double *)malloc(sizeof(double)*row);
	double* Ax = mat_vect(A,x,row,col,'c');
	for(int i=0; i<row; i++)
	{
		res[i] = b[i] - Ax[i];
	}
	return res;
}

int check_error(double* xi,double* x,double errlimit,int row)
{
	double *err;
	int check_flag = 0;
	err = (double *)malloc(sizeof(double)*row);
	for(int i=0;i<row;i++)
	{
		err[i] = xi[i] - x[i];
		if(err[i] > errlimit)
		{
			check_flag = 1;
		}
	}
	return check_flag;
}

static void calculate(struct steep *this)
{
	this->A = read_mat(this->row,this->col,this->matfile);
	this->b = read_vect(this->row,this->arrfile);
	this->x = read_vect(this->row,this->initfile);
	FILE *fp;
	fp = fopen("steep-coords.txt","w");
	while(true)
	{
		double alp = 0.0;
		double* res = residual(this->A,this->b,this->x,this->row,this->col);
		alp = alpha(res,this->A,this->row,this->col);
		double* newx = xupdate(this->x,alp,res,this->row);

		if(check_error(newx,this->x,this->errlimit,this->row) == 0) break;
		write_vect(this->x,this->row,fp);	
		fprintf(fp,"\n");

		for(int i=0;i<this->row;i++)
		{
			this->x[i] = newx[i];
		}
	}

	fclose(fp);

}

static struct steep new(const char* matfile, const char* arrfile,const char* initfile, int row, int col, double errlimit)
{
	return (struct steep){.matfile=matfile, .arrfile=arrfile, .initfile=initfile, .row=row, .col=col, .errlimit=errlimit, .calculate=&calculate};
}

const struct steepClass steep={.new=&new};
