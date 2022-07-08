#include "qform.h"
#include<stdio.h>
#include "matrix.h"
#include<stdlib.h>
#include<string.h>

void _vect_generate(double *values,int width,int cur_col,double max,double ds,double** A,double* b,double c,FILE* fp,FILE *fp2)
{
	if(cur_col == width)
	{

		double xtax,bx,f;
	        double* xtA = mat_vect(A,values,width,width,'c');
		xtax = vect_prod(xtA,values,width);
		bx = vect_prod(b,values,width);
		f = 0.5*xtax-bx+c;
		for(int i=0;i<width;i++)
		{
			fprintf(fp2,"%lf ",xtA[i]-b[i]);
		}
		
		
		
		fprintf(fp,"%lf ",f);
		for(int i=0;i<width;i++)
		{
			fprintf(fp,"%lf%c",values[i],(i<width-1)?' ':'\n');
			fprintf(fp2,"%lf%c",values[i],(i<width-1)?' ':'\n');
		}
		free(xtA);
	}
	else
	{
		double i = -1*max;
		while(i<=max)
		{
			values[cur_col]=i;
			_vect_generate(values,width,cur_col+1,max,ds,A,b,c,fp,fp2);
			i = i + ds;
		}
	}
}
static void calculate(struct qform *this)
{
	this->A=read_mat(this->row,this->col,this->matfile);	
	this->b=read_vect(this->row,this->arrfile);		
  	FILE *fp,*fp2;
	fp = fopen("qf-output.txt","w");	
	fp2 = fopen("qf-vect-output.txt","w");
	double start = -1*this->nx;
	double values[this->row] ;
	memset(values,start,this->row*sizeof(double));
	_vect_generate(values,this->row,0,this->nx,this->ds,this->A,this->b,this->c,fp,fp2);
	fclose(fp);
	fclose(fp2);
}

static struct qform new(const char* matfile,const char* arrfile,double c,int row,int col,int nx,int ny,double ds){
	return (struct qform){.matfile=matfile,.arrfile=arrfile,.c=c,.row=row,.col=col,.nx=nx,.ny=ny,.ds=ds,.calculate=&calculate};
}


const struct qformClass qform={.new=&new};

