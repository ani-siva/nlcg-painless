#include"gram_schmidt_cg.h"
#include<stdlib.h>
#include<stdio.h>
#include"matrix.h"
#include<math.h>

double** gram_schmidt_orthogonal(double** mat,double** u, int row,int col)
{

	double** gso = init_mat(row,col);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
			gso[i][j]=u[i][j];
	}
	void normalize(double *x)
	{
		double norm=sqrt(vect_prod(x,x,row));
		for(int i=0;i<row;i++)
			x[i]/=norm;
	}
	       		
	double scale(double* x1 ,double** mat,double* x2,int row,int col)
	{	
		double num = vect_prod(mat_vect(mat,x1,row,col,'r'),x2,row);
		double denom = vect_prod(mat_vect(mat,x2,row,col,'r'),x2,row);
		return num/denom;
	}
	for(int i=1;i<row;++i)
	{
		for(int j=0;j<i;++j)
		{
			double sc  = scale(gso[i],mat,gso[j],row,col);
			for(int k=0;k<col;++k)
				gso[i][k] -= sc * gso[j][k];
		}
	}
	
		
	return gso;
}


static void calculate(struct gram_schmidt_cg *this)
{
	double* residual(double** A,double* b,double* x,int row,int col)
	{
		double* res = init_vect(row);
		double* Ax = mat_vect(A,x,row,col,'c');
		for(int i=0;i<row;i++)
			res[i] = b[i]-Ax[i];
		return res;
	}
	double alpha(double* res,double* d,double** A,int row,int col)
	{
		double alpha;
		double num = vect_prod(d,res,row);
		double denom = vect_prod(mat_vect(A,d,row,col,'c'),d,row);
		return num/denom;
	}
	this->A=read_mat(this->row,this->col,this->matfile);
	this->b=read_vect(this->row,this->arrfile);
	this->x=read_vect(this->row,this->initfile);
	this->u=read_mat(this->row,this->col,this->ufile);
	
	
	double** v = gram_schmidt_orthogonal(this->A,this->u,this->row, this->col);
	
	FILE *fp ;
	fp = fopen("gso-coords.txt","w");	

	for(int i=0;i<this->row;++i)
	{
		double alp = 0.0;
		double* res = residual(this->A,this->b,this->x,this->row,this->col);
		alp = alpha(res,v[i],this->A,this->row,this->col);
		double* newx = init_vect(this->row);
		for(int j=0;j<this->row;j++)
		{
			newx[j] = this->x[j] + alp*v[i][j];
		}
		write_vect(this->x,this->row,fp);		
		fprintf(fp,"\n");	
		
		for(int i=0;i<this->row;i++)
		{
			this->x[i] = newx[i];
		}

	}


		write_vect(this->x,this->row,fp);		
		fprintf(fp,"\n");
		fclose(fp);	
}

static struct gram_schmidt_cg new(const char* matfile,const char* arrfile,const char* initfile,const char* ufile,int row,int col)
{
	return (struct gram_schmidt_cg){.matfile=matfile,.arrfile=arrfile,.initfile=initfile,.ufile=ufile,.row=row,.col=col,.calculate=&calculate};

}

const struct gram_schmidt_cgClass gram_schmidt_cg={.new=&new};

