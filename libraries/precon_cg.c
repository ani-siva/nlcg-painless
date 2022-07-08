#include<stdio.h>
#include"precon_cg.h"
#include"matrix.h"
#include<stdbool.h>

static void calculate(struct precon_cg *this)
{

	bool check_error(double* newx,double* x,int row,double errlimit)
	{
		for(int i=0;i<row;i++)
			if((newx[i]-x[i])>errlimit) return false;
		return true;
	}
	double alpha(double* r,double** Minv, double* d, double** A ,int row,int col)
	{
		double num = vect_prod(r,mat_vect(Minv,r,row,col,'c'),row);
		double denom = vect_prod(d,mat_vect(A,d,row,col,'c'),row);
		return num/denom;
	}
	double beta(double* newr,double** Minv,double* r,int row,int col)
	{
		double num = vect_prod(newr,mat_vect(Minv,newr,row,col,'c'),row);
		double denom = vect_prod(r,mat_vect(Minv,r,row,col,'c'),row);
		return num/denom;
	}
		
	this->A=read_mat(this->row,this->col,this->matfile);
	this->b=read_vect(this->row,this->arrfile);
	this->x=read_vect(this->row,this->initfile);
	this->M=read_mat(this->row,this->col,this->preconfile);
	double* r=init_vect(this->row);
	double* Ax=mat_vect(this->A,this->x,this->row,this->col,'c');
	double** Minv=lu_inverse_mat(this->M,this->row,this->col);
	this->MInvA = mat_prod(Minv,this->A,this->row,this->col);

	for(int i=0;i<this->row;i++)
		r[i]=this->b[i]-Ax[i];

	double* d = mat_vect(Minv,r,this->row,this->col,'c');
	double* newx = init_vect(this->row);
	double* newr = init_vect(this->row);
	double* newd = init_vect(this->row);
	while(true)
	{
		print_vect(this->x,this->row);
		printf("\n");
		double alp = alpha(r,Minv,d,this->A,this->row,this->col);
		double* Ad = mat_vect(this->A,d,this->row,this->col,'c');
		for(int i=0;i<this->row;i++)
		{
			newx[i] = this->x[i] + alp * d[i];
			newr[i] = r[i] - alp*Ad[i];	
		}
		if(check_error(newx,this->x,this->row,this->errlimit)) break;
		double bet = beta(newr,Minv,r,this->row,this->col);
		double* Minvr = mat_vect(Minv,newr,this->row,this->col,'c');
		for(int i=0;i<this->row;i++)
			newd[i] = Minvr[i] + bet * d[i];
		for(int k=0;k<this->row;k++)
		{
			this->x[k] = newx[k];	
			r[k] = newr[k];
			d[k] = newd[k];	
		}
	}
}

static struct precon_cg new(const char* matfile,const char* arrfile,const char* initfile,const char* preconfile, int row,int col,double errlimit)
{
	return (struct precon_cg){.matfile=matfile,.arrfile=arrfile,.initfile=initfile,.preconfile=preconfile,.row=row,.col=col,.errlimit=errlimit,.calculate=&calculate};
}

const struct precon_cgClass precon_cg={.new=&new};

