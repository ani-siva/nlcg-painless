#include"cg_method.h"
#include<stdlib.h>
#include<stdio.h>
#include"matrix.h"
#include<math.h>
#include<stdbool.h>


static void calculate(struct cg_method *this)
{

	double alpha(double* r,double* d,double** A,int row,int col)
	{
		double num = vect_prod(r,r,row);
		double denom = vect_prod(mat_vect(A,d,row,col,'c'),d,row);
		return num/denom;
	}
	double beta(double* newr,double* r,int row)
	{
		double num = vect_prod(newr,newr,row);
		double denom = vect_prod(r,r,row);
		return num/denom;
	}
	bool check_error(double* newx,double* x,double errlimit,int row)
	{
		for(int i=0;i<this->row;i++)
			if((newx[i]-x[i])>errlimit) return false;
		return true;
	}
	this->A = read_mat(this->row,this->col,this->matfile);
	this->b = read_vect(this->row,this->arrfile);
	this->x = read_vect(this->row,this->initfile);

	double* d = init_vect(this->row);
	double* r = init_vect(this->row);
	double* Ax = mat_vect(this->A,this->x,this->row,this->col,'c');
	for(int i=0;i<this->row;i++)
	{
		r[i] = this->b[i] - Ax[i];
		d[i] = this->b[i] - Ax[i];
	}
	free(Ax);
	free(this->b);
	FILE *fp;
	fp = fopen("cg-method-coords.txt","w");
	while(true)
	{
		write_vect(this->x,this->row,fp);
		fprintf(fp,"\n");
		double alp = alpha(r,d,this->A,this->row,this->col);
		double* newx = init_vect(this->row);
		double* newr = init_vect(this->row);
		double* newd = init_vect(this->row);
		double* Ad = mat_vect(this->A,d,this->row,this->col,'c');
		for(int j=0;j<this->row;j++)
		{
			newx[j] = this->x[j] + alp * d[j];	
			newr[j] = r[j] - alp * Ad[j];
		}

		if(check_error(newx,this->x,this->errlimit,this->row)) break;
		
		double bet = beta(newr,r,this->row);
		for(int j=0;j<this->row;j++)
			newd[j] = newr[j] + bet * d[j];
		for(int k=0;k<this->row;k++)
		{
			this->x[k] = newx[k];
			r[k] = newr[k];
			d[k] = newd[k];
		}
		free(newx);
		free(newr);
		free(newd);
		free(Ad);
	}
	fclose(fp);
}

static struct cg_method new(const char* matfile,const char* arrfile,const char* initfile,int row,int col,double errlimit)
{
	return (struct cg_method){.matfile=matfile,.arrfile=arrfile,.initfile=initfile,.row=row,.col=col,.errlimit=errlimit,.calculate=&calculate};
}

const struct cg_methodClass cg_method={.new=&new};

