#include"jacobi_iter.h"
#include"matrix.h"
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>

bool check_err(double* newx,double* x,int row,double errlimit)
{
	double* err = init_vect(row);
	for(int i=0;i<row;i++)
	{	
		err[i] = fabs(newx[i]-x[i]);
		if(err[i] > errlimit)
			return false;
	}

	return true;
}
		
static void calculate(struct jacobi_iter *this)
{
	this->A = read_mat(this->row,this->col,this->matfile);
	this->b = read_vect(this->row,this->arrfile);
	this->x = read_vect(this->row,this->initfile);

	FILE *fp;
	fp = fopen("jacobi-iter-coords.txt","w");

	double** D = init_mat(this->row,this->col);
	double** E = init_mat(this->row,this->col);
		
	for(int i=0;i<this->row;i++)
	{
		for(int j=0;j<this->col;j++)
		{
			if(i==j)
				D[i][j] = this->A[i][j];
			if(i!=j)
				E[i][j] = this->A[i][j];
		}
	}

	double** Dinv = lu_inverse_mat(D,this->row,this->col);
	double** B = mat_prod(Dinv,E,this->row,this->col);
	double* z = mat_vect(Dinv,this->b,this->row,this->col,'r');
	while(true)
	{
		double* Bx = mat_vect(B,this->x,this->row,this->col,'r');
		
		write_vect(this->x,this->row,fp);
		fprintf(fp,"\n");
		double* newx = init_vect(this->row);
		for(int i=0;i<this->row;i++)
		{
			newx[i] = -1*Bx[i] + z[i];
		}
		
		if(check_err(newx,this->x,this->row,this->errlimit))
			break;		
		
		for(int i=0;i<this->row;i++)
		{
			this->x[i] = newx[i];
		}

		//print_vect(newx,this->row);
		//printf("\n");
		free(newx);
		free(Bx);
	}
	fclose(fp);
}

static struct jacobi_iter new(const char* matfile,const char* arrfile,const char* initfile,int row,int col,double errlimit)
{

	return (struct jacobi_iter){.matfile=matfile,.arrfile=arrfile,.initfile=initfile,.row=row,.col=col,.errlimit=errlimit,.calculate=&calculate};

}

const struct jacobi_iterClass jacobi_iter={.new=&new};


