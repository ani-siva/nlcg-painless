#include"nlcg.h"
#include<stdio.h>
#include<math.h>
#include"matrix.h"
#include<stdbool.h>
#include<stdlib.h>
static void calculate(struct nlcg *this)
{


	double func(double** A,double* b,double c,double* x,int row,int col)
	{
		double f;		
		f = 0.5*vect_prod(mat_vect(A,x,row,col,'c'),x,row)-vect_prod(b,x,row)+c;
		return f;
	}

	double* fd(double** A,double* b,double c,double* x,int row,int col,double ds)
	{
		double* fp = init_vect(row);
		double* newx = init_vect(row);
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<row;j++)
				newx[j]=x[j];
			newx[i]=x[i]+ds;
			fp[i] = -1*(func(A,b,c,newx,row,col)-func(A,b,c,x,row,col))/ds;
		}
		return fp;
	}

	double** hessian(double** A,double* b,double c,double* x,int row,int col,double ds)
	{
		double** hess = init_mat(row,col);
		double* newx = init_vect(row);	
		
		double* newx2 = init_vect(row);
		double* fp = fd(A,b,c,x,row,col,ds);
		for(int i=0;i<row;i++)
		{
			for(int k=0;k<row;k++)
			{
				newx[k]=x[k];
				newx2[k]=x[k];
			}
			
			newx[i]=x[i]+ds;
			newx2[i]=x[i]-ds;
			for(int j=0;j<row;j++)
			{
				if(i!=j)
				{		
					newx[j]=x[j]+ds;
					double d1 = func(A,b,c,newx,row,col)/(ds*ds);
					double d2=0;
					for(int k=0;k<row;k++)
						d2 += (fp[k]/ds);
					double d3 = func(A,b,c,x,row,col)/(ds*ds);
					hess[i][j] = d1+d2-d3;
				}
				else
				{
					
					for(int k=0;k<row;k++)
					{
						newx[k]=x[k];
						newx2[k]=x[k];
					}
						newx[i]=x[i]+ds;
						newx2[i]=x[i]-ds;
					hess[i][j]=(func(A,b,c,newx,row,col)-2*(func(A,b,c,x,row,col))+func(A,b,c,newx2,row,col))/(ds*ds);
				}
			}
		}
		return hess;
	}
	double** A = read_mat(this->row,this->col,this->matfile);
	double* b = read_vect(this->row,this->arrfile);	
	double* x = read_vect(this->row,this->initfile);
	double** M = read_mat(this->row,this->col,this->preconfile);	
	double c = 0;
	
	double** h = hessian(A,b,c,x,this->row,this->col,this->sstep);

	int i=0,k=0;
	int imax = 2,jmax=2;
	double* r = fd(A,b,c,x,this->row,this->col,this->sstep);
	double* d = fd(A,b,c,x,this->row,this->col,this->sstep);
	
	double* fp = fd(A,b,c,x,this->row,this->col,this->sstep);
	for(int i=0;i<this->row;i++)
		fp[i] = -1*fp[i];

	double delta_new = vect_prod(r,r,this->row);
	double delta_0 = delta_new;
	double* newx = init_vect(this->row);
	double* newd = init_vect(this->row);
	print_vect(x,this->row);
	printf("\n");
	while(i<imax && delta_new > pow(this->errlimit,2)*delta_0)
	{
		int j=0;
		double alpha;
		
		double delta_d = vect_prod(d,d,this->row);
		do
		{
			alpha = -1*(vect_prod(fp,d,this->row))/(vect_prod(mat_vect(h,d,this->row,this->col,'c'),d,this->row));
			for(int i=0;i<this->row;i++)
				newx[i] = x[i]+alpha*d[i];
			j++;
		}while(j<jmax && pow(alpha,2)*delta_d > pow(this->errlimit,2));
		
		print_vect(newx,this->row);
		printf("\n");
		double* newr = fd(A,b,c,newx,this->row,this->col,this->sstep);
		double delta_old = delta_new;
		delta_new = vect_prod(newr,newr,this->row);	
		double beta = delta_new/delta_old;
		for(int i=0;i<this->row;i++)
			newd[i] = newr[i] + beta * d[i]; 
		k++;
		if(k==10 || vect_prod(newr,d,this->row)<=0)
		{
			for(int i=0;i<this->row;i++)
				newd[i] = newr[i];
			k=0;
		}
		i++;
		for(int i=0;i<this->row;i++)
		{
			d[i]=newd[i];
			r[i]=newr[i];
			x[i]=newx[i];
		}
	}
}

static struct nlcg new(const char* matfile,const char* arrfile,const char* initfile,const char* preconfile,int row,int col,double errlimit,double sstep)
{
	return (struct nlcg){.matfile=matfile,.arrfile=arrfile,.initfile=initfile,.preconfile=preconfile,.row=row,.col=col,.errlimit=errlimit,.sstep=sstep,.calculate=&calculate};
}

const struct nlcgClass nlcg={.new=&new};
