#include "matrix.h"
#include<stddef.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

//Initializes a matrix of given dimensions rxc
double** init_mat(int r,int c)
{
	
	double **mat;
	mat = malloc(r*sizeof(double *));

	for(int i=0;i<r;i++)
	{
		mat[i] = malloc(c*sizeof(double));
	}
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			mat[i][j]=0;
		}
	}
	return mat;
}


//reads matrix values from a text file and stores it in memory.
double** read_mat(size_t r,size_t c,const char* filename)
{
	double** mat = init_mat(r,c);
	FILE *pf;
	pf = fopen(filename,"rt");
	if(pf == NULL)
		return 0;
	for(size_t i=0;i<r;++i)
	{
		for(size_t j=0;j<c;j++)
			fscanf(pf,"%lf",mat[i]+j);
	}
	fclose(pf);
	return mat;
}

//prints a matrix onto the shell
void print_mat(double** A,int r,int c)
{
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			printf("%lf ",A[i][j]);
		}
		printf("\n");
	}

}

//inverts a matrix using the LU decomposition method

double** lu_inverse_mat(double** A,int r,int c)
{
	
	double l[r][c], u[r][c] ; //initializing the upper and lower triangular matrices
	memset(l,0,sizeof(l));
	memset(u,0,sizeof(u));

	double** x = init_mat(r,c);
	//Decomposing matrix A into upper and lower triangular matrices//
	for(int i=0;i<r;i++) 
	{
		for(int k=i;k<c;k++)  //Reading through each column and finding the sum
		{
			double sum=0;
			for(int j=0;j<i;j++)
				sum += (l[i][j]*u[j][k]);
			u[i][k] = A[i][k] - sum; //Calculating the upper triangular matrix
		}
		for(int k=i;k<c;k++)
		{
			if(i==k)
				l[i][i] = 1;
			else
			{
				double sum = 0;
				for(int j=0;j<i;j++)
					sum += (l[k][j] * u[j][i]);
				l[k][i] = (A[k][i]-sum)/u[i][i]; //Calculating lower triangular matrix
			}
		}
	}
		

	//Decomposition Done//
	
	double I[r][c]; //Identitiy Matrix
	double z[r]; //Solution of L.z = b, which will be used to solve for the inverse matrix x, by U.x = z

	//Initializng the identity matrix	
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			if(i==j)
				I[i][j]=1;
			
			else
				I[i][j]=0;
		}
	}
	int k=0;
	while(k<c) //Solving for x column wise
	{
		memset(z,0,sizeof(z)); //Restoring z value after every column of X is solved
		//Using forward substitution to solve L.z =b
		z[0]=I[k][0];
		for(int i=0;i<r;i++)
		{
			double sum=0;
			for(int j=0;j<=i-1;j++)
			{
				sum += l[i][j]*z[j];
				z[i] = (I[k][i]-sum)/l[i][i];
			}
		}
		//Using backward substitution to solve U.x = z
		x[k][r] = z[r];
		for(int i=r;i>=0;i--)
		{
			double sum=0;
			for(int j=i;j<r;j++)
			{
				sum += u[i][j]*x[k][j];
				x[k][i] = (z[i]-sum)/u[i][i];
			}
		}
			k++; //Moving to the next column

	}

	return x;
}

//Calculates the product of two matrices

double** mat_prod(double** mat1,double** mat2,int row,int col)
{
	double** prod = init_mat(row,col);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			for(int k=0;k<row;k++)
				prod[i][j] += mat2[i][k]*mat1[k][j];
		}
	}

	return prod;
}

//Performs Gram-Schmidt Orthogonalization on a set of chosen vectors and outputs the orthogonal vectors in the form of a matrix.

double** gram_schmidt_orthog(double** q,int row,int col)
{
	double** gso = init_mat(row,col);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
			gso[i][j] = q[i][j];
	}
	void normalize(double *x)
	{
		double norm = sqrt(vect_prod(x,x,row));
		for(int i=0;i<row;++i)
			x[i] /= norm;
	
	}
	for(int i=1 ; i<row ; ++i)
	{
		for(int j=0;j<i;++j)
		{
			double scale = vect_prod(gso[j],gso[i],row)/vect_prod(gso[j],gso[j],row);
			for(int k=0;k<col;++k)
				gso[i][k] -= scale*gso[j][k];
		}
	}
	for(int i=0;i<col;++i)
		normalize(gso[i]);
	return gso;
}

//Writes a matrix onto a file
void write_mat(double** A,int row,int col,FILE *fp)
{
	for(int i=0;i<row;i++)
	{	for(int j=0;j<col;j++)
			fprintf(fp,"%lf ",A[i][j]);
		fprintf(fp,"\n");
	}

}

//prints a vector onto the shell
void print_vect(double* vect,int r)
{
	for(int i=0;i<r;i++)
	{
		printf("%lf \n",vect[i]);
	}
}


//calculates the scalar product of two vectors
double vect_prod(double* vect1,double* vect2,int r)
{
	double prod=0.0;
	for(int i=0;i<r;i++)
	{
		prod += vect1[i]*vect2[i];
	}
	return prod;
}
//Initializes a vector of length r
double* init_vect(int r)
{
	double* vec = malloc(r*sizeof(double));
	for(int i=0;i<r;i++)
	{
		vec[i] = 0.0;
	}
	return vec;
}

//read vector from a file and store it in memory
double* read_vect(size_t r,const char* filename)
{
	double* vec = init_vect(r);
	FILE *pf;
	pf = fopen(filename,"rt");
	if(pf == NULL)
		return 0;
	for(size_t i=0;i<r;++i)
	{
		fscanf(pf,"%lf",&vec[i]);
	}
	fclose(pf);
	return vec;
}

//Writes a vector to a file

void write_vect(double* vect, int row, FILE *fp)
{
	for(int i =0 ;i<row; i++)
	{
		fprintf(fp,"%lf ",vect[i]);
	}
}



//Does a row-wise or column-wise product of a matrix and vector
double* mat_vect(double** mat,double* vect,int r,int c,char choice)
{
	double* prodr = malloc(r*sizeof(double));
	double* prodc = malloc(c*sizeof(double));
	if(choice=='r')
	{
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			prodr[i] += mat[j][i] * vect[j];
		}
	}
	return prodr;
	}
	else
	{
	for(int j=0;j<c;j++)
	{
		for(int i=0;i<r;i++)
		{
			prodc[j] += mat[i][j] * vect[i];
		}
	}
	return prodc;
	}
}	
