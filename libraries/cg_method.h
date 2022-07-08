#include<stdio.h>
 #ifndef _CG_METHOD_H
	#define _CG_METHOD_H
	
		struct cg_method{
			int row,col;
			const char* matfile;
			const char* arrfile;
			const char* initfile;

			double** A;
			double* b;
			double* x;
			double** u;
			double errlimit;
			void(*calculate)(struct cg_method *this);

		};
	
		extern const struct cg_methodClass{
			struct cg_method(*new)(const char* matfile,const char* arrfile,const char* initfile,int row,int col,double errlimit);
		}cg_method;

#endif
