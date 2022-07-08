#include<stdio.h>
#ifndef _NLCG_H
	#define _NLCG_H
		struct nlcg{
			int row,col;
			const char* matfile;
			const char* arrfile;
			const char* initfile;
			const char* preconfile;

			double errlimit;
			double sstep;
			void(*calculate)(struct nlcg *new);
		};

		extern const struct nlcgClass{
			struct nlcg(*new)(const char* matfile,const char* arrfile,const char* initfile,const char* preconfile,int row,int col,double errlimit,double sstep);
		}nlcg;
#endif
