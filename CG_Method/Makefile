FILENAME=main
LIBBPATH=$(HOME)/Documents/Work_Bench/Codes/NCG_Mine/libraries
LIBSOURCES= cg_method.o qform.o matrix.o gram_schmidt_cg.o 
FLAGS=-lm
run : libs.o $(FILENAME).o exec.o 
libs.o : $(LIBSOURCES) 

cg_method.o:
	gcc -c $(LIBBPATH)/cg_method.c
steep.o:
	gcc -c $(LIBBPATH)/steep.c
jacobi_iter.o:
	gcc -c $(LIBBPATH)/jacobi_iter.c
gram_schmidt_cg.o:
	gcc -c $(LIBBPATH)/gram_schmidt_cg.c
qform.o:
	gcc -c $(LIBBPATH)/qform.c
matrix.o:
	gcc -c $(LIBBPATH)/matrix.c 
$(FILENAME).o:
	gcc -I$(LIBBPATH) -c $(FILENAME).c
exec.o:	
	gcc -o exec.out  $(FILENAME).o $(LIBSOURCES) $(FLAGS) 
	rm *.o
.PHONY : clean full-clean

full-clean: clean 
	rm *.txt

clean:
	rm *.out
