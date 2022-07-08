FILENAME=main
LIBBPATH=$(HOME)/Documents/Work_Bench/Codes/NCG_Mine/libraries
LIBSOURCES =  matrix.o jacobi_iter.o qform.o steep.o 

run : libs.o $(FILENAME).o exec.o 
libs.o : $(LIBSOURCES) 

jacobi_iter.o:
	gcc -c $(LIBBPATH)/jacobi_iter.c
qform.o:
	gcc -c $(LIBBPATH)/qform.c
matrix.o:
	gcc -c $(LIBBPATH)/matrix.c
steep.o:
	gcc -c $(LIBBPATH)/steep.c
$(FILENAME).o:
	gcc -I$(LIBBPATH) -c $(FILENAME).c
exec.o:	
	gcc -o exec.out  $(FILENAME).o $(LIBSOURCES) 

.PHONY : clean full-clean

full-clean: clean 
	rm *.txt

clean:
	rm *.o
	rm *.out
