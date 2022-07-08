FILENAME=main
LIBBPATH = lib
#LIBBPATH=$(HOME)/Documents/Work_Bench/Codes/NCG_Mine/libraries
LIBSOURCES= nlcg.o  matrix.o  
FLAGS=-lm
run : libs.o $(FILENAME).o exec.o 
libs.o : $(LIBSOURCES) 

nlcg.o:
	gcc -c $(LIBBPATH)/nlcg.c
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
