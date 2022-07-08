FILENAME=main
LIBRARY=matrix
FLASGS=-lm
LIBBPATH=$(HOME)/Documents/Work_Bench/Codes/NCG_Mine/libraries

run : $(LIBRARY).o $(FILENAME).o qform.o exec.o 

qform.o:
	gcc -c $(LIBBPATH)/qform.c
$(LIBRARY).o:
	gcc -c $(LIBBPATH)/$(LIBRARY).c
$(FILENAME).o:
	gcc -I$(LIBBPATH) -c $(FILENAME).c
exec.o:	
	gcc -o exec.out $(LIBRARY).o $(FILENAME).o $(FLAGS)qform.o
.PHONY : clean
clean:
	rm *.out
	rm *.o
