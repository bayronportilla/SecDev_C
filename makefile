############################################################
#
# variables
#
############################################################

CC       = gcc
CFLAGS   = -I/usr/local/include
CCFLAGS  = -L/usr/local/lib
TARGET   = main 

#OBJS = FetchInfo.o ModQuad.o Params.o Units.o \
#allvars.o bulk.o interpol.o main.o matrix.o

OBJGS = main generic

#LIBS = -lgsl -lgslcblas -lm -lconfig


main:main.o generic.o matrix.o Units.o bulk.o Params.o FetchInfo.o interpol.o ModQuad.o
	$(CC) main.o generic.o matrix.o Units.o bulk.o Params.o FetchInfo.o interpol.o ModQuad.o -lgsl -lgslcblas -lm -lconfig -o main.x
	rm *.o



main.o:main.c
	gcc -c main.c

generic.o:generic.c
	gcc -c generic.c

#matrix.o:matrix.c
#	gcc -c matrix.c

#Units.o:Units.c
#	gcc -c Units.c

#bulk.o:bulk.c
#	gcc -c bulk.c



clean:
	rm *~
	rm *.x

run:
	./main.x


