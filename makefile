CC = gcc
CFLAGS = -Wall -I/usr/local/inlcude
CCFLAGS = -L/usr/local/lib
TARGET = main


$(TARGET):$(TARGET).c
	$(CC) $(CFLAGS) -c $@.c
	$(CC) $(CCFLAGS) $@.o -lgsl -lgslcblas -lm -o $@.x


clean:
	rm -rf *#
	rm -rf *~

run:
	./main.x
