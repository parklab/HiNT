CC = gcc
CFLAGS = -g -O2 -lm -Wall

bootstrap: bootstrap.o read.o gamma.o
	$(CC) $(CFLAGS) bootstrap.o read.o gamma.o -o bootstrapTest

bootstrap.o: read.h
	$(CC) $(CFLAGS) -c bootstrap.c read.c

read.o: read.h
	$(CC) $(CFLAGS) -c read.c

gamma.o: gamma.h
	$(CC) $(CFLAGS) -c gamma.c

clean:
	rm -rf *.o

