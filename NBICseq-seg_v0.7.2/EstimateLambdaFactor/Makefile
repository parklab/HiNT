CC = gcc
CFLAGS = -g -O2 -lm -Wall

EstLambdaFct: EstLambdaFct.o read.o
	$(CC) $(CFLAGS) EstLambdaFct.o read.o -o EstLambdaFct

EstLambdaFct.o: read.h
	$(CC) $(CFLAGS) -c EstLambdaFct.c read.c

read.o: read.h
	$(CC) $(CFLAGS) -c read.c

clean:
	rm -rf *.o

