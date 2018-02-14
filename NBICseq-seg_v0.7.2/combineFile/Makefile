CC = gcc
CFLAGS = -g -O2 -lm -Wall

combineFile: combine.o read.o
	$(CC) $(CFLAGS) combine.o read.o -o combineFile

combine.o: read.h
	$(CC) $(CFLAGS) -c combine.c read.c

read.o: read.h
	$(CC) $(CFLAGS) -c read.c

clean:
	rm -rf *.o

