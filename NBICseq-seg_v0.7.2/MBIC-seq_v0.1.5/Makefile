CC = gcc
CFLAGS = -g -O2 -lm -Wall


MBICseq: bin.o mbic-seq.o rbtree.o read.o
	$(CC) $(CFLAGS) bin.o mbic-seq.o rbtree.o read.o -o MBICseq

mbic-seq.o: bin.h rbtree.h read.h
	$(CC) $(CFLAGS) -c mbic-seq.c bin.c read.c rbtree.c

bin.o: bin.h rbtree.h read.h
	$(CC) $(CFLAGS) -c bin.c read.c rbtree.c

rbtree.o: rbtree.h
	$(CC) $(CFLAGS) -c rbtree.c

read.o: read.h
	$(CC) $(CFLAGS) -c read.c

clean:
	rm -rf *.o
