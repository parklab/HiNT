CC = gcc
CFLAGS = -g -O2 -lm -Wall

all: countRead genotype

countRead: countRead.o read.o ITree.o rbtree.o
	$(CC) $(CFLAGS) countRead.o read.o ITree.o rbtree.o -o countRead

genotype: genotype.o read.o gamma.o
	$(CC) $(CFLAGS) genotype.o read.o gamma.o -o genotype

genotype.o: read.h
	$(CC) $(CFLAGS) -c genotype.c read.c

countRead.o: read.h
	$(CC) $(CFLAGS) -c countRead.c read.c ITree.c rbtree.c

read.o: read.h
	$(CC) $(CFLAGS) -c read.c

ITree.o: ITree.h
	$(CC) $(CFLAGS) -c ITree.c rbtree.c

rbtree.o: rbtree.h
	$(CC) $(CFLAGS) -c rbtree.c

gamma.o: gamma.h
	$(CC) $(CFLAGS) -c gamma.c

clean:
	rm -rf *.o

