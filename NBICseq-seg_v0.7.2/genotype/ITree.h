/*Written by Ruibin Xi, March 23 2011*/

#ifndef ITREE_XI
#define ITREE_XI

#include <stdlib.h>
#include <stdio.h>
#include "../lib/rbtree.h"


typedef double BElem; /*To make the code reusable, I define a new data type BElem (basic element) as int, later this maybe changed to double or float*/


typedef struct ITreeInterval_t{
		BElem l, r; /*left and right end points of the interval*/
		BElem max; /* the maximum of the end point of the intervals in the subtree of the node corresponding to this interval*/
		void *pt; /*A pointer that users can use to link the interval to other data types*/
		} ITreeInterval;

typedef struct ITreeIntervalSet_t{
		ITreeInterval **Is;
		int size, ll;
		} ITreeIntervalSet; /*Auxilary variable used for store the result of a query*/

typedef struct ITree_t{
		rbtree rbtr; /*the red-black tree*/
		ITreeInterval *Is; /*All intervals in the tree; once the tree is constructed, this should not be changed*/
		int size_Is, ll_Is; /*the size and length of Is and ll_endpt*/
		} *ITree;


int cmpInterval_left(void *a, void *b);/*compare two intervals (a and b should be of type ITreeInterval *)*/

ITreeIntervalSet * ITreeIntervalSet_create();
void ITreeIntervalSet_destroy(ITreeIntervalSet *);
ITree ITree_create(int size); /*create a new tree with specified size (allocate memory of Is)*/
int ITree_addInterval(ITree tr, BElem l, BElem r, void *pt); 

void ITreeQuery(ITree tr, BElem l, BElem r, ITreeIntervalSet *Iset); 
/* Given a interval [l,r], find all intervals that overalap with the given interval
   The query result will be stored in Iset; make sure initialize Iset before calling this function;
 */ 


void ITreeConstruct(ITree tr); /*construct an interval tree from its intervals Is; tr should be initialized.*/

void ITreeDestroy(ITree tr); 
/*The user is responsible to free the memory to pt of ITreeInterval;
  Other than the memory to pt of ITreeInterval, this function will free all memory assigned to the interval tree (including tr itselt)*/
#endif



