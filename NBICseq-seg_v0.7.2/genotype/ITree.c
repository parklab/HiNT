#include "ITree.h"

static BElem find_max_subrbtree(rbtree_node node); 
/*Given a node in the ITree (used a rbtree data structure)
  find the maximum of the endpoints of the intervals rooted at the given node;
  Here I use recursion to achieve the goal. Later, I may want to use non-recursion method
*/

static void ITreeQuery_node(rbtree_node node, BElem l, BElem r, ITreeIntervalSet *Iset);
/* given a node node, find all intervals of the subtree rooted at node that overlap with the given interval [l, r]
   recursive method is used.
   here I assume r>=l, this funciton will NOT check if this is true
 */

ITree ITree_create(int size_in)
	{ int size;
	  ITree tr = malloc(sizeof(struct ITree_t));
	  if(tr==NULL){
		fprintf(stderr,"Error in ITree_create: memory allocation failed\n");
		exit(1);
		}

	  if(size_in>0) size = size_in+10; else size = 1010;

	  tr->Is = (ITreeInterval *) malloc(sizeof(ITreeInterval)*size);
	  tr->rbtr = rbtree_create();

	  if(tr->Is==NULL||tr->rbtr==NULL){
                fprintf(stderr,"Error in ITree_create: memory allocation failed\n");
                exit(1);		
		}

	  tr->size_Is = size;
	  tr->ll_Is = 0;

	  return tr;
	}


int ITree_addInterval(ITree tr, BElem l, BElem r, void *pt)
	{ ITreeInterval *I;
	  if(tr==NULL) {
		fprintf(stderr,"Error in ITree_addInterval: the tree tr is not initialized\n");
		exit(1);
		}
	  
	  if(l>r){
		fprintf(stderr,"Error in ITree_addInterval: left end point (%g) is larger than the right end point (%g)\n", (double) l, (double) r);
		exit(1);
		}

	  /*not enough space, allocate more*/
	  if(tr->ll_Is+10 > tr->size_Is){
		tr->size_Is += 1000;
		tr->Is = (ITreeInterval *) realloc(tr->Is,sizeof(ITreeInterval)*(tr->size_Is));
		if(tr->Is==NULL){
			fprintf(stderr,"Error in Tree_addInterval: memory allocation failed\n");
			exit(1);	
			}
		}

	   tr->Is[tr->ll_Is].l = l;
	   tr->Is[tr->ll_Is].r = r;
	   tr->Is[tr->ll_Is].max = l-1;/*When adding this interval to the tree, I don't know what is the max is yet, so set it as l-1*/
	   tr->Is[tr->ll_Is].pt = pt;
	   I = &(tr->Is[tr->ll_Is]);
	   rbtree_insert(tr->rbtr, (void *)I, (void *)I, cmpInterval_left); /*add a new node to the rbtree*/
	   tr->ll_Is++;	   

	   return 0;
	}

int cmpInterval_left(void *a, void *b)
	{ ITreeInterval *tmp1, *tmp2;
	  tmp1 = (ITreeInterval *) a;
	  tmp2 = (ITreeInterval *) b;

	  if(tmp1->l < tmp2->l) return -1;
	  else if(tmp1->l > tmp2->l) return 1;
	  else { if(tmp1->r < tmp2-> r) return -1;
		 else if(tmp1->r > tmp2->r) return 1;
		 else { if(tmp1->pt < tmp2->pt) return -1;
			else if(tmp1->pt > tmp2->pt) return 1;
			else return 0;
			}
		}
	}



void ITreeDestroy(ITree tr)
	{ ITreeInterval *I;
	  int i;
	  
	  if(tr==NULL) return;
	  
	  if(tr->ll_Is>0){
		for(i=0;i<tr->ll_Is;i++){
			I = &(tr->Is[i]);
			rbtree_delete(tr->rbtr, (void *)I, cmpInterval_left);
			}
		if(tr->rbtr->root!=NULL) {fprintf(stderr,"Warining in ITreeDestroy: tree is still not empty.\n");}
		}
	  //printf("rbtree freed\n");getchar();
	  if(tr->Is!=NULL){free(tr->Is); tr->Is=NULL; tr->ll_Is=0; tr->size_Is = 0;}
	  if(tr->rbtr!=NULL){free(tr->rbtr); tr->rbtr=NULL;}

	  return;
	}


static BElem find_max_subrbtree(rbtree_node node)
	{ BElem max_left, max_right;
	  ITreeInterval *I;
	  
	  if(node==NULL) {fprintf(stderr,"Error in 'find_max_subrbtree': NULL node\n");exit(1);}
	  
	  I = (ITreeInterval *) node->value;
	  if(I->max<I->l){ /*I->max is unknown*/
		if(node->left==NULL&&node->right==NULL){
			I->max = I->r;
			}else{
				if(node->left!=NULL){
				max_left = find_max_subrbtree(node->left);
				}
				if(node->right!=NULL){
				max_right = find_max_subrbtree(node->right);
				}
			 I->max = max_right > max_left ? max_right : max_left;
			}
		}

	 return I->max;
	}

void ITreeConstruct(ITree tr)
	{ if(tr==NULL||tr->rbtr==NULL||tr->rbtr->root==NULL) return;
	  find_max_subrbtree(tr->rbtr->root);
	  return;
	}


static void ITreeQuery_node(rbtree_node node, BElem l, BElem r, ITreeIntervalSet *Iset)
	{ ITreeInterval *I;
	  if(node!=NULL){
		I = (ITreeInterval *) node->value;
		if(l<=I->r&&I->l<=r){ /*overlap with this interval*/
			if(Iset->ll+2 > Iset->size){ /*not enough space*/
				Iset->size += 100;
				Iset->Is = (ITreeInterval **) realloc(Iset->Is, sizeof(ITreeInterval *)*Iset->size);
				if(Iset->Is==NULL){
					fprintf(stderr,"Error in ITreeQuery_node: memory allocation failed\n");
					exit(1);
					}
				}
			Iset->Is[Iset->ll] = I;
			Iset->ll++;
			}

		if(r>=I->l){ /*it is possible that [l,r] could overlap with intervals on the right subtree of node*/
			ITreeQuery_node(node->right,l,r, Iset);
			}
		if(l<=I->max){ /*it is possible that [l,r] could overlap with intervals on the left subtree of node*/
			ITreeQuery_node(node->left,l,r, Iset);
			}
		}
	 return;
	}

void ITreeQuery(ITree tr, BElem l, BElem r, ITreeIntervalSet *Iset)
	{ if(l > r) {
		fprintf(stderr,"Error in TreeQuery: invalid interval [%g,%g]",(double) l,(double)r);
		exit(1);
		}
	  if(Iset->ll!=0) {Iset->ll=0;}
	  if(Iset->size<0){fprintf(stderr,"Error in TreeQuery: invalid size (%g) of IntervalSet; Iset uninitialized?",(double) Iset->size); exit(1);}

	  if(tr!=NULL&&tr->rbtr!=NULL&&tr->rbtr->root!=NULL&&tr->ll_Is>0){
	  	ITreeQuery_node(tr->rbtr->root,l,r,Iset);
		}
	  return; 
	}


ITreeIntervalSet * ITreeIntervalSet_create()
	{ ITreeIntervalSet *Iset = malloc(sizeof(ITreeIntervalSet));
	  if(Iset==NULL) {
		fprintf(stderr,"Error in ITreeIntervalSet_create: memory allocation failed\n");
		exit(1);
		}
	  Iset->Is = NULL;
	  Iset->size = 0;
	  Iset->ll = 0;
	  return Iset;
	}


void ITreeIntervalSet_destroy(ITreeIntervalSet * Iset)
	{ if(Iset!=NULL){
		if(Iset->Is!=NULL&&Iset->size>0) free(Iset->Is);
		Iset->Is = NULL;
		Iset->ll=0;
		Iset->size = 0; 
		free(Iset);
		}
	  return;
	}
