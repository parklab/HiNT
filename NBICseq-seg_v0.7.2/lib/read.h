/* This is the read data program written by Ruibin Xi Oct 12, 2009.
*/

#ifndef READ_HX_XI
#define READ_HX_XI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>


#define LR_BUFFER 1000
/*The minimum memory allocated to ll_buffer*/

char * lr_read(FILE *); /*read a line; if error, return NULL*/

double * read_table(FILE *, int *nrow, int *ncol, int lines,int skip);
/*read a file into a double array (a matrix).
  nrow will be assigned as the number of rows of the matrix;
  ncol will be assigned as the number of cols of the matrix;
  lines tells maximum number of lines that need to be read, -1 means that read until the end of the file;
  If the end of file is reached before lines number of lines, then the returned matrix will have less than lines number of rows.
  ignore any row that has different columns from the first row.
  skip is the number of lines that need to be skipped
*/

#endif

