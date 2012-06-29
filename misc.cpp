#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"misc.h"
#include"ProgConst.h"

extern int DEBUG;

int d2[13] = {0,13,0,14,0,0,0,15,0,0,0,0,16};

double R = 0.001987;
double T = 310.00;

/* 
  d1 breakdown

  00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 
  AU UA GC CG GU UG AC AG CA CU GA UC GG AA CC UU

  00 01 02 03 04 05 06 07
  0A 0U 0C 0G A0 U0 C0 G0
*/
                  // 0 1  2 3  4 5 6  7 8 9 10 11  12 13 14 15 16  17 18 19  20 21  22  23 24  25  26  27  28
int d1[29] = {1,4,-1,0,-1,2,2, 9,6,0,10, 3, -1, 7, 7, 4, 8, 11, 3, 6, -1, 1, -1,  5, 5, 12, 13, 14, 15};

//cannonical basepairs indexed by d1[getIndex(i,j)]
int cbp[29] = {1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//guau il penalty indexed by d1
int guau[29] = {1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

// function takes sequence and converts to integer array through dictionary G->1, A->3, C->7, U->12
// 1-indexed.
int * TransformSeq(char * sequence)
{
  int * integer_seq;
  int i;
  char * start;
  integer_seq = (int *)malloc(strlen(sequence)*sizeof(int));
  i = 0; 
  start = sequence;
  while(*start != '\0'){
 
    if (*start == '0'){
      integer_seq[i] = strlen(sequence)-1;
    }
    else if(*start == 'G'){
      integer_seq[i] = 1;
    }
    else if(*start == 'A'){
      integer_seq[i] = 3;
    }
    else if(*start == 'C'){
      integer_seq[i] = 7;
    }
    else if(*start == 'U'){
      integer_seq[i] = 12;
    }
   
    start = start+ sizeof(char);
    i++;
  }
  return integer_seq;
}

// reason we have d2 is to account for AA,GG,CC,UU, +12 because integer value for U is 12 
int GetIndex(int i, int j)
{
  int index;
  index = i-j+12+((i==j)*d2[i]);
  return index;
}

void PrintStr(char *ptr)
{
  char *start;
  start = ptr;
  start = start+ sizeof(char);
  while(*start != '\0'){
    printf("%c",*start);
    start = start + sizeof(char);
  }
  printf("\n");
}

void PrintInt(int *ptr, int len)
{
  int i;
  for(i = 0; i < len+1; i++){
    printf("%d", ptr[i]);
  }
  printf("\n");
}

char GetNucl(int trans){

  char nucl;

  if (trans == 1){
    nucl = 'G';
  }

  else if(trans == 3){
    nucl = 'A';
  }

  else if(trans == 7){
    nucl = 'C';
  }

  else if(trans == 12){
    nucl = 'U';
  }

  return nucl;
}

int CountChar(char * str, char c){

  int count;
  count = 0;

  while(*str != '\0'){
    if(*str == c){
      count++;
    }
    str++;
  }
  
  return count;
}

int **Allocate2Dmatrix(int a, int b){
  int i = 0;
  int j = 0;
  int **Matrix;
  

  Matrix = (int **) calloc(a,sizeof(int *));
  if(Matrix == NULL){
    printf("out of memory\n");
    exit(1);
  }
  for (i = 0; i <a; i++){
    Matrix[i] = (int *) calloc(b,sizeof(int));
    if(Matrix[i] == NULL){
      printf("out of memory\n");
      exit(1);
    }
  }

  return Matrix;
}

double **Allocate2DMatrixDouble(int a, int b){
  int i = 0;
  int j = 0;
  double **Matrix;

  Matrix = (double **) calloc(a,sizeof(double *));
  if(Matrix == NULL){
    printf("out of memory\n");
    exit(1);
  }
  for (i = 0; i <a; i++){
    Matrix[i] = (double *) calloc(b,sizeof(double));
    if(Matrix[i] == NULL){
      printf("out of memory\n");
      exit(1);
    }
  }
  return Matrix;
}

// 1-indexed.
int **BPBetweenBool(int** Matrix, int * bp, int len){

  int pos,t,x,i,j,k,l,i_index,k_index;

  for(i_index = 0; i_index < len; i_index++){

    if(bp[i_index] != -1){
   
      i = i_index;
      j = bp[i];

      for(k_index=i_index+1;k_index < j; k_index++){
        if (bp[k_index] != -1){
          k = k_index;
          l = bp[k];
          break;
        }
      }
  
      if(k > i && j > l){
        Matrix[i][j] = 1;
      }
    
      for(t=l+1;t<j;t++){
        x = bp[t];
        if(x!=-1){
          Matrix[l][j] = 1;
          break;
        }
      }
    }
  }
  
  if(DEBUG){
    for(i = 0; i < len; i++){
      for(j = 0; j < len; j++){
        printf("%d ",Matrix[i][j]);
      }
      printf("\n");
    }
  }  
  
  return Matrix;
}

// 1-indexed.
int * GetBPList(char * str, int len){

  int i = 0;
  int pointer = 0;
  int * left;
  int * bp;
  int leftbp;

  left = (int *)malloc(len*sizeof(int));
  bp = (int *)malloc(len*sizeof(int));

  for (i = 0; i < len; i++){
    bp[i] = -1;
  }
  
  for (i = 0; i < len; i++){

    if (str[i] == '('){
      left[pointer] = i;
      pointer++;
    }
  
    else if(str[i] == ')'){
      if (pointer == 0){
        printf("Unbalanced Secondary Structure\n");
        exit(1);
      }
      leftbp = left[pointer-1];
      bp[leftbp] = i;
      bp[i] = leftbp;
      pointer--;
    }
  }

  if(DEBUG){
   
    for(i = 0; i < len; i++){
      printf("%d %d\n",i,bp[i]);
    }
  
  }

  free(left);

  return bp;
}

int linearize2D(int row, int column){

  int index;
  index = row*16+column;

  return index;
}
