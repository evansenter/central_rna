#include "ProgConst.h"

/* declaration of functions for initializaing, transformation */

extern int d2[13];
extern int d1[29];
extern int cbp[29];
extern int guau[29];
extern double R;
extern double T;

int * TransformSeq(char sequence[MAXSIZE]); // converts sequence to list of int through dictionary G->1, A->3, C->7, U->12
int GetIndex(int i, int j); // returns value in  d1 through index from formula (i-j)+12+(i==j)[d2[i]]
void PrintStr(char *ptr); // prints string given ptr, end character must be \0
void PrintInt(int *ptr, int len);
int CountChar(char * str, char c);// counts the number of char c in str
int **Allocate2Dmatrix(int a, int b);// Allocate 2D Matrix of axb
double **Allocate2DMatrixDouble(int a, int b);// Allocate 2D double Matrix of axb
int **BPBetweenBool(int** Matrix, int * bp, int len);//fills in matrix such that matrix[i][j] returns 1 if BP between i,j, otherwise 0
int * GetBPList(char * str,int len);
char GetNucl(int trans);
int linearize2D(int row, int column);
