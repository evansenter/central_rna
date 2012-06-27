#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double calculateZBpfns(int *seq_int, int i, int j, double **ZB, double **ZM, double **ZM1, int n);
double calculateZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double calculateZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double calculateZpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double GetPartFunc(int * trans_seq, int n); // seq and ss  
