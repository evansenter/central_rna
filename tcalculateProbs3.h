#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double calculateSP3ZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1);
double calculateSP3ZBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1);
double calculateSP3ZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double calculateSP3ZMpfns(int *seq_int, int i, int j, double **ZM, double **ZM1);
double calculateSP3Zpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double calculateProbs3(int * trans_seq, int n); // seq and ss  
