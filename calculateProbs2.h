#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double calculateSP2ZBpfns(int *seq_int, int i, int j, double **ZB, double **ZM, double **ZM1);
double calculateSP2ZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double calculateSP2ZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double calculateSP2Zpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double calculateProbs2(int * trans_seq, int n, int samples); // seq and ss  
int backtrackSP2ZB(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, int samples, int *seq_int);
int backtrackSP2ZM1(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, int samples, int *seq_int);
int backtrackSP2ZM(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, int samples, int *seq_int);
int backtrackSP2Z(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, int samples, int *seq_int);
