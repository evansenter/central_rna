#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double calculateSP1ZBpfns(int *seq_int, int i, int j, double **ZB, double **ZM, double **ZM1);
double calculateSP1ZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double calculateSP1ZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double calculateSP1Zpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double calculateProbs1(int * seq_int, int n, int samples); // seq and ss  
int backtrackSP1ZB(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int);
int backtrackSP1ZM1(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int);
int backtrackSP1ZM(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int);
int backtrackSP1Z(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int);
