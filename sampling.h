#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double calculateSZBpfns(int *seq_int, int i, int j, double **ZB, double **ZM, double **ZM1);
double calculateSZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double calculateSZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double calculateSZpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double sample(int * seq_int, int n, int samples); // seq and ss 
int sampleZB(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, char *str, long *seed, int *seq_int);
int sampleZM1(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, char *str, long *seed, int *seq_int);
int sampleZM(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, char *str, long *seed, int *seq_int);
int sampleZ(int i, int j, double **Z, double **ZB,  double **ZM, double **ZM1, char *str, long *seed, int *seq_int); 
