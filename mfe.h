#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double calculateEB(int *seq_int, int i, int j, double **EB, double **EM, double **EM1, int n);
double calculateEM1(int *seq_int, int i, int j, double **EB, double **EM1);
double calculateEM(int *seq_int, int i, int j, double **EM1, double **EM);
double calculateE(int *seq_int, int i, int j, double **EB, double **E);
void  GetMFE(int * trans_seq, int n, char * mfe_seq); // seq and ss 
int backtrackE(int i, int j, double **E, double **EB, double **EM1, double **EM, char *ss, int n);
int backtrackEB(int i, int j, double **EB, double **EM1, double **EM, char *ss, int n);
int backtrackEM(int i, int j, double **EB, double **EM1, double **EM, char *ss, int n);
int backtrackEM1(int i, int j, double **EB, double **EM1, double **EM, char *ss, int n);
