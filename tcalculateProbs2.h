#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double tcalculateSP2ZBBLpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP2ZBBRpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP2ZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP2ZBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP2ZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double tcalculateSP2ZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double tcalculateSP2Zpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double tcalculateProbs2(int * trans_seq, int n, int samples); // seq and ss  
int tbacktrackSP2ZBBL(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
int tbacktrackSP2ZBBR(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
int tbacktrackSP2ZBB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
int tbacktrackSP2ZB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
int tbacktrackSP2ZM1(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
int tbacktrackSP2ZM(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
int tbacktrackSP2Z(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, int samples, int *seq_int);
