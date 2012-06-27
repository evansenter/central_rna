#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double tcalculateSP1ZBBLpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP1ZBBRpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP1ZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP1ZBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSP1ZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double tcalculateSP1ZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double tcalculateSP1Zpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double tcalculateProbs1(int * seq_int, int n, int samples); // seq and ss  
int tbacktrackSP1ZBBL(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
int tbacktrackSP1ZBBR(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
int tbacktrackSP1ZBB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
int tbacktrackSP1ZB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
int tbacktrackSP1ZM1(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
int tbacktrackSP1ZM(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
int tbacktrackSP1Z(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int);
