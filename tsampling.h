#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double tcalculateSZBBLpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSZBBRpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSZBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateSZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double tcalculateSZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double tcalculateSZpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double tsample(int * seq_int, int n, int samples); // seq and ss 
int tsampleZBBL(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int);
int tsampleZBBR(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int);
int tsampleZBB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int);
int tsampleZB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int);
int tsampleZM1(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int);
int tsampleZM(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int);
int tsampleZ(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, char *str, long *seed, int *seq_int); 
