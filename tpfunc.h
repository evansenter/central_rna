#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double tcalculateZBBLpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateZBBRpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateZBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR);
double tcalculateZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1);
double tcalculateZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM);
double tcalculateZpfns(int *seq_int, int i, int j, double **ZB, double **Z);
double GetTPartFunc(int * trans_seq, int n); // seq and ss  
