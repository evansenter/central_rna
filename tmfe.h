#include "misc.h"
//#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */
double tcalculateEBBL(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n);
double tcalculateEBBR(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n);
double tcalculateEBB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n);
double tcalculateEB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n);
double tcalculateEM1(int *seq_int, int i, int j, double **EB, double **EM1);
double tcalculateEM(int *seq_int, int i, int j, double **EM1, double **EM);
double tcalculateE(int *seq_int, int i, int j, double **EB, double **E);
void GetTMFE(int * trans_seq, int n, char *mfe_seq); // seq and ss 
int tbacktrackE(int i, int j, double **E, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n);
int tbacktrackEB(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n);
int tbacktrackEBBL(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n);
int tbacktrackEBBR(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n);
int tbacktrackEBB(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n);
int tbacktrackEM(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n);
int tbacktrackEM1(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR,char *ss, int n);
