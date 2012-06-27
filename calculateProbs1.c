#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"energy_par.h"
#include"energy_func.h"
//#include"tpfunc.h"
#include <math.h>
#include "random.h"
#include "tcalculateProbs1.h"
 
//minimum unpaired nts in a loop
#define SIGMA 3
#define MAX_PREC 1000000
/*double calculateSP1ZBBLpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR){
  if(j-i<=SIGMA+3)
    return 0;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  int ij_index2 = d1[GetIndex(seq_int[i+2],seq_int[j-1])];
  if((cbp[ij_index]==0)||(cbp[ij_index2]==0))
    return 0;

  double RT = R*T;

  ij_index = d1[GetIndex(seq_int[i+3],seq_int[j-2])];
  int kj_index, ik_index,lr_index;
  double QSns = cbp[ij_index]*ZBB[j-1][i+2];
  double gte = GetTripletEnergy(i,j,i+2,j-1,i+3,j-2,seq_int);
  double egte = exp(-1*(gte/RT));
  double QSpf = cbp[ij_index]*ZBB[i+2][j-1]*egte;

  double QHns = 1;
  double QHpf = exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetHairpinEnergy(i+2,j-1,seq_int))/RT));

  kj_index = d1[GetIndex(seq_int[i+4],seq_int[j-2])];
  double QLBns = cbp[kj_index]*ZBBL[j-1][i+2];
  double QLBpf = cbp[kj_index]*ZBBL[i+2][j-1]*exp(-1*((GetTripletEnergy(i,j,i+2,j-1,i+4,j-2,seq_int)+GetLeftBulgeEnergy(i+2,j-1,i+4,seq_int))/RT));

  int k,l,r;
  for(k = i+5; (k <= j-SIGMA-3); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-2])];
    QLBns += cbp[kj_index]*ZB[j-2][k];
    QLBpf += cbp[kj_index]*ZB[k][j-2]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetLeftBulgeEnergy(i+2,j-1,k,seq_int))/RT));
  }

  ik_index = d1[GetIndex(seq_int[i+3],seq_int[j-3])];
  double QRBns = cbp[ik_index]*ZBBR[j-1][i+2];
  double QRBpf = cbp[ik_index]*ZBBR[i+2][j-1]*exp(-1*((GetTripletEnergy(i,j,i+2,j-1,i+3,j-3,seq_int)+GetRightBulgeEnergy(i+2,j-1,j-3,seq_int))/RT));

  for(k = j-4; (k>=i+SIGMA+4); k--){
    ik_index = d1[GetIndex(seq_int[i+3],seq_int[k])];
    QRBns += cbp[ik_index]*ZB[k][i+3];
    QRBpf += cbp[ik_index]*ZB[i+3][k]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetRightBulgeEnergy(i+2,j-1,k,seq_int))/RT));
  }

  double QIns = 0;
  double QIpf = 0;

for(l = i+4; (l <= j-SIGMA-4)&&(l-i-1<=30); l++){
    for(r = j-3; (r>=l+SIGMA+1)&&(j-r+l-i-1<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      QIns += cbp[lr_index]*ZB[r][l];
      QIpf += cbp[lr_index]*ZB[l][r]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetILEnergy(i+2,j-1,l,r,seq_int))/RT));
    }
  }
  double QMns = 0;//ZM1[j-2][i+2];
  double QMpf = 0;//ZM1[i+2][j-2];
  for(k = i+SIGMA+5; k <= j-SIGMA-3; k++){
    QMns += ZM[k-1][i+3]*ZM1[j-2][k];
    QMpf += ZM[i+3][k-1]*ZM1[k][j-2];
  }
  QMpf = QMpf*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+MultiloopA+MultiloopB)/RT));

  ZBBL[i][j] = QSpf + QHpf + QLBpf + QRBpf + QIpf + QMpf;

  return (QSns+QHns+QLBns+QRBns+QIns+QMns);
}


int backtrackSP1ZBBL(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int){
  if(j-i<=SIGMA+3)
    return 0;
  //printf("ZBB[%d][%d] = %f\n",i,j,ZBB[i][j]);
  p[i+2][j-1]++;
  p[j-1][i+2]++;
 
  double RT = R*T;
  int cumulative = 0;
  int rd = floor(ran2(seed)*MAX_PREC);
  int k,l,r;
  int ij_index = d1[GetIndex(seq_int[i+3],seq_int[j-2])];
  int kj_index, ik_index, lr_index;

  cumulative = ceil((cbp[ij_index]*ZBB[i+2][j-1]*(exp(-1*(GetTripletEnergy(i,j,i+2,j-1,i+3,j-2,seq_int)/RT))))/ZBBL[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBB(i+2,j-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  cumulative += ceil(( exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetHairpinEnergy(i+2,j-1,seq_int))/RT)))/ZBBL[i][j]*MAX_PREC);
  if(rd<cumulative)
    return 0;

  kj_index = d1[GetIndex(seq_int[i+4],seq_int[j-2])];
  cumulative += ceil((cbp[kj_index]*ZBBL[i+2][j-1]*exp(-1*((GetTripletEnergy(i,j,i+2,j-1,i+4,j-2,seq_int)+GetLeftBulgeEnergy(i+2,j-1,i+4,seq_int))/RT)))/ZBBL[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBBL(i+2,j-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  for(k = i+5; (k <= j-SIGMA-3); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-2])];
    cumulative += ceil((cbp[kj_index]*ZB[k][j-2]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetLeftBulgeEnergy(i+2,j-1,k,seq_int))/RT)))/ZBBL[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(k,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
  }

  ik_index = d1[GetIndex(seq_int[i+3],seq_int[j-3])];
  cumulative += ceil((cbp[ik_index]*ZBBR[i+2][j-1]*exp(-1*((GetTripletEnergy(i,j,i+2,j-1,i+3,j-3,seq_int)+GetRightBulgeEnergy(i+2,j-1,j-3,seq_int))/RT)))/ZBBL[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBBR(i+2,j-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  for(k = j-4; (k>=i+SIGMA+4); k--){
    ik_index = d1[GetIndex(seq_int[i+3],seq_int[k])];
    cumulative += ceil((cbp[ik_index]*ZB[i+3][k]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetRightBulgeEnergy(i+2,j-1,k,seq_int))/RT)))/ZBBL[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(i+3,k,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
  }

  for(l = i+4; (l <= j-SIGMA-4)&&(l-i-1<=30); l++){
    for(r = j -3; (r>=l+SIGMA+1)&&(j-r+l-i-1<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      cumulative += ceil((cbp[lr_index]*ZB[l][r]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+GetILEnergy(i+2,j-1,l,r,seq_int))/RT)))/ZBBL[i][j]*MAX_PREC);
      if(rd<cumulative)
       return backtrackSP1ZB(l,r,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
    }
  }

  for(k = i+SIGMA+5; k <= j-SIGMA-3; k++){
    cumulative += ceil((ZM[i+3][k-1]*ZM1[k][j-2]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int)+MultiloopA+MultiloopB)/RT)))/ZBBL[i][j]*MAX_PREC);
    if(rd<cumulative){
       backtrackSP1ZM(i+3,k-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
       backtrackSP1ZM1(k,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
       return 0;
    }
  }
  //printf("Nothing? rd = %f\n",rd);
}                      

double calculateSP1ZBBRpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR){
  if(j-i<=SIGMA+3)
    return 0;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  int ij_index2 = d1[GetIndex(seq_int[i+1],seq_int[j-2])];
  if((cbp[ij_index]==0)||(cbp[ij_index2]==0))
    return 0;

  double RT = R*T;

  ij_index = d1[GetIndex(seq_int[i+2],seq_int[j-3])];
  int kj_index, ik_index,lr_index;
  double QSns = cbp[ij_index]*ZBB[j-2][i+1];
  double gte = GetTripletEnergy(i,j,i+1,j-2,i+2,j-3,seq_int);
  double egte = exp(-1*(gte/RT));
  double QSpf = cbp[ij_index]*ZBB[i+1][j-2]*egte;

  double QHns = 1;
  double QHpf = exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetHairpinEnergy(i+1,j-2,seq_int))/RT));

  kj_index = d1[GetIndex(seq_int[i+3],seq_int[j-3])];
  double QLBns = cbp[kj_index]*ZBBL[j-2][i+1];
  double QLBpf = cbp[kj_index]*ZBBL[i+1][j-2]*exp(-1*((GetTripletEnergy(i,j,i+1,j-2,i+3,j-3,seq_int)+GetLeftBulgeEnergy(i+1,j-2,i+3,seq_int))/RT));

  int k,l,r;
  for(k = i+4; (k <= j-SIGMA-4); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-3])];
    QLBns += cbp[kj_index]*ZB[j-3][k];
    QLBpf += cbp[kj_index]*ZB[k][j-3]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetLeftBulgeEnergy(i+1,j-2,k,seq_int))/RT));
  }

  ik_index = d1[GetIndex(seq_int[i+2],seq_int[j-4])];
  double QRBns = cbp[ik_index]*ZBBR[j-2][i+1];
  double QRBpf = cbp[ik_index]*ZBBR[i+1][j-2]*exp(-1*((GetTripletEnergy(i,j,i+1,j-2,i+2,j-4,seq_int)+GetRightBulgeEnergy(i+1,j-2,j-4,seq_int))/RT));

  for(k = j-5; (k>=i+SIGMA+3); k--){
    ik_index = d1[GetIndex(seq_int[i+2],seq_int[k])];
    QRBns += cbp[ik_index]*ZB[k][i+2];
    QRBpf += cbp[ik_index]*ZB[i+2][k]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetRightBulgeEnergy(i+1,j-2,k,seq_int))/RT));
  }

  double QIns = 0;
  double QIpf = 0;

  for(l = i+3; (l <= j-SIGMA-5)&&(l-i-2<=30); l++){
    for(r = j-4; (r>=l+SIGMA+1)&&(j-r+l-i-3<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      QIns += cbp[lr_index]*ZB[r][l];
      QIpf += cbp[lr_index]*ZB[l][r]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetILEnergy(i+1,j-2,l,r,seq_int))/RT));
    }
  }
  double QMns = 0;//ZM1[j-2][i+2];
  double QMpf = 0;//ZM1[i+2][j-2];
  for(k = i+SIGMA+4; k <= j-SIGMA-4; k++){
    QMns += ZM[k-1][i+2]*ZM1[j-3][k];
    QMpf += ZM[i+2][k-1]*ZM1[k][j-3];
  }
  QMpf = QMpf*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+MultiloopA+MultiloopB)/RT));

  ZBBR[i][j] = QSpf + QHpf + QLBpf + QRBpf + QIpf + QMpf;

  return (QSns+QHns+QLBns+QRBns+QIns+QMns);
}

int backtrackSP1ZBBR(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int){
  if(j-i<=SIGMA+3)
    return 0;
  //printf("ZBB[%d][%d] = %f\n",i,j,ZBB[i][j]);
  p[i+1][j-2]++;
  p[j-2][i+1]++;
 
  if(j-i<=SIGMA+4)
    return 0;
  double RT = R*T;
  int cumulative = 0;
  int rd = floor(ran2(seed)*MAX_PREC);
  int k,l,r;
  int ij_index = d1[GetIndex(seq_int[i+2],seq_int[j-3])];
  int kj_index, ik_index, lr_index;

  cumulative = ceil((cbp[ij_index]*ZBB[i+1][j-2]*(exp(-1*(GetTripletEnergy(i,j,i+1,j-2,i+2,j-2,seq_int)/RT))))/ZBBR[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBB(i+1,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  cumulative += ceil(( exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetHairpinEnergy(i+1,j-2,seq_int))/RT)))/ZBBR[i][j]*MAX_PREC);
  if(rd<cumulative)
    return 0;

  kj_index = d1[GetIndex(seq_int[i+3],seq_int[j-3])];
  cumulative += ceil((cbp[kj_index]*ZBBL[i+1][j-2]*exp(-1*((GetTripletEnergy(i,j,i+1,j-2,i+3,j-3,seq_int)+GetLeftBulgeEnergy(i+1,j-2,i+3,seq_int))/RT)))/ZBBR[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBBL(i+1,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  for(k = i+4; (k <= j-SIGMA-4); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-3])];
    cumulative += ceil((cbp[kj_index]*ZB[k][j-3]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetLeftBulgeEnergy(i+1,j-2,k,seq_int))/RT)))/ZBBR[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(k,j-3,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
  }

  ik_index = d1[GetIndex(seq_int[i+2],seq_int[j-4])];
  cumulative += ceil((cbp[ik_index]*ZBBR[i+1][j-2]*exp(-1*((GetTripletEnergy(i,j,i+1,j-2,i+2,j-4,seq_int)+GetRightBulgeEnergy(i+1,j-2,j-4,seq_int))/RT)))/ZBBR[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBBR(i+1,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  for(k = j-5; (k>=i+SIGMA+3); k--){
    ik_index = d1[GetIndex(seq_int[i+2],seq_int[k])];
    cumulative += ceil((cbp[ik_index]*ZB[i+2][k]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetRightBulgeEnergy(i+1,j-2,k,seq_int))/RT)))/ZBBR[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(i+2,k,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
  }

  for(l = i+3; (l <= j-SIGMA-5)&&(l-i-2<=30); l++){
    for(r = j -4; (r>=l+SIGMA+1)&&(j-r+l-i-3<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      cumulative += ceil((cbp[lr_index]*ZB[l][r]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+GetILEnergy(i+1,j-2,l,r,seq_int))/RT)))/ZBBR[i][j]*MAX_PREC);
      if(rd<cumulative)
       return backtrackSP1ZB(l,r,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
    }
  }

  for(k = i+SIGMA+4; k <= j-SIGMA-4; k++){
    cumulative += ceil((ZM[i+2][k-1]*ZM1[k][j-3]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int)+MultiloopA+MultiloopB)/RT)))/ZBBR[i][j]*MAX_PREC);
    if(rd<cumulative){
       backtrackSP1ZM(i+2,k-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
       backtrackSP1ZM1(k,j-3,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
       return 0;
    }
  }
  //printf("Nothing? rd = %f\n",rd);
}                      

double calculateSP1ZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR){
  if(j-i<=SIGMA+2)
    return 0;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  int ij_index2 = d1[GetIndex(seq_int[i+1],seq_int[j-1])];
  if((cbp[ij_index]==0)||(cbp[ij_index2]==0))
    return 0;

  double RT = R*T;
 
  ij_index = d1[GetIndex(seq_int[i+2],seq_int[j-2])];
  int kj_index, ik_index,lr_index;
  //double QSns = cbp[ij_index]*ZBB[j-1][i+1];
  double gte = GetTripletEnergy(i,j,i+1,j-1,i+2,j-2,seq_int);
  double egte = exp(-1*(gte/RT));
  double QSpf = cbp[ij_index]*ZBB[i+1][j-1]*egte;

  //double QHns = 1;
  double QHpf = exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetHairpinEnergy(i+1,j-1,seq_int))/RT));

  kj_index = d1[GetIndex(seq_int[i+3],seq_int[j-2])];
  //double QLBns = 0;
  double QLBpf = cbp[kj_index]*ZBBL[i+1][j-1]*exp(-1*((GetTripletEnergy(i,j,i+1,j-1,i+3,j-2,seq_int)+GetLeftBulgeEnergy(i+1,j-1,i+3,seq_int))/RT));
  int k,l,r;
  for(k = i+4; (k <= j-SIGMA-3); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-2])];
    //QLBns += cbp[kj_index]*ZB[j-2][k];
    QLBpf += cbp[kj_index]*ZB[k][j-2]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetLeftBulgeEnergy(i+1,j-1,k,seq_int))/RT));
  }  

  ik_index = d1[GetIndex(seq_int[i+2],seq_int[j-3])];
  //double QRBns = 0;
  double QRBpf = cbp[ik_index]*ZBBR[i+1][j-1]*exp(-1*((GetTripletEnergy(i,j,i+1,j-1,i+2,j-3,seq_int)+GetRightBulgeEnergy(i+1,j-1,j-3,seq_int))/RT));
  for(k = j-4; (k>=i+SIGMA+3); k--){
    ik_index = d1[GetIndex(seq_int[i+2],seq_int[k])];
    //QRBns += cbp[ik_index]*ZB[k][i+2];
    QRBpf += cbp[ik_index]*ZB[i+2][k]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetRightBulgeEnergy(i+1,j-1,k,seq_int))/RT));
  }

  //double QIns = 0;
  double QIpf = 0;

  for(l = i+3; (l <= j-SIGMA-4)&&(l-i-2<=30); l++){
    for(r = j-3; (r>=l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      //QIns += cbp[lr_index]*ZB[r][l];
      QIpf += cbp[lr_index]*ZB[l][r]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetILEnergy(i+1,j-1,l,r,seq_int))/RT));
    }
  }
  //double QMns = 0;//ZM1[j-2][i+2];
  double QMpf = 0;//ZM1[i+2][j-2];
  for(k = i+SIGMA+4; k <= j-SIGMA-3; k++){
    //QMns += ZM[k-1][i+2]*ZM1[j-2][k];
    QMpf += ZM[i+2][k-1]*ZM1[k][j-2];
  }
  QMpf = QMpf*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+MultiloopA+MultiloopB)/RT));
 
  ZBB[i][j] = QSpf + QHpf + QLBpf + QRBpf + QIpf + QMpf;

  return 0;// (QSns+QHns+QLBns+QRBns+QIns+QMns);
}

int backtrackSP1ZBB(int i, int j, double **Z, double **ZB, double **ZBB, double **ZM, double **ZM1, double **ZBBL, double **ZBBR, double **p, long *seed, int *seq_int){
  if(j-i<=SIGMA+2)
    return 0;
  //printf("ZBB[%d][%d] = %f\n",i,j,ZBB[i][j]);
  p[i+1][j-1]++;
  p[j-1][i+1]++;
 
  if(j-i<=SIGMA+4)
    return 0;
  double RT = R*T;
  int cumulative = 0;
  int rd = floor(ran2(seed)*MAX_PREC);
  int k,l,r;
  int ij_index = d1[GetIndex(seq_int[i+2],seq_int[j-2])];
  int kj_index, ik_index, lr_index;

  cumulative = ceil((cbp[ij_index]*ZBB[i+1][j-1]*(exp(-1*(GetTripletEnergy(i,j,i+1,j-1,i+2,j-2,seq_int)/RT))))/ZBB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBB(i+1,j-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  cumulative += ceil(( exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetHairpinEnergy(i+1,j-1,seq_int))/RT)))/ZBB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return 0;

  kj_index = d1[GetIndex(seq_int[i+3],seq_int[j-2])];
  cumulative += ceil((cbp[kj_index]*ZBBL[i+1][j-1]*exp(-1*((GetTripletEnergy(i,j,i+1,j-1,i+3,j-2,seq_int)+GetLeftBulgeEnergy(i+1,j-1,i+3,seq_int))/RT)))/ZBB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBBL(i+1,j-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  for(k = i+4; (k <= j-SIGMA-3); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-2])];
    cumulative += ceil((cbp[kj_index]*ZB[k][j-2]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetLeftBulgeEnergy(i+1,j-1,k,seq_int))/RT)))/ZBB[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(k,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
  }

  ik_index = d1[GetIndex(seq_int[i+2],seq_int[j-3])];
  cumulative += ceil((cbp[ik_index]*ZBBR[i+1][j-1]*exp(-1*((GetTripletEnergy(i,j,i+1,j-1,i+2,j-3,seq_int)+GetRightBulgeEnergy(i+1,j-1,j-3,seq_int))/RT)))/ZBB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZBBR(i+1,j-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);

  for(k = j-4; (k>=i+SIGMA+3); k--){
    ik_index = d1[GetIndex(seq_int[i+2],seq_int[k])];
    cumulative += ceil((cbp[ik_index]*ZB[i+2][k]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetRightBulgeEnergy(i+1,j-1,k,seq_int))/RT)))/ZBB[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(i+2,k,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
  }

  for(l = i+3; (l <= j-SIGMA-4)&&(l-i-2<=30); l++){
    for(r = j -3; (r>=l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      cumulative += ceil((cbp[lr_index]*ZB[l][r]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetILEnergy(i+1,j-1,l,r,seq_int))/RT)))/ZBB[i][j]*MAX_PREC);
      if(rd<cumulative)
       return backtrackSP1ZB(l,r,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
    }
  }

  for(k = i+SIGMA+4; k <= j-SIGMA-3; k++){
    cumulative += ceil((ZM[i+2][k-1]*ZM1[k][j-2]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+MultiloopA+MultiloopB)/RT)))/ZBB[i][j]*MAX_PREC);
    if(rd<cumulative){
       backtrackSP1ZM(i+2,k-1,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
       backtrackSP1ZM1(k,j-2,Z,ZB,ZBB,ZM,ZM1,ZBBL,ZBBR,p,seed,seq_int);
       return 0;
    }
  }
  //printf("Nothing? rd = %f\n",rd);
}                                                                                                                                        
*/

double calculateSP1ZBpfns(int *seq_int, int i, int j, double **ZB, double **ZM, double **ZM1){
  if(j-i<=SIGMA)
    return 0;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  if(cbp[ij_index]==0)
    return 0;

  double RT = R*T;

  int k,l,r;
  ij_index = d1[GetIndex(seq_int[i+1],seq_int[j-1])];
  int kj_index, ik_index, lr_index;
  

  //double PSns = cbp[ij_index]*ZBB[j][i];
  double PSpf = cbp[ij_index]*ZB[i+1][j-1]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int))/RT));
  
  //double PHns = 1;
  double PHpf = exp(-1*(GetHairpinEnergy(i,j,seq_int)/RT));

  kj_index =  d1[GetIndex(seq_int[i+2],seq_int[j-1])];
  //double PLBns = 0;
  double PLBpf = cbp[kj_index]*ZB[i+2][j-1]*exp(-1*(GetStackEnergy(i,j,i+2,j-1,seq_int) + GetLeftBulgeEnergy(i,j,i+2,seq_int)/RT));
  for(k = i+3; (k <= j-SIGMA-2); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-1])];
    //PLBns += cbp[kj_index]*ZB[j-1][k];
    PLBpf += cbp[kj_index]*ZB[k][j-1]*exp(-1*(GetLeftBulgeEnergy(i,j,k,seq_int)/RT));
  }

  ik_index = d1[GetIndex(seq_int[i+1],seq_int[j-2])];
  //double PRBns = 0;
  double PRBpf = cbp[ik_index]*ZB[i+1][j-2]*exp(-1*(GetStackEnergy(i,j,i+1,j-2,seq_int) + GetRightBulgeEnergy(i,j,j-2,seq_int)/RT));
  for(k = j-3; (k>=i+SIGMA+2); k--){
    ik_index = d1[GetIndex(seq_int[i+1],seq_int[k])];
    //PRBns += cbp[ik_index]*ZB[k][i+1];
    PRBpf += cbp[ik_index]*ZB[i+1][k]*exp(-1*(GetRightBulgeEnergy(i,j,k,seq_int)/RT));
  }
  //double PIns = 0;
  double PIpf = 0;
  for(l = i+2; (l <= j-SIGMA-3)&&(l-i-2<=30); l++){
    for(r = j -2; (r>=l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
    //  PIns += cbp[lr_index]*ZB[r][l];
      PIpf += cbp[lr_index]*ZB[l][r]*exp(-1*(GetILEnergy(i,j,l,r,seq_int)/RT));
    }
  }
  //double PMns = 0;//ZM1[j-2][i+2];
  double PMpf = 0;//ZM1[i+2][j-2];
  for(k = i+SIGMA+3; k <= j-SIGMA-2; k++){
    //PMns += ZM[k-1][i+1]*ZM1[j-1][k];
    PMpf += ZM[i+1][k-1]*ZM1[k][j-1];
  }
  PMpf = PMpf*exp(-1*((MultiloopA+2*MultiloopB)/RT));


  ZB[i][j] = PSpf + PHpf + PLBpf + PRBpf + PIpf + PMpf;
  return 0;//(PSns+PHns+PLBns+PRBns+PIns+PMns);
}

int backtrackSP1ZB(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int){
  if(j-i<=SIGMA)
    return 0;
  //printf("ZB[%d][%d] = %f\n",i,j,ZB[i][j]);
  p[i][j]++;
  p[j][i]++;

  if(j-i<=SIGMA+2)
    return 0;

  double RT = R*T;
  int cumulative = 0;  
  int rd = floor(ran2(seed)*MAX_PREC);
  int k,l,r;
  int ij_index = d1[GetIndex(seq_int[i+1],seq_int[j-1])];
  int kj_index, ik_index, lr_index;
  
  cumulative = ceil((cbp[ij_index]*ZB[i+1][j-1]*exp(-1*(GetStackEnergy(i,j,i+1,j-1,seq_int)/RT)))/ZB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZB(i+1,j-1,Z,ZB,ZM,ZM1,p,seed,seq_int);

  cumulative += ceil((exp(-1*(GetHairpinEnergy(i,j,seq_int)/RT)))/ZB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return 0;
  
  kj_index =  d1[GetIndex(seq_int[i+2],seq_int[j-1])];
  cumulative += ceil((cbp[kj_index]*ZB[i+2][j-1]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int) + GetLeftBulgeEnergy(i,j,i+2,seq_int))/RT)))/ZB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZB(i+2,j-1,Z,ZB,ZM,ZM1,p,seed,seq_int);

  for(k = i+3; (k <= j-SIGMA-2); k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-1])];
    cumulative += ceil((cbp[kj_index]*ZB[k][j-1]*exp(-1*(GetLeftBulgeEnergy(i,j,k,seq_int)/RT)))/ZB[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(k,j-1,Z,ZB,ZM,ZM1,p,seed,seq_int);
  }

  ik_index = d1[GetIndex(seq_int[i+1],seq_int[j-2])];
  cumulative += ceil((cbp[ik_index]*ZB[i+1][j-2]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int) + GetRightBulgeEnergy(i,j,j-2,seq_int))/RT)))/ZB[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZB(i+1,j-2,Z,ZB,ZM,ZM1,p,seed,seq_int);
  
  for(k = j-3; (k>=i+SIGMA+2); k--){
    ik_index = d1[GetIndex(seq_int[i+1],seq_int[k])];
    cumulative += ceil((cbp[ik_index]*ZB[i+1][k]*exp(-1*(GetRightBulgeEnergy(i,j,k,seq_int)/RT)))/ZB[i][j]*MAX_PREC);
    if(rd<cumulative)
       return backtrackSP1ZB(i+1,k,Z,ZB,ZM,ZM1,p,seed,seq_int);
  }
 
  for(l = i+2; (l <= j-SIGMA-3)&&(l-i-2<=30); l++){
    for(r = j -2; (r>=l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      cumulative += ceil((cbp[lr_index]*ZB[l][r]*exp(-1*(GetILEnergy(i,j,l,r,seq_int)/RT)))/ZB[i][j]*MAX_PREC);
      if(rd<cumulative)
       return backtrackSP1ZB(l,r,Z,ZB,ZM,ZM1,p,seed,seq_int);
    }
  }
  
  for(k = i+SIGMA+3; k <= j-SIGMA-2; k++){
    cumulative += ceil((ZM[i+1][k-1]*ZM1[k][j-1]*exp(-1*((MultiloopA+2*MultiloopB)/RT)))/ZB[i][j]*MAX_PREC);
    if(rd<cumulative){
       backtrackSP1ZM(i+1,k-1,Z,ZB,ZM,ZM1,p,seed,seq_int);
       backtrackSP1ZM1(k,j-1,Z,ZB,ZM,ZM1,p,seed,seq_int);
       return 0;
    }
  }
  //printf("Nothing? rd = %f\n",rd);
}

double calculateSP1ZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1){
  if(j-i<=SIGMA)
    return 0;
  int k, ik_index;

  //double ZM1ns = 0;
  double ZM1pf = 0;
  double RT = R*T;

  for(k = i+SIGMA+1; k <= j; k++){
    ik_index = d1[GetIndex(seq_int[i],seq_int[k])];
    //ZM1ns += cbp[ik_index]*ZB[k][i];
    ZM1pf += cbp[ik_index]*ZB[i][k]*exp(-1*((MultiloopC*(j-k))/RT)); 
  }

  ZM1[i][j] = ZM1pf;
  return 0;//ZM1ns;
}

int backtrackSP1ZM1(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int){
   if(j-i<=SIGMA)
    return 0;
  int k, ik_index;
  double RT = R*T;
  int cumulative = 0;  
  int rd = floor(ran2(seed)*MAX_PREC);


  for(k = i+SIGMA+1; k <= j; k++){
    ik_index = d1[GetIndex(seq_int[i],seq_int[k])];
    cumulative += ceil((cbp[ik_index]*ZB[i][k]*exp(-1*((MultiloopC*(j-k))/RT)))/ZM1[i][j]*MAX_PREC);
    if(rd<cumulative)
      return backtrackSP1ZB(i,k,Z,ZB,ZM,ZM1,p,seed,seq_int);
  }

}

double calculateSP1ZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM){
  if((i<=j)&&(j-i<=SIGMA))
    return 0;
  
  int k;
  //double ZMns = 0;
  double ZMpf = 0;
  double RT = R*T;

  for(k = i; k <= j-SIGMA-2; k++){
    //ZMns += ZM1[j][k] + ZM[k][i]*ZM1[j][k+1];
    ZMpf += ZM1[k][j]*exp(-1*((MultiloopB+MultiloopC*(k-i))/RT)) + ZM[i][k]*ZM1[k+1][j]*exp(-1*(MultiloopB/RT));
  }
  //ZMns += ZM1[j][j-SIGMA-1];
  ZMpf += ZM1[j-SIGMA-1][j]*exp(-1*((MultiloopB+MultiloopC*(j-SIGMA-1-i))/RT));

  ZM[i][j] = ZMpf;

  return 0;//ZMns;
}

int backtrackSP1ZM(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1,  double **p, long *seed, int *seq_int){
   if((i<=j)&&(j-i<=SIGMA))
    return 0;

  int k;
  double RT = R*T;
  int cumulative = 0;  
  int rd = floor(ran2(seed)*MAX_PREC);
  float aux;

  for(k = i; k <= j-SIGMA-1; k++){
    aux = ceil((ZM1[k][j]*exp(-1*((MultiloopB+MultiloopC*(k-i))/RT)))/ZM[i][j]*MAX_PREC);
    cumulative += aux;
    if(rd<cumulative)
      return backtrackSP1ZM1(k,j,Z,ZB,ZM,ZM1,p,seed,seq_int);
  }

  for(k = i; k <= j-SIGMA-1; k++){
    aux = ceil((ZM[i][k]*ZM1[k+1][j]*exp(-1*(MultiloopB/RT)))/ZM[i][j]*MAX_PREC);
    cumulative += aux;
    if(rd<cumulative){
      backtrackSP1ZM(i,k,Z,ZB,ZM,ZM1,p,seed,seq_int);
      backtrackSP1ZM1(k+1,j,Z,ZB,ZM,ZM1,p,seed,seq_int);
      return 0;
    }
  }
}

double calculateSP1Zpfns(int *seq_int, int i, int j, double **ZB, double **Z){
  if(j-i<=SIGMA)
    return 1;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  /*if(cbp[ij_index]==1)
    return 2;
  */

  int k, kj_index;
  //double Zns = Z[j-1][i] + (cbp[ij_index])*ZB[j][i];
  double Zpf = Z[i][j-1] + (cbp[ij_index])*ZB[i][j];
  for(k = i+1; k <= j-SIGMA-1; k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j])];
    //Zns += cbp[kj_index]*Z[k-1][i]*ZB[j][k];
    Zpf += cbp[kj_index]*Z[i][k-1]*ZB[k][j];
  }

  Z[i][j] = Zpf;
 
  return 0;//Zns;
}

int backtrackSP1Z(int i, int j, double **Z, double **ZB, double **ZM, double **ZM1, double **p, long *seed, int *seq_int){
  if(j-i<=SIGMA)
    return 1;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  int k, kj_index;
  float cumulative = 0;  
  int rd = floor(ran2(seed)*MAX_PREC);
  
  cumulative = ceil(Z[i][j-1]/Z[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1Z(i,j-1,Z,ZB,ZM,ZM1,p,seed,seq_int);
   
  cumulative+= ceil(cbp[ij_index]*ZB[i][j]/Z[i][j]*MAX_PREC);
  if(rd<cumulative)
    return backtrackSP1ZB(i,j,Z,ZB,ZM,ZM1,p,seed,seq_int);

  for(k = i+1;k<=j-SIGMA-1;k++){ 
    kj_index = d1[GetIndex(seq_int[k],seq_int[j])];
    cumulative += ceil(cbp[kj_index]*Z[i][k-1]*ZB[k][j]/Z[i][j]*MAX_PREC);
    if(rd<cumulative){
      backtrackSP1Z(i,k-1,Z,ZB,ZM,ZM1,p,seed,seq_int);
      backtrackSP1ZB(k,j,Z,ZB,ZM,ZM1,p,seed,seq_int);
      return 0;
    }
  }

}


double calculateProbs1(int * seq_int, int n, int backtrackSP1s){ 
  double **Z;
  double **ZB;
  double **ZM;
  double **ZM1;
  double **p;
  double Z1n,Zn1;
  int i,j,d;
  long seed = get_local_seed();

  //Allocate MAtrices
  Z = (double**)malloc((n+1)*sizeof(double*));
  ZB = (double**)malloc((n+1)*sizeof(double*));
  ZM = (double**)malloc((n+1)*sizeof(double*));
  ZM1 = (double**)malloc((n+1)*sizeof(double*));
  p = (double**)malloc((n+1)*sizeof(double*));
  for(i=0;i<n+1;i++){
    Z[i] = (double*) calloc(n+1,sizeof(double));
    ZB[i] = (double*) calloc(n+1,sizeof(double));
    ZM[i] = (double*) calloc(n+1,sizeof(double));
    ZM1[i] = (double*) calloc(n+1,sizeof(double));
    p[i] = (double*) calloc(n+1,sizeof(double));
  }

  //Initialize Z
  for(d=0;d<=SIGMA;d++)
    for(i=1;i<=n-d;i++){
      j=i+d;
      Z[i][j] = 1;
      Z[j][i] = 1;
    }

  //Start Recursions
  for (d = SIGMA+1; d <= n; d++){
    for(i=1; (i+d<=n); i++){
      j=i+d;
  
      ZB[j][i] = calculateSP1ZBpfns(seq_int,i,j,ZB,ZM,ZM1);

      ZM1[j][i] = calculateSP1ZM1pfns(seq_int,i,j,ZB,ZM1);

      ZM[j][i] = calculateSP1ZMpfns(seq_int,i,j,ZM1,ZM);
   
      Z[j][i] = calculateSP1Zpfns(seq_int,i,j,ZB,Z);
    }
  }

  //Save results  
  //Z1n = Z[1][n];//partition function
  //Zn1 = Z[n][1];//number of puctures

  //Sampling
  for(i = 0; i<backtrackSP1s;i++){
      backtrackSP1Z(1,n,Z,ZB,ZM,ZM1,p,&seed,seq_int);
  }

  //Calculte probs
  for(i = 1; i<=n;i++)
    for(j=1; j<=n;j++)
      p[i][j] = p[i][j]/backtrackSP1s;

  float aux = 0;
  for(i = 1; i<=n; i++){
    aux = 0;
    for(j=1;j<=n;j++)
      aux += p[i][j];
    p[i][i] = 1 - aux;
  }

  //Print probs
  for(i=1;i<=n;i++){
    for(j=i;j<=n;j++)
      printf("%f ",p[i][j]);
    printf("\n");
  }

  //Free Matrices
  for(i=0; i<n+1; i++){
    free(Z[i]);
    free(ZB[i]);
    free(ZM[i]);
    free(ZM1[i]);
    free(p[i]);
  }
  free(Z);
  free(ZB);
  free(ZM);
  free(ZM1);
  free(p);

  return 0;
}
