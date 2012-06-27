#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"energy_par.h"
#include"energy_func.h"
#include"tpfunc.h"
#include <math.h>
 
//minimum unpaired nts in a loop
#define SIGMA 3

double calculateZBBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1){
  if(j-i<=SIGMA+2)
    return 0;

  double RT = R*T;
 
  int ij_index = d1[GetIndex(seq_int[i+2],seq_int[j-2])];
  int kj_index, ik_index,lr_index;
  double QSns = cbp[ij_index]*ZBB[j-1][i+1];
  double gte = GetTripletEnergy(i,j,i+1,j-1,i+2,j-2,seq_int);
  double egte = exp(-1*(gte/RT));
  double QSpf = cbp[ij_index]*ZBB[i+1][j-1]*egte;

  double QHns = 1;
  double QHpf = exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int)+GetHairpinEnergy(i+1,j-1,seq_int))/RT));

  double QLBns = 0;
  double QLBpf = 0;
  int k,l,r;
  for(k = i+3; k <= j-SIGMA-3; k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-2])];
    QLBns += cbp[kj_index]*ZB[j-2][k];
    QLBpf += cbp[kj_index]*ZB[k][j-2]*exp(-1*(GetLeftBulgeEnergy(i+1,j-1,k,seq_int)/RT));
  }  

  double QRBns = 0;
  double QRBpf = 0;
  for(k = i+SIGMA+3; k <= j-3; k++){
    ik_index = d1[GetIndex(seq_int[i+2],seq_int[k])];
    QRBns += cbp[ik_index]*ZB[k][i+2];
    QRBpf += cbp[ik_index]*ZB[i+1][k]*exp(-1*(GetRightBulgeEnergy(i+1,j-1,k,seq_int)/RT));
  }

  double QIns = 0;
  double QIpf = 0;

  for(l = i+3; (l <= j-SIGMA-4)&&(l-i-2<=2); l++){
    for(r = l+SIGMA+1; (r<=j-3)&&(j-r+l-i-1<=30); r++){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      QIns += cbp[lr_index]*ZB[r][l];
      QIpf += cbp[lr_index]*ZB[l][r]*exp(-1*(GetILEnergy(i+1,j-1,l,r,seq_int)/RT));
    }
  }
  double QMns = ZM1[j-2][i+2];
  double QMpf = ZM1[i+2][j-2];
  for(k = i+SIGMA+4; k <= j-SIGMA-3; j++){
    QMns += ZM[k-1][i+2]*ZM1[j-2][k];
    QMpf += ZM[i+2][k-1]*ZM1[k][j-2];
  }
  QMpf = QMpf*exp(-1*((MultiloopA+MultiloopB)/RT));
 
  ZBB[i][j] = QSpf + QHpf + QLBpf + QRBpf + QIpf + QMpf;

  return (QSns+QHns+QLBns+QRBns+QIns+QMns);
}

double calculateZBpfns(int *seq_int, int i, int j, double **ZBB, double **ZB, double **ZM, double **ZM1){
  if(j-i<=SIGMA)
    return 0;

  double RT = R*T;

  int k,l,r;
  int ij_index = d1[GetIndex(seq_int[i+1],seq_int[j-1])];
  int kj_index, ik_index, lr_index;
  double PSns = cbp[ij_index]*ZBB[j][i];
  double PSpf = cbp[ij_index]*ZBB[i][j];
  
  double PHns = 1;
  double PHpf = exp(-1*(GetHairpinEnergy(i,j,seq_int)/RT));

  double PLBns = 0;
  double PLBpf = 0;
  for(k = i+2; k <= j-SIGMA-2; k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j-1])];
    PLBns += cbp[kj_index]*ZB[j-1][k];
    PLBpf += cbp[kj_index]*ZB[k][j-1]*exp(-1*(GetLeftBulgeEnergy(i,j,k,seq_int)/RT));
  }
  double PRBns = 0;
  double PRBpf = 0;
  for(k = i+SIGMA+2; k <= j-2; k++){
    ik_index = d1[GetIndex(seq_int[i+1],seq_int[k])];
    PRBns += cbp[ik_index]*ZB[k][i+1];
    PRBpf += cbp[ik_index]*ZB[i+1][k]*exp(-1*(GetRightBulgeEnergy(i,j,k,seq_int)/RT));
  }
  double PIns = 0;
  double PIpf = 0;
  for(l = i+2; (l <= j-SIGMA-3)&&(l-i-2<=2); l++){
    for(r = l+SIGMA+1; (r<=j-2)&&(j-r+l-i-1<=30); r++){
      lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
      PIns += cbp[lr_index]*ZB[r][l];
      PIpf += cbp[lr_index]*ZB[l][r]*exp(-1*(GetILEnergy(i,j,l,r,seq_int)/RT));
    }
  }
  double PMns = ZM1[j-2][i+2];
  double PMpf = ZM1[i+2][j-2];
  for(k = i+3; k <= j-SIGMA-2; k++){
    PMns += ZM[k-1][i+1]*ZM1[j-1][k];
    PMpf += ZM[i+1][k-1]*ZM1[k][j-1];
  }
  PMpf = PMpf*exp(-1*((MultiloopA+MultiloopB)/RT));


  ZB[i][j] = PSpf + PHpf + PLBpf + PRBpf + PIpf + PMpf;
  return (PSns+PHns+PLBns+PRBns+PIns+PMns);
}

double calculateZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1){
  if(j-i<=SIGMA)
    return 0;
  int k, ik_index;

  double ZM1ns = 0;
  double ZM1pf = 0;
  double RT = R*T;

  for(k = i+SIGMA+1; k <= j; k++){
    ik_index = d1[GetIndex(seq_int[i],seq_int[k])];
    ZM1ns += cbp[ik_index]*ZB[k][i];
    ZM1pf += cbp[ik_index]*ZB[i][k]*exp(-1*((MultiloopC*(j-k))/RT)); 
  }

  ZM1[i][j] = ZM1pf;
  return ZM1ns;
}

double calculateZMpfns(int *seq_int, int i, int j, double **ZM, double **ZM1){
  if((i<=j)&&(j-i<=SIGMA))
    return 0;
  
  int k;
  double ZMns = 0;
  double ZMpf = 0;
  double RT = R*T;

  for(k = i; k <= j-SIGMA-2; k++){
    ZMns += ZM1[j][k] + ZM[k][i]*ZM1[j][k+1];
    ZMpf += ZM1[k][j]*exp(-1*((MultiloopB+MultiloopC*(k-i))/RT)) + ZM[i][k]*ZM1[k+1][j]*exp(-1*(MultiloopB/RT));
  }
  ZMns += ZM1[j][j-SIGMA-1];
  ZMpf += ZM1[j-SIGMA-1][j]*exp(-1*((MultiloopB+MultiloopC*(j-SIGMA-1-i))/RT));

  ZM[i][j] = ZMpf;

  return ZMns;
}

double calculateZpfns(int *seq_int, int i, int j, double **ZB, double **Z){
  if(j-i<=SIGMA)
    return 1;

  int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
  /*if(cbp[ij_index]==1)
    return 2;
  */

  int k, kj_index;
  double Zns = Z[j-1][i] + (cbp[ij_index])*ZB[j][i];
  double Zpf = Z[i][j-1] + (cbp[ij_index])*ZB[i][j];
  for(k = i+1; k <= j-SIGMA-1; k++){
    kj_index = d1[GetIndex(seq_int[k],seq_int[j])];
    Zns += cbp[kj_index]*Z[k-1][i]*ZB[j][k];
    Zpf += cbp[kj_index]*Z[i][k-1]*ZB[k][j];
  }

  Z[i][j] = Zpf;
 
  return Zns;
}


double GetPartFunc(int *seq_int, int n){
  
  double **Z;
  double **ZB;
  double **ZBB;
  double **ZM;
  double **ZM1;
  double Z1n,Zn1;
  int i,j,d;


  //Allocate MAtrices
  Z = (double**)malloc((n+1)*sizeof(double*));
  ZB = (double**)malloc((n+1)*sizeof(double*));
  ZBB = (double**)malloc((n+1)*sizeof(double*));
  ZM = (double**)malloc((n+1)*sizeof(double*));
  ZM1 = (double**)malloc((n+1)*sizeof(double*));
  for(i=0;i<n+1;i++){
    Z[i] = (double*) calloc(n+1,sizeof(double));
    ZB[i] = (double*) calloc(n+1,sizeof(double));
    ZBB[i] = (double*) calloc(n+1,sizeof(double));
    ZM[i] = (double*) calloc(n+1,sizeof(double));
    ZM1[i] = (double*) calloc(n+1,sizeof(double));
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
      ZBB[j][i] = calculateZBBpfns(seq_int,i,j,ZBB,ZB,ZM,ZM1);

      ZB[j][i] = calculateZBpfns(seq_int,i,j,ZBB,ZB,ZM,ZM1);

      ZM1[j][i] = calculateZM1pfns(seq_int,i,j,ZB,ZM1);

      ZM[j][i] = calculateZMpfns(seq_int,i,j,ZM1,ZM);
   
      Z[j][i] = calculateZpfns(seq_int,i,j,ZB,Z);
    }
  }

  //Save results  
  Z1n = Z[1][n];//partition function
  Zn1 = Z[n][1];//number of structures

  //Free Matrices
  for(i=0; i<n+1; i++){
    free(Z[i]);
    free(ZB[i]);
    free(ZBB[i]);
    free(ZM[i]);
    free(ZM1[i]);
  }
  free(Z);
  free(ZB);
  free(ZBB);
  free(ZM);
  free(ZM1);

  printf("Number of Structures: %f\n",Zn1);

  return Z1n;
}
