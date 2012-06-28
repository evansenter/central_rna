#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"energy_par.h"
#include"energy_func.h"
#include"tmfe.h"
#include <math.h>
#include <limits.h>
 
//minimum unpaired nts in a loop
#define SIGMA 3
#define INF INT_MAX/2

double tcalculateEBBL(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
  if(j-i<=SIGMA+3)
    return INF;


  double RT = R*T;

  double EBBmfe;
  double gte = GetTripletEnergy(i,j,i+2,j-1,i+3,j-2,seq_int);
  double min = EBB[i+2][j-1] + gte;
  int bp = i;

  EBBmfe = GetStackEnergy(i,j,i+2,j-1,seq_int)+GetHairpinEnergy(i+2,j-1,seq_int);
  if(EBBmfe<min){
      min = EBBmfe;
      bp = 0;
  }

  EBBmfe = EBBL[i+2][j-1] + GetTripletEnergy(i,j,i+2,j-1,i+4,j-2,seq_int) + GetLeftBulgeEnergy(i+2,j-1,i+4,seq_int);
  if(EBBmfe<min){
    min = EBBmfe;
    bp = i+4;
  }

  int k,l,r;
  for(k = i+5; k <= j-SIGMA-3; k++){
    EBBmfe = EB[k][j-2] + GetStackEnergy(i,j,i+2,j-1,seq_int) + GetLeftBulgeEnergy(i+2,j-1,k,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = k;
    }
  }

  EBBmfe = EBBR[i+2][j-1] + GetTripletEnergy(i,j,i+2,j-1,i+3,j-3,seq_int) + GetRightBulgeEnergy(i+2,j-1,j-3,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -(j-3);
    }


  for(k = i+SIGMA+4; k <= j-4; k++){
    EBBmfe = EB[i+3][k] + GetStackEnergy(i,j,i+2,j-1,seq_int) + GetRightBulgeEnergy(i+2,j-1,k,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -k;
    }
  }


 for(l = i+4; (l <= j-SIGMA-4)&&(l-i-1<=30); l++){
    for(r = j-3; (r>=l+SIGMA+1)&&(j-r+l-i-1<=30); r--){
      EBBmfe = EB[l][r] + GetStackEnergy(i,j,i+2,j-1,seq_int) + GetILEnergy(i+2,j-1,l,r,seq_int);
      if(EBBmfe<min){
        min = EBBmfe;
        bp = r*n + l;
      }
    }
  }

  //EBBmfe = EM1[i+2][j-2] + /*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
  //if(EBBmfe<min){
  //    min = EBBmfe;
  //    bp = -1; 
  //}

  for(k = i+SIGMA+5; k <= j-SIGMA-3; k++){
    EBBmfe = EM[i+3][k-1] + EM1[k][j-2] + GetStackEnergy(i,j,i+2,j-1,seq_int) +  MultiloopA+2*MultiloopB;
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -n - k;
    }
  }

  EBBL[j][i] = bp;

  return min;
}


double tcalculateEBBR(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
  if(j-i<=SIGMA+3)
    return INF;

  double RT = R*T;

  double EBBmfe;
  double gte = GetTripletEnergy(i,j,i+1,j-2,i+2,j-3,seq_int);
  double min = EBB[i+1][j-2] + gte;
  int bp = i;

  EBBmfe = GetStackEnergy(i,j,i+1,j-2,seq_int)+GetHairpinEnergy(i+1,j-2,seq_int);
  if(EBBmfe<min){
      min = EBBmfe;
      bp = 0;
  }

  EBBmfe = EBBL[i+1][j-2] + GetTripletEnergy(i,j,i+1,j-2,i+3,j-3,seq_int) + GetLeftBulgeEnergy(i+1,j-2,i+3,seq_int);
  if(EBBmfe<min){
    min = EBBmfe;
    bp = i+3;
  }

  int k,l,r;
  for(k = i+4; k <= j-SIGMA-4; k++){
    EBBmfe = EB[k][j-3] + GetStackEnergy(i,j,i+1,j-2,seq_int) + GetLeftBulgeEnergy(i+1,j-2,k,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = k;
    }
  }

  EBBmfe = EBBR[i+1][j-2] + GetTripletEnergy(i,j,i+1,j-2,i+2,j-4,seq_int) + GetRightBulgeEnergy(i+1,j-2,j-4,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -(j-4);
    }


  for(k = i+SIGMA+3; k <= j-5; k++){
    EBBmfe = EB[i+2][k] + GetStackEnergy(i,j,i+1,j-2,seq_int) + GetRightBulgeEnergy(i+1,j-2,k,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -k;
    }
  }

  for(l = i+3; (l <= j-SIGMA-5)&&(l-i-2<=30); l++){
    for(r = j-4; (r>=l+SIGMA+1)&&(j-r+l-i-3<=30); r--){
      EBBmfe = EB[l][r] + GetStackEnergy(i,j,i+1,j-2,seq_int) + GetILEnergy(i+1,j-2,l,r,seq_int);
      if(EBBmfe<min){
        min = EBBmfe;
        bp = r*n + l;
      }
    }
  }

  //EBBmfe = EM1[i+2][j-2] + /*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
  //if(EBBmfe<min){
  //    min = EBBmfe;
  //    bp = -1; 
  //}

  for(k = i+SIGMA+4; k <= j-SIGMA-4; k++){
    EBBmfe = EM[i+2][k-1] + EM1[k][j-3] + GetStackEnergy(i,j,i+1,j-2,seq_int) +  MultiloopA+2*MultiloopB;
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -n - k;
    }
  }

  EBBR[j][i] = bp;

  return min;
}


double tcalculateEBB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
  if(j-i<=SIGMA+2)
    return INF;


  double RT = R*T;
 
  double EBBmfe;
  double gte = GetTripletEnergy(i,j,i+1,j-1,i+2,j-2,seq_int);
  double min = EBB[i+1][j-1] + gte;
  int bp = i;

  EBBmfe = GetStackEnergy(i,j,i+1,j-1,seq_int)+GetHairpinEnergy(i+1,j-1,seq_int);
  if(EBBmfe<min){
      min = EBBmfe;
      bp = 0;
  }

  EBBmfe = EBBL[i+1][j-1] + GetTripletEnergy(i,j,i+1,j-1,i+3,j-2,seq_int) + GetLeftBulgeEnergy(i+1,j-1,i+3,seq_int);
  if(EBBmfe<min){
    min = EBBmfe;
    bp = i+3;
  }

  int k,l,r;
  for(k = i+4; k <= j-SIGMA-3; k++){
    EBBmfe = EB[k][j-2] + GetStackEnergy(i,j,i+1,j-1,seq_int) + GetLeftBulgeEnergy(i+1,j-1,k,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = k;
    }
  }  

  EBBmfe = EBBR[i+1][j-1] + GetTripletEnergy(i,j,i+1,j-1,i+2,j-3,seq_int) + GetRightBulgeEnergy(i+1,j-1,j-3,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -(j-3);
    }


  for(k = i+SIGMA+3; k <= j-4; k++){
    EBBmfe = EB[i+2][k] + GetStackEnergy(i,j,i+1,j-1,seq_int) + GetRightBulgeEnergy(i+1,j-1,k,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -k;
    }
  }


  for(l = i+3; (l <= j-SIGMA-4)&&(l-i-2<=30); l++){
    for(r = j-3; (r>=l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
      EBBmfe = EB[l][r] + GetStackEnergy(i,j,i+1,j-1,seq_int) + GetILEnergy(i+1,j-1,l,r,seq_int);
      if(EBBmfe<min){
        min = EBBmfe;
        bp = r*n + l;
      }
    }
  }

  //EBBmfe = EM1[i+2][j-2] + /*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
  //if(EBBmfe<min){
  //    min = EBBmfe;
  //    bp = -1; 
  //}

  for(k = i+SIGMA+4; k <= j-SIGMA-3; k++){
    EBBmfe = EM[i+2][k-1] + EM1[k][j-2] + GetStackEnergy(i,j,i+1,j-1,seq_int) +  MultiloopA+2*MultiloopB;
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -n - k;
    }
  }

  EBB[j][i] = bp;
   
  return min;
}

double tcalculateEB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
  if(j-i<=SIGMA)
    return INF;

  double RT = R*T;

  int k,l,r;
  double min = EBB[i][j];
  double EBmfe;
  int bp = i;
  
  EBmfe = GetHairpinEnergy(i,j,seq_int);
  if(EBmfe<min){
    min = EBmfe;
    bp = 0;
  }
 
  EBmfe = EBBL[i][j] + GetLeftBulgeEnergy(i,j,i+2,seq_int);
  if(EBmfe<min){
    min = EBmfe;
    bp = i+2;
  } 

  for(k = i+3; k <= j-SIGMA-2; k++){
    EBmfe = EB[k][j-1] + GetLeftBulgeEnergy(i,j,k,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = k;
    }
  }
  
  EBmfe = EBBR[i][j] + GetRightBulgeEnergy(i,j,j-2,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = -(j-2);
    }

  for(k = i+SIGMA+2; k <= j-3; k++){
    EBmfe = EB[i+1][k] + GetRightBulgeEnergy(i,j,k,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = -k;
    }
  }
  
  for(l = i+2; (l <= j-SIGMA-3)&&(l-i-2<=30); l++){
    for(r = j-2; (r>= l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
      EBmfe = EB[l][r] + GetILEnergy(i,j,l,r,seq_int);
      if(EBmfe<min){
        min = EBmfe;
        bp = r*n + l;
      }
    }
  }
  
  /*EBmfe = EM1[i+2][j-2]+MultiloopA+MultiloopB;
  if(EBmfe<min){
      min = EBmfe;
      bp = -1;
  }*/

  for(k = i+SIGMA+3; k <= j-SIGMA-2; k++){
    EBmfe = EM[i+1][k-1] + EM1[k][j-1] + MultiloopA+2*MultiloopB;
    if(EBmfe<min){
      min = EBmfe;
      bp = -n - k;
    }
  }

  EB[j][i] = bp;

  return min;
}

double tcalculateEM1(int *seq_int, int i, int j, double **EB, double **EM1){
  if(j-i<=SIGMA)
    return INF;
  int k, ik_index;

  double min = INF;
  double EM1mfe = 0;
  double RT = R*T;
  int bp = -1;

  for(k = i+SIGMA+1; k <= j; k++){
    EM1mfe = EB[i][k] + MultiloopC*(j-k);
    if(EM1mfe < min){
      min = EM1mfe;
      bp = k;
    } 
  }

  EM1[j][i] = bp;
  return min;
}

double tcalculateEM(int *seq_int, int i, int j, double **EM1, double **EM){
  if((i<=j)&&(j-i<=SIGMA))
    return INF;

  int bp = -1;
  int k;
  double min = INF;
  double EMmfe = 0;
  double RT = R*T;

  for(k = i; k <= j-SIGMA-2; k++){
    EMmfe = EM1[k][j] + MultiloopB+MultiloopC*(k-i);
    if(EMmfe < min){
      min = EMmfe;
      bp = k;
    }
    EMmfe = EM[i][k] + EM1[k+1][j] + MultiloopB;
    if(EMmfe < min){
      min = EMmfe;
      bp = -k;
    }
  }
  
  EMmfe = EM1[j-SIGMA-1][j] + MultiloopB+MultiloopC*(j-SIGMA-1-i);
  if(EMmfe < min){
      min = EMmfe;
      bp = j-SIGMA-1;
    }
  EM[j][i] = bp;

  return min;
}

double tcalculateE(int *seq_int, int i, int j, double **EB, double **E){
  if(j-i<=SIGMA)
    return 0;

  int k;
  double min = E[i][j-1];//, kj_index;
  double Emfe = E[i][j-1];
  int bp = -1;
  if(EB[i][j]<min){
    min = EB[i][j];
    bp = i;
  }
  
  for(k = i+1; k <= j-SIGMA-1; k++){
    Emfe = E[i][k-1] + EB[k][j];
    if(Emfe<min){
      min = Emfe;
      bp = k;
    }
  }

  E[j][i] = bp;
 
  return min;
}

int tbacktrackEM(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EM[j][i];
  if(k>0)//EM1[k,j]
    return tbacktrackEM1(k,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<0){//EM[i,k] EM1[k+1,j]
    tbacktrackEM(i,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    tbacktrackEM1(-k+1,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    return 0;
  }

  return 0;
}

int tbacktrackEM1(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EM1[j][i];
  if(k<=-1)
    return 0;
  if(k>0)
    return tbacktrackEB(i,k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  return 0;
}

int tbacktrackEBBL(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EBBL[j][i];
  ss[i+2] = '(';
  ss[j-1] = ')';
  if(k==i)//EBB[i][j]
    return tbacktrackEBB(i+2,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return tbacktrackEM1(i+2,j-2,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return tbacktrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-3)
        return tbacktrackEBBR(i+2,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(i+3,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      tbacktrackEM(i+3,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      tbacktrackEM1(-k-n,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+4)
        return tbacktrackEBBL(i+2,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(k,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;

}

int tbacktrackEBBR(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EBBR[j][i];
  ss[i+1] = '(';
  ss[j-2] = ')';
  if(k==i)//EBB[i][j]
    return tbacktrackEBB(i+1,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return tbacktrackEM1(i+2,j-2,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return tbacktrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-4)
        return tbacktrackEBBR(i+1,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(i+2,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      tbacktrackEM(i+2,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      tbacktrackEM1(-k-n,j-3,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+3)
        return tbacktrackEBBL(i+1,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(k,j-3,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;

}

int tbacktrackEBB(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EBB[j][i];
  ss[i+1] = '(';
  ss[j-1] = ')';
  if(k==i)//EBB[i][j]
    return tbacktrackEBB(i+1,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return tbacktrackEM1(i+2,j-2,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return tbacktrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-3)
        return tbacktrackEBBR(i+1,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(i+2,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      tbacktrackEM(i+2,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      tbacktrackEM1(-k-n,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+3)
        return tbacktrackEBBL(i+1,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(k,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;

}


int tbacktrackEB(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EB[j][i];
  ss[i] = '(';
  ss[j] = ')';
  if(k==i)//EBB[i][j]
    return tbacktrackEBB(i,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return tbacktrackEM1(i+1,j-1,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return tbacktrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-2)
        return tbacktrackEBBR(i,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(i+1,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      tbacktrackEM(i+1,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      tbacktrackEM1(-k-n,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+2)
        return tbacktrackEBBL(i,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return tbacktrackEB(k,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;
  
}

int tbacktrackE(int i, int j, double **E, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
   int k;
   if (j-i>SIGMA){ 
     k = E[j][i];
     if(k != -1)  {
       if(k==i)
         return tbacktrackEB(i,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
       /*if( SIGMA <= (j-1)-(k+1))
         tbacktrackEB(k+1,j-1,EB,ss);
       if (SIGMA <= k-1-i  )
         tbacktrackEB(i,k-1,EB,ss);*/
       else{
         tbacktrackE(i,k-1,E,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
         tbacktrackEB(k,j,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
         return 0;
       }
     } 
     else{ // k==-1
       if( SIGMA <= j-1-i ){
          return tbacktrackE(i,j-1,E,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
          }
       else 
          return 0;
     }
  }
  return 0;  
   
}// endBacktrack


void GetTMFE(int *seq_int, int n, char *mfe_seq){
  
  double **E;
  double **EB;
  double **EBB;
  double **EBBL;
  double **EBBR;
  double **EM;
  double **EM1;
  double E1n,En1;
  int i,j,d;


  //Allocate MAtrices
  E = (double**)malloc((n+1)*sizeof(double*));
  EB = (double**)malloc((n+1)*sizeof(double*));
  EBB = (double**)malloc((n+1)*sizeof(double*));
  EBBL = (double**)malloc((n+1)*sizeof(double*));
  EBBR = (double**)malloc((n+1)*sizeof(double*));
  EM = (double**)malloc((n+1)*sizeof(double*));
  EM1 = (double**)malloc((n+1)*sizeof(double*));
  for(i=0;i<n+1;i++){
    E[i] = (double*) calloc(n+1,sizeof(double));
    EB[i] = (double*) calloc(n+1,sizeof(double));
    EBB[i] = (double*) calloc(n+1,sizeof(double));
    EBBL[i] = (double*) calloc(n+1,sizeof(double));
    EBBR[i] = (double*) calloc(n+1,sizeof(double));
    EM[i] = (double*) calloc(n+1,sizeof(double));
    EM1[i] = (double*) calloc(n+1,sizeof(double));
  }

  //Initialize E
  for(i=1;i<=n;i++)
    for(j=i+1;j<=n;j++){
      E[j][i] = -1;
    }

  //Start Recursions
  for (d = SIGMA+1; d <= n; d++){
    for(i=1; (i+d<=n); i++){
      j=i+d;
      EBBL[i][j] = tcalculateEBBL(seq_int,i,j,EBB,EB,EM,EM1,EBBL,EBBR,n);

      EBBR[i][j] = tcalculateEBBR(seq_int,i,j,EBB,EB,EM,EM1,EBBL,EBBR,n);

      EBB[i][j] = tcalculateEBB(seq_int,i,j,EBB,EB,EM,EM1,EBBL,EBBR,n);

      EB[i][j] = tcalculateEB(seq_int,i,j,EBB,EB,EM,EM1,EBBL,EBBR,n);

      EM1[i][j] = tcalculateEM1(seq_int,i,j,EB,EM1);

      EM[i][j] = tcalculateEM(seq_int,i,j,EM1,EM);
   
      E[i][j] = tcalculateE(seq_int,i,j,EB,E);
    }
  }

  //Save results  
  E1n = E[1][n];//minimum free energy

  /*char* mfe_seq = (char*)malloc((n+2)*sizeof(char));
  for(i=1;i<=n;i++)
    mfe_seq[i] = '.';
  mfe_seq[n+1] = '\0';
*/
  tbacktrackE(1,n,E,EB,EBB,EM1,EM,EBBL,EBBR,mfe_seq,n);
  for(i=1;i<=n;i++)
    printf("%c",mfe_seq[i]);
  printf("\n");

  /*for(i=1;i<=n;i++){
    for(j=i+1;j<=n;j++)
      printf("%f ",E[i][j]);
    printf("\n");
  }*/

  //Free Matrices
  for(i=0; i<n+1; i++){
    free(E[i]);
    free(EB[i]);
    free(EBB[i]);
    free(EBBL[i]);
    free(EBBR[i]);
    free(EM[i]);
    free(EM1[i]);
  }
  free(E);
  free(EB);
  free(EBB);
  free(EBBL);
  free(EBBR);
  free(EM);
  free(EM1);
  //free(mfe_seq);

  printf("\n%f\n",E1n);
  //return mfe_seq;
}
