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

/*double calculateEBBL(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
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

  EBBmfe = EBBL[i+2][j-1] + GetStackEnergy(i,j,i+2,j-1,seq_int) + GetLeftBulgeEnergy(i+2,j-1,i+4,seq_int);
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

  EBBmfe = EBBR[i+2][j-1] + GetStackEnergy(i,j,i+2,j-1,seq_int) + GetRightBulgeEnergy(i+2,j-1,j-3,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -(j-3);
    }


  for(k = i+SIGMA+4; k <= j-3; k++){
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

  //if(EBBmfe<min){
  //    min = EBBmfe;
  //    bp = -1; 
  //}

  for(k = i+SIGMA+5; k <= j-SIGMA-3; k++){
    EBBmfe = EM[i+3][k-1] + EM1[k][j-2] + GetStackEnergy(i,j,i+2,j-1,seq_int) +  MultiloopA+MultiloopB;
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -n - k;
    }
  }

  EBBL[j][i] = bp;

  return min;
}


double calculateEBBR(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
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

  EBBmfe = EBBL[i+1][j-2] + GetStackEnergy(i,j,i+1,j-2,seq_int) + GetLeftBulgeEnergy(i+1,j-2,i+3,seq_int);
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

  EBBmfe = EBBR[i+1][j-2] + GetStackEnergy(i,j,i+1,j-2,seq_int) + GetRightBulgeEnergy(i+1,j-2,j-4,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -(j-4);
    }


  for(k = i+SIGMA+3; k <= j-4; k++){
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

  //if(EBBmfe<min){
  //    min = EBBmfe;
  //    bp = -1; 
  //}

  for(k = i+SIGMA+4; k <= j-SIGMA-4; k++){
    EBBmfe = EM[i+2][k-1] + EM1[k][j-3] + GetStackEnergy(i,j,i+1,j-2,seq_int) +  MultiloopA+MultiloopB;
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -n - k;
    }
  }

  EBBR[j][i] = bp;

  return min;
}


double calculateEBB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1, double **EBBL, double **EBBR, int n){
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

  EBBmfe = EBBL[i+1][j-1] + GetStackEnergy(i,j,i+1,j-1,seq_int) + GetLeftBulgeEnergy(i+1,j-1,i+3,seq_int);
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

  EBBmfe = EBBR[i+1][j-1] + GetStackEnergy(i,j,i+1,j-1,seq_int) + GetRightBulgeEnergy(i+1,j-1,j-3,seq_int);
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -(j-3);
    }


  for(k = i+SIGMA+3; k <= j-3; k++){
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

  //if(EBBmfe<min){
  //    min = EBBmfe;
  //    bp = -1; 
  //}

  for(k = i+SIGMA+4; k <= j-SIGMA-3; k++){
    EBBmfe = EM[i+2][k-1] + EM1[k][j-2] + GetStackEnergy(i,j,i+1,j-1,seq_int) +  MultiloopA+MultiloopB;
    if(EBBmfe<min){
      min = EBBmfe;
      bp = -n - k;
    }
  }

  EBB[j][i] = bp;
   
  return min;
}
*/
double calculateEB(int *seq_int, int i, int j, double **EB, double **EM, double **EM1, int n){
  if(j-i<=SIGMA)
    return INF;

  double RT = R*T;

  int k,l,r;
  double min = INF;

  double EBmfe;
  int bp = 0;
  
  if(j-i>SIGMA+2){
    EBmfe = EB[i+1][j-1]  + GetStackEnergy(i,j,i+1,j-1,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = i;
    }
  }
  
  EBmfe = GetHairpinEnergy(i,j,seq_int);
  if(EBmfe<min){
    min = EBmfe;
    bp = 0;
  }
 
  if(j-1-i-2>SIGMA){
    EBmfe = EB[i+2][j-1] + GetLeftBulgeEnergy(i,j,i+2,seq_int) + GetStackEnergy(i,j,i+2,j-1,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = i+2;
    }  
  }

  for(k = i+3; k <= j-SIGMA-2; k++){
    EBmfe = EB[k][j-1] + GetLeftBulgeEnergy(i,j,k,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = k;
    }
  }
 
  if(j-2-i-1>SIGMA){  
    EBmfe = EB[i+1][j-2] + GetRightBulgeEnergy(i,j,j-2,seq_int) + GetStackEnergy(i,j,i+1,j-2,seq_int);
    if(EBmfe<min){
      min = EBmfe;
      bp = -(j-2);
    }
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

double calculateEM1(int *seq_int, int i, int j, double **EB, double **EM1){
  if(j-i<=SIGMA)
    return INF;
  int k, ik_index;

  /*if((i==21)&&(j==43))
    printf("Here\n");*/

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

double calculateEM(int *seq_int, int i, int j, double **EM1, double **EM){
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

double calculateE(int *seq_int, int i, int j, double **EB, double **E){
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

int backtrackEM(int i, int j, double **EB, double **EM1, double **EM, char *ss, int n){
  int k = EM[j][i];
  if(k>0)//EM1[k,j]
    return backtrackEM1(k,j,EB,EM1,EM,ss,n);
  if(k<0){//EM[i,k] EM1[k+1,j]
    backtrackEM(i,-k,EB,EM1,EM,ss,n);
    backtrackEM1(-k+1,j,EB,EM1,EM,ss,n);
    return 0;
  }

  return 0;
}

int backtrackEM1(int i, int j, double **EB, double **EM1, double **EM, char *ss, int n){
  int k = EM1[j][i];
  if(k<=-1)
    return 0;
  if(k>0)
    return backtrackEB(i,k,EB,EM1,EM,ss,n);
  return 0;
}
/*
int backtrackEBBL(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EBBL[j][i];
  ss[i+2] = '(';
  ss[j-1] = ')';
  if(k==i)//EBB[i][j]
    return backtrackEBB(i+2,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return backtrackEM1(i+2,j-2,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return backtrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-3)
        return backtrackEBBR(i+2,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return backtrackEB(i+3,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      backtrackEM(i+3,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      backtrackEM1(-k-n,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+4)
        return backtrackEBBL(i+2,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return backtrackEB(k,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;

}

int backtrackEBBR(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EBBR[j][i];
  ss[i+1] = '(';
  ss[j-2] = ')';
  if(k==i)//EBB[i][j]
    return backtrackEBB(i+1,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return backtrackEM1(i+2,j-2,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return backtrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-4)
        return backtrackEBBR(i+1,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return backtrackEB(i+2,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      backtrackEM(i+2,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      backtrackEM1(-k-n,j-3,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+3)
        return backtrackEBBL(i+1,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return backtrackEB(k,j-3,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;

}

int backtrackEBB(int i, int j, double **EB, double **EBB, double **EM1, double **EM, double **EBBL, double **EBBR, char *ss, int n){
  int k = EBB[j][i];
  ss[i+1] = '(';
  ss[j-1] = ')';
  if(k==i)//EBB[i][j]
    return backtrackEBB(i+1,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return backtrackEM1(i+2,j-2,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return backtrackEB(k%n,k/n,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-3)
        return backtrackEBBR(i+1,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return backtrackEB(i+2,-k,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      backtrackEM(i+2,-k-n-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      backtrackEM1(-k-n,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+3)
        return backtrackEBBL(i+1,j-1,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
      return backtrackEB(k,j-2,EB,EBB,EM1,EM,EBBL,EBBR,ss,n);
    }
  }
  return 0;

}

*/
int backtrackEB(int i, int j, double **EB, double **EM1, double **EM, char *ss, int n){
  int k = EB[j][i];
  ss[i] = '(';
  ss[j] = ')';
  if(j-i<=SIGMA+2)
    return 0;

  if(k==i)//EBB[i][j]
    return backtrackEB(i+1,j-1,EB,EM1,EM,ss,n);
  if(k==0)//Hairpin
    return 0;
  //if(k==-1)//EM1[i+1][j-1]
  //  return backtrackEM1(i+1,j-1,EB,EBB,EM1,EM,ss,n);
  if(k>n)//IL
    return backtrackEB(k%n,k/n,EB,EM1,EM,ss,n);
  if(k<-1){
    if(k>=-n){//RB
      if(-k==j-2)
        return backtrackEB(i+1,j-2,EB,EM1,EM,ss,n);
      return backtrackEB(i+1,-k,EB,EM1,EM,ss,n);
    }
    if(k<-n){//EM[i+1,k-1]+EM1[k,j-1]
      backtrackEM(i+1,-k-n-1,EB,EM1,EM,ss,n);
      backtrackEM1(-k-n,j-1,EB,EM1,EM,ss,n);
      return 0;
    }
  }else{
    if(k<=n){//LB
      if(k==i+2)
        return backtrackEB(i+2,j-1,EB,EM1,EM,ss,n);
      return backtrackEB(k,j-1,EB,EM1,EM,ss,n);
    }
  }
  return 0;
  
}

int backtrackE(int i, int j, double **E, double **EB, double **EM1, double **EM, char *ss, int n){
   int k;
   if (j-i>SIGMA){ 
     k = E[j][i];
     if(k != -1)  {
       if(k==i)
         return backtrackEB(i,j,EB,EM1,EM,ss,n);
       /*if( SIGMA <= (j-1)-(k+1))
         backtrackEB(k+1,j-1,EB,ss);
       if (SIGMA <= k-1-i  )
         backtrackEB(i,k-1,EB,ss);*/
       else{
         backtrackE(i,k-1,E,EB,EM1,EM,ss,n);
         backtrackEB(k,j,EB,EM1,EM,ss,n);
         return 0;
       }
     } 
     else{ // k==-1
       if( SIGMA <= j-1-i ){
          return backtrackE(i,j-1,E,EB,EM1,EM,ss,n);
          }
       else 
          return 0;
     }
  }
  return 0;  
   
}// endBacktrack


void GetMFE(int *seq_int, int n, char * mfe_seq){
  
  double **E;
  double **EB;
  double **EM;
  double **EM1;
  double E1n,En1;
  int i,j,d;


  //Allocate MAtrices
  E = (double**)malloc((n+1)*sizeof(double*));
  EB = (double**)malloc((n+1)*sizeof(double*));
  EM = (double**)malloc((n+1)*sizeof(double*));
  EM1 = (double**)malloc((n+1)*sizeof(double*));
  for(i=0;i<n+1;i++){
    E[i] = (double*) calloc(n+1,sizeof(double));
    EB[i] = (double*) calloc(n+1,sizeof(double));
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

      EB[i][j] = calculateEB(seq_int,i,j,EB,EM,EM1,n);

      EM1[i][j] = calculateEM1(seq_int,i,j,EB,EM1);

      EM[i][j] = calculateEM(seq_int,i,j,EM1,EM);
   
      E[i][j] = calculateE(seq_int,i,j,EB,E);
    }
  }

  //Save results  
  E1n = E[1][n];//minimum free energy

/*  char* mfe_seq = (char*)malloc((n+1)*sizeof(char));
  for(i=1;i<=n;i++)
    mfe_seq[i] = '.';*/

  backtrackE(1,n,E,EB,EM1,EM,mfe_seq,n);


  for(i=1;i<=n;i++)
    printf("%c",mfe_seq[i]);
  printf("\n");

  /*printf("E:\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      printf("%f ",E[i][j]);
    printf("\n");
  }
  printf("EB:\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      printf("%f ",EB[i][j]);
    printf("\n");
  }
  printf("EM:\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      printf("%f ",EM[i][j]);
    printf("\n");
  }
  printf("EM1:\n");
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++)
      printf("%f ",EM1[i][j]);
    printf("\n");
  }*/


  //Free Matrices
  for(i=0; i<n+1; i++){
    free(E[i]);
    free(EB[i]);
    free(EM[i]);
    free(EM1[i]);
  }
  free(E);
  free(EB);
  free(EM);
  free(EM1);
  //free(mfe_seq);

  printf("\n%f\n",E1n);
  //return E1n;
}
