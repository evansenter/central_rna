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

double calculateEBB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1){
  if(j-i<=SIGMA+2)
    return INF;

  double RT = R*T;
 
  double EBBmfe;
  double gte = GetTripletEnergy(i,j,i+1,j-1,i+2,j-2,seq_int);
  double min = EBB[i+1][j-1] + gte;

  EBBmfe = /*exp(-1*((*/GetStackEnergy(i,j,i+1,j-1,seq_int)+GetHairpinEnergy(i+1,j-1,seq_int)/*)/RT))*/;
  if(EBBmfe<min)
      min = EBBmfe;

  int k,l,r;
  for(k = i+3; k <= j-SIGMA-3; k++){
    EBBmfe = EB[k][j-2] + /*exp(-1*(*/GetLeftBulgeEnergy(i+1,j-1,k,seq_int)/*/RT))*/;
    if(EBBmfe<min)
      min = EBBmfe;
  }  

  for(k = i+SIGMA+3; k <= j-3; k++){
    EBBmfe = EB[i+1][k] + /*exp(-1*(*/GetRightBulgeEnergy(i+1,j-1,k,seq_int)/*/RT))*/;
    if(EBBmfe<min)
      min = EBBmfe;
  }


  for(l = i+3; (l <= j-SIGMA-4)&&(l-i-2<=2); l++){
    for(r = l+SIGMA+1; (r<=j-3)&&(j-r+l-i-1<=30); r++){
      EBBmfe = EB[l][r] + /*exp(-1*(*/GetILEnergy(i+1,j-1,l,r,seq_int)/*/RT))*/;
      if(EBBmfe<min)
        min = EBBmfe;
    }
  }

  EBBmfe = EM1[i+2][j-2] + /*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
  if(EBBmfe<min)
      min = EBBmfe;

  for(k = i+SIGMA+4; k <= j-SIGMA-3; j++){
    EBBmfe = EM[i+2][k-1] + EM1[k][j-2] + /*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
    if(EBBmfe<min)
      min = EBBmfe;
  }
   
  return min;
}

double calculateEB(int *seq_int, int i, int j, double **EBB, double **EB, double **EM, double **EM1){
  if(j-i<=SIGMA)
    return INF;

  double RT = R*T;

  int k,l,r;
  double min = EBB[i][j];
  double EBmfe;
  
  
  EBmfe = /*exp(-1*(*/GetHairpinEnergy(i,j,seq_int)/*/RT))*/;
  if(/*(cbp[ij_index])&&*/(EBmfe<min)){
    min = EBmfe;
  }
  
  for(k = i+2; k <= j-SIGMA-2; k++){
    EBmfe = EB[k][j-1] + /*exp(-1*(*/GetLeftBulgeEnergy(i,j,k,seq_int)/*/RT))*/;
    if(EBmfe<min)
      min = EBmfe;
  }
  
  
  for(k = i+SIGMA+2; k <= j-2; k++){
    EBmfe = EB[i+1][k] + /*exp(-1*(*/GetRightBulgeEnergy(i,j,k,seq_int)/*/RT))*/;
    if(EBmfe<min)
      min = EBmfe;
  }
  
  for(l = i+2; (l <= j-SIGMA-3)&&(l-i-2<=2); l++){
    for(r = l+SIGMA+1; (r<=j-2)&&(j-r+l-i-1<=30); r++){
      EBmfe = EB[l][r] + /*exp(-1*(*/GetILEnergy(i,j,l,r,seq_int)/*/RT))*/;
      if(EBmfe<min)
        min = EBmfe;
    }
  }
  
  EBmfe = EM1[i+2][j-2]+/*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
  if(EBmfe<min)
      min = EBmfe;

  for(k = i+3; k <= j-SIGMA-2; k++){
    EBmfe = EM[i+1][k-1] + EM1[k][j-1] + /*exp(-1*((*/MultiloopA+MultiloopB/*)/RT))*/;
    if(EBmfe<min)
      min = EBmfe;
  }

  return min;
}

double calculateEM1(int *seq_int, int i, int j, double **EB, double **EM1){
  if(j-i<=SIGMA)
    return INF;
  int k, ik_index;

  double min = INF;
  double EM1mfe = 0;
  double RT = R*T;
  //int bp = -1;

  for(k = i+SIGMA+1; k <= j; k++){
    EM1mfe = EB[i][k] + /*exp(-1*((*/MultiloopC*(j-k)/*)/RT))*/;
    if(EM1mfe < min){
      min = EM1mfe;
      //bp = k;
    } 
  }

  //EM1[j][i] = bp;
  return min;
}

double calculateEM(int *seq_int, int i, int j, double **EM, double **EM1){
  if((i<=j)&&(j-i<=SIGMA))
    return INF;
  
  //int bp = -1;
  int k;
  double min = INF;
  double EMmfe = 0;
  double RT = R*T;

  for(k = i; k <= j-SIGMA-2; k++){
    EMmfe = EM1[k][j] + /*exp(-1*((*/MultiloopB+MultiloopC*(k-i)/*)/RT))*/;
    if(EMmfe < min){
      min = EMmfe;
    //  bp = k;
    }
    EMmfe = EM[i][k] + EM1[k+1][j] + /*exp(-1*(*/MultiloopB/RT/*))*/;
    if(EMmfe < min){
      min = EMmfe;
      //bp = k;
    }
  }
  
  EMmfe = EM1[j-SIGMA-1][j] + /*exp(-1*((*/MultiloopB+MultiloopC*(j-SIGMA-1-i)/*)/RT))*/;
  if(EMmfe < min){
      min = EMmfe;
      //bp = j-SIGMA-1;
    }
  //EM[j][i] = bp;

  return min;
}

double calculateE(int *seq_int, int i, int j, double **EB, double **E){
  if(j-i<=SIGMA)
    return INF;

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

int backtrack(int i, int j, double **E, char *ss){
   int k;
   if (j-i>SIGMA){ 
     k = E[j][i];
     if(k != -1)  {
       ss[k] = '('; 
       ss[j] = ')';
       if( SIGMA <= (j-1)-(k+1))
         backtrack(k+1,j-1,E,ss);
       if (SIGMA <= k-1-i  )
         backtrack(i,k-1,E,ss);
     } 
     else{ // k==-1
       if( SIGMA <= j-1-i ){
          backtrack(i,j-1,E,ss);
          }
       else 
          return 0;
     }
  }  
   
}// endBacktrack


double GetTMFE(int *seq_int, int n){
  
  double **E;
  double **EB;
  double **EBB;
  double **EM;
  double **EM1;
  double E1n,En1;
  int i,j,d;


  //Allocate MAtrices
  E = (double**)malloc((n+1)*sizeof(double*));
  EB = (double**)malloc((n+1)*sizeof(double*));
  EBB = (double**)malloc((n+1)*sizeof(double*));
  EM = (double**)malloc((n+1)*sizeof(double*));
  EM1 = (double**)malloc((n+1)*sizeof(double*));
  for(i=0;i<n+1;i++){
    E[i] = (double*) calloc(n+1,sizeof(double));
    EB[i] = (double*) calloc(n+1,sizeof(double));
    EBB[i] = (double*) calloc(n+1,sizeof(double));
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
      EBB[i][j] = calculateEBB(seq_int,i,j,EBB,EB,EM,EM1);

      EB[i][j] = calculateEB(seq_int,i,j,EBB,EB,EM,EM1);

      EM1[i][j] = calculateEM1(seq_int,i,j,EB,EM1);

      EM[i][j] = calculateEM(seq_int,i,j,EM1,EM);
   
      E[i][j] = calculateE(seq_int,i,j,EB,E);
    }
  }

  //Save results  
  E1n = E[1][n];//minimum free energy

  char* mfe_seq = (char*)malloc((n+1)*sizeof(char));
  for(i=1;i<=n;i++)
    mfe_seq[i] = '.';

  backtrack(1,n,E,mfe_seq);

  for(i=1;i<=n;i++)
    printf("%c",mfe_seq[i]);
  printf("\n");

  //Free Matrices
  for(i=0; i<n+1; i++){
    free(E[i]);
    free(EB[i]);
    free(EBB[i]);
    free(EM[i]);
    free(EM1[i]);
  }
  free(E);
  free(EB);
  free(EBB);
  free(EM);
  free(EM1);
  free(mfe_seq);

  return E1n;
}
