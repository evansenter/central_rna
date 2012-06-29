// #include<stdio.h>
// #include<string.h>
// #include<stdlib.h>
// #include"energy_par.h"
// #include"energy_func.h"
// #include"tpfunc.h"
// #include <math.h>
//  
// //minimum unpaired nts in a loop
// #define SIGMA 3
// #define MIN(x,y) ((x)<=(y) ? (x) : (y) )
// #define MAX(x,y) ((x)>(y)  ? (x) : (y) )
// 
// double calculateZBpfns(int *seq_int, int i, int j, double **ZB, double **ZM, double **ZM1, int n){
//   if(j-i<=SIGMA)
//     return 0;
// 
//   double RT = R*T;
//   int k,l,r,k1,k2,x,y;
//   int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
//   double PIns = 0;
//   double PIpf = 0;
//   double z0ns = 0;
//   double z0pf = 0;
//   double zns = 0;
//   double zpf = 0;
//   double zpf1;
//   int assymetry;
//   int i0,j0,i1,j1,b,b0,ok;
//   int ij0_index, ij1_index, ijb_index;
//   int kj_index, ik_index, lr_index, xy_index;
// 
//   if (i+2+SIGMA+1<=j-2){ 
//     //condition on i,j for there to be possible internal loop
//     /*-------------------------------------------
//      Here we consider 2,3,4 loops; i.e. k1xk2 loops where
//      1<=k1,k2<=3 and k1+k2 in [2,3,4]. Concretely this means all
//      2 loops (1x1), 3 loops (1x2,2x1), 4 loops (1x3,2x2,3x1).
//     -------------------------------------------*/
//     for (k1=1;k1<=3;k1++)
//       for (k2=1;k2<=3;k2++)
//         if (k1+k2==2 || k1+k2==3 || k1+k2==4){
//           x=i+k1+1; y=j-k2-1;
//           xy_index = d1[GetIndex(seq_int[x],seq_int[y])];
//           if (x+SIGMA+1<=y && cbp[xy_index]){ 
//             z0ns += ZB[y][x];
//             z0pf += ZB[x][y]*exp(-cbp[ij_index]*GetILEnergy(i,j,x,y,seq_int)/RT);
//             if((i==1)&&(j==15))
//               printf("IL(1,15,%d,%d) = %f * %f = %f\n",x,y,ZB[x][y],exp(-cbp[ij_index]*GetILEnergy(i,j,x,y,seq_int)/RT),ZB[x][y]*exp(-cbp[ij_index]*GetILEnergy(i,j,x,y,seq_int)/RT));
//           }
//         }
//     PIpf += z0pf;
//     PIns += z0ns; 
//     //ZB[i][j] now contains contributions from all 2,3,4 loops
//     //closed by (i,j).
//     /*-------------------------------------------
//       Now we add to ZB[i][j] the contributions from all
//       k loops, where k=5,6. We then have a for loop with loop
//       control variable b, where we (1) extend, (2) introduce
//       1xb0 and b0x1 loops, where 1+b0=k+2*b. 
// 
//       Problem is that there are different assymetries for 5,6 loops.
//       2x3,3x2: k=5, assymetry=1
//       1x4,4x1: k=5, assymetry=3
//       3x3:     k=6, assymetry=0
//       2x4,4x2: k=6, assymetry=2
//       1x5,5x1: k=6, assymetry=4
//       A subtle point is that we want to add 1xb0 and b0x1 loops,
//       where 1+b0=k+2b only ONCE for each value of k. Thus we will
//       consider 1xb0 and b0x1 loops ONLY for the values
//       assymetry=0  k=6
//       assymetry=1  k=5 
// 
//       Structure of following loop is approximately:
//       for assymetry = 0 to 4
//         for k1 = 1 to 5
//           for k2 = 1 to 5
//             k = k1+k2
//               if abs(k1-k2)==assymetry
//                 x=i+k1+1; y=j-k2-1
//                   for b = 1 to min(i-1,N-j)
//                     b0 = k+2b-1
//                     extend by b closed by (i-b,j+b)
//                     if (assymetry=0 and k=6 or assymetry=1 and k=5)
//                       create 1xb0 and b0x1 loops closed by (i-b,j+b)
//           Code now follows.
//           -------------------------------------------*/
//     for (assymetry=0;assymetry<=4;assymetry++){
//       zns = 0; //reinitialization for new value of assymetry 
//       zpf = 0;
//       for (k1=1;k1<=5;k1++)
//         for (k2=1;k2<=5;k2++){
//           ok = (assymetry==0&&k1==3&&k2==3);
//           ok+= (assymetry==1&&(k1==2&&k2==3 ||k1==3&&k2==2));
//           ok+= (assymetry==2&&(k1==2&&k2==4 ||k1==4&&k2==2));
//           ok+= (assymetry==3&&(k1==1&&k2==4 ||k1==4&&k2==1));
//           ok+= (assymetry==4&&(k1==1&&k2==5 ||k1==5&&k2==1));
//           if ( ok ){
//             k=k1+k2;
//             x=i+k1+1; y=j-k2-1;
//             xy_index = d1[GetIndex(seq_int[x],seq_int[y])];
//             if (x+SIGMA+1<=y && cbp[xy_index]){
//               zns = ZB[y][x];
//               zpf = ZB[x][y]*exp(-cbp[ij_index]*GetILEnergy(i,j,x,y,seq_int)/RT);
//               if((i==1)&&(j==15))
//                 printf("IL(1,15,%d,%d) = %f * %f = %f\n",x,y,ZB[x][y],exp(-cbp[ij_index]*GetILEnergy(i,j,x,y,seq_int)/RT),zpf);
//            
//               PIns += zns;
//               PIpf += zpf;
//            }      
//         for (b=1;b<=MIN(i-1,n-j);b++) {
//         /*-------------------------------------------------------
//            It is ESSENTIAL that one NOT test for basePair(i-b,j+b,a)
//            since such "bogus" computations must be performed for the
//            extension phase.
//            first term -- unnecessary when counting structures,
//            but added here for clarity
//          -------------------------------------------------------*/
//         i0 = i-(b-1); j0 = j+(b-1);
//         //last "rung" in ladder descending (i,j)
//         i1 = i-b ; j1 = j+b ; //next "rung" in ladder descending (i,j)
//         //next "rung" in ladder descending (i,j) at i0-1, j0+1 
//         ij0_index = d1[GetIndex(seq_int[i0],seq_int[j0])];
//         ij1_index = d1[GetIndex(seq_int[i1],seq_int[j1])];
//         zpf  = zpf*exp(cbp[ij0_index]*GetTMMEnergy(i0,j0,i0+1,j0-1,seq_int)/RT);  //remove last rung
//         zpf  = zpf*exp(cbp[ij0_index]*GetILSizeCost(k+2*(b-1))/RT);  //remove last rung
//         zpf  = zpf*exp(cbp[ij0_index]*guau[ij0_index]*AUGUinternal/RT); //remove last rung
//         zpf  = zpf*exp(-cbp[ij1_index]*GetTMMEnergy(i1,j1,i1+1,j1-1,seq_int)/RT); //add new rung
//         zpf  = zpf*exp(-cbp[ij1_index]*GetILSizeCost(k+2*b)/RT); //add new rung
//         zpf  = zpf*exp(-cbp[ij1_index]*guau[ij1_index]*AUGUinternal/RT);
//         //since assymetry is fixed, don't need to account for this
//         if((i1==1)&&(j1==15))
//               printf("from (%d,%d) -- IL(1,15,%d,%d) = %f\n",i,j,x,y,zpf);
//         b0 = k+2*b-1; //size of left/right portion of 1xb0,b0x1 loop
//         /*-------------------------------------------------------
//           second term: (i-b,j+b;x,y) is 1xb0 loop, where 1+b0=k+2*b 
//             only computed when (assymetry=0 and k=6) or when
//             (assymetry=1 and k=5)
//         -------------------------------------------------------*/
//         x  = (i-b)+1+1; y = (j+b)-b0-1;
//         ok = (x+SIGMA+1<=y);
//         xy_index = d1[GetIndex(seq_int[x],seq_int[y])];
//         ok*= (cbp[xy_index]);
//         ok*= (assymetry==0&&k==6||assymetry==1&&k==5);
//         /*-------------------------------------------------------
//           WARNING:do NOT test for (basePair(i-b,j+b,a) && basePair(x,y,a))
//          -------------------------------------------------------*/
//         if (ok){
//           ZB[j1][i1] += ZB[y][x];
//           ZB[i1][j1] += ZB[x][y]*exp(-cbp[ij1_index]*GetILEnergy(i-b,j+b,x,y,seq_int)/RT);
//           if((i1==1)&&(j1==15))
//               printf("from (%d,%d) -- IL(1,15,%d,%d) = %f * %f = %f\n",i,j,x,y,ZB[x][y],exp(-cbp[ij1_index]*GetILEnergy(i-b,j+b,x,y,seq_int)/RT),ZB[x][y]*exp(-cbp[ij1_index]*GetILEnergy(i-b,j+b,x,y,seq_int)/RT));
//         }
//         /*-------------------------------------------------------
//          third term: (i-b,j+b;x,y) is b0x1 loop, where 1+b0=k+2*b 
//          only computed when (assymetry=0 and k=6) or when
//          (assymetry=1 and k=5)
//         -------------------------------------------------------*/
//         x  = (i-b)+b0+1; y = (j+b)-1-1;
//         ok = (x+SIGMA+1<=y);
//         xy_index = d1[GetIndex(seq_int[x],seq_int[y])];
//         ok*= (cbp[xy_index]);
//         ok*= (assymetry==0&&k==6||assymetry==1&&k==5);
//         /*-------------------------------------------------------
//          WARNING:do NOT test for (basePair(i-b,j+b,a) && basePair(x,y,a))
//         -------------------------------------------------------*/
//         if (ok){
//           ZB[j1][i1] += ZB[y][x];
//           ZB[i1][j1] += ZB[x][y]*exp(-cbp[ij1_index]*GetILEnergy(i-b,j+b,x,y,seq_int)/RT);
//           if((i1==1)&&(j1==15))
//               printf("from (%d,%d) -- IL(1,15,%d,%d) = %f * %f = %f\n",i,j,x,y,ZB[x][y],exp(-cbp[ij1_index]*GetILEnergy(i-b,j+b,x,y,seq_int)/RT),ZB[x][y]*exp(-cbp[ij1_index]*GetILEnergy(i-b,j+b,x,y,seq_int)/RT));
//         }
//         //if(cbp[ij1_index]){
//         ZB[j+b][i-b] += zns;
//         ZB[i-b][j+b] += zpf;
//         //}
//       }// end of for loop with loop control variable b 
//             
//           }
//         }//end of for k1,k2 loop
//     }// end of for loop over assymetry values
//   }
// 
//   if(cbp[ij_index]==0){
//     ZB[i][j] += PIpf; 
//     return PIns;
//   }
//   
// 
//   ij_index = d1[GetIndex(seq_int[i+1],seq_int[j-1])];
//   double PSns = cbp[ij_index]*ZB[j-1][i+1];
//   double PSpf = cbp[ij_index]*ZB[i+1][j-1]*exp(-1*((GetStackEnergy(i,j,i+1,j-1,seq_int))/RT));
//   
//   double PHns = 1;
//   double PHpf = exp(-1*(GetHairpinEnergy(i,j,seq_int)/RT));
// 
//   kj_index =  d1[GetIndex(seq_int[i+2],seq_int[j-1])];
//   double PLBns = cbp[kj_index]*ZB[j-1][i+2];
//   double PLBpf = cbp[kj_index]*ZB[i+2][j-1]*exp(-1*((GetStackEnergy(i,j,i+2,j-1,seq_int) + GetLeftBulgeEnergy(i,j,i+2,seq_int))/RT));
// 
//   for(k = i+3; (k <= j-SIGMA-2); k++){
//     kj_index = d1[GetIndex(seq_int[k],seq_int[j-1])];
//     PLBns += cbp[kj_index]*ZB[j-1][k];
//     PLBpf += cbp[kj_index]*ZB[k][j-1]*exp(-1*(GetLeftBulgeEnergy(i,j,k,seq_int)/RT));
//   }
// 
//   ik_index = d1[GetIndex(seq_int[i+1],seq_int[j-2])];
//   double PRBns = cbp[ik_index]*ZB[j-2][i+1];
//   double PRBpf = cbp[ik_index]*ZB[i+1][j-2]*exp(-1*((GetStackEnergy(i,j,i+1,j-2,seq_int) + GetRightBulgeEnergy(i,j,j-2,seq_int))/RT));
//   for(k = j-3; (k>=i+SIGMA+2); k--){
//     ik_index = d1[GetIndex(seq_int[i+1],seq_int[k])];
//     PRBns += cbp[ik_index]*ZB[k][i+1];
//     PRBpf += cbp[ik_index]*ZB[i+1][k]*exp(-1*(GetRightBulgeEnergy(i,j,k,seq_int)/RT));
//   }
// 
//   /*double PIns = 0;
//   double PIpf = 0;
//   for(l = i+2; (l <= j-SIGMA-3)&&(l-i-2<=30); l++){
//     for(r = j -2; (r>=l+SIGMA+1)&&(j-r+l-i-2<=30); r--){
//       lr_index = d1[GetIndex(seq_int[l],seq_int[r])];
//       PIns += cbp[lr_index]*ZB[r][l];
//       PIpf += cbp[lr_index]*ZB[l][r]*exp(-1*(GetILEnergy(i,j,l,r,seq_int)/RT));
//     }
//   }*/
// 
//   double PMns = 0;//ZM1[j-2][i+2];
//   double PMpf = 0;//ZM1[i+2][j-2];
//   for(k = i+SIGMA+3; k <= j-SIGMA-2; k++){
//     PMns += ZM[k-1][i+1]*ZM1[j-1][k];
//     PMpf += ZM[i+1][k-1]*ZM1[k][j-1];
//   }
//   PMpf = PMpf*exp(-1*((MultiloopA+2*MultiloopB)/RT));
// 
// 
//   ZB[i][j] += PSpf + PHpf + PLBpf + PRBpf + PIpf + PMpf;
//   return (PSns+PHns+PLBns+PRBns+PIns+PMns);
// }
// 
// double calculateZM1pfns(int *seq_int, int i, int j, double **ZB, double **ZM1){
//   if(j-i<=SIGMA)
//     return 0;
//   int k, ik_index;
// 
//   double ZM1ns = 0;
//   double ZM1pf = 0;
//   double RT = R*T;
// 
//   for(k = i+SIGMA+1; k <= j; k++){
//     ik_index = d1[GetIndex(seq_int[i],seq_int[k])];
//     ZM1ns += cbp[ik_index]*ZB[k][i];
//     ZM1pf += cbp[ik_index]*ZB[i][k]*exp(-1*((MultiloopC*(j-k))/RT)); 
//   }
// 
//   ZM1[i][j] = ZM1pf;
//   return ZM1ns;
// }
// 
// double calculateZMpfns(int *seq_int, int i, int j, double **ZM1, double **ZM){
//   if((i<=j)&&(j-i<=SIGMA))
//     return 0;
//   
//   int k;
//   double ZMns = 0;
//   double ZMpf = 0;
//   double RT = R*T;
// 
//   for(k = i; k <= j-SIGMA-2; k++){
//     ZMns += ZM1[j][k] + ZM[k][i]*ZM1[j][k+1];
//     ZMpf += ZM1[k][j]*exp(-1*((MultiloopB+MultiloopC*(k-i))/RT)) + ZM[i][k]*ZM1[k+1][j]*exp(-1*(MultiloopB/RT));
//   }
//   ZMns += ZM1[j][j-SIGMA-1];
//   ZMpf += ZM1[j-SIGMA-1][j]*exp(-1*((MultiloopB+MultiloopC*(j-SIGMA-1-i))/RT));
// 
//   ZM[i][j] = ZMpf;
// 
//   return ZMns;
// }
// 
// double calculateZpfns(int *seq_int, int i, int j, double **ZB, double **Z){
//   if(j-i<=SIGMA)
//     return 1;
// 
//   int ij_index = d1[GetIndex(seq_int[i],seq_int[j])];
//   /*if(cbp[ij_index]==1)
//     return 2;
//   */
// 
//   int k, kj_index;
//   double Zns = Z[j-1][i] + (cbp[ij_index])*ZB[j][i];
//   double Zpf = Z[i][j-1] + (cbp[ij_index])*ZB[i][j];
//   for(k = i+1; k <= j-SIGMA-1; k++){
//     kj_index = d1[GetIndex(seq_int[k],seq_int[j])];
//     Zns += cbp[kj_index]*Z[k-1][i]*ZB[j][k];
//     Zpf += cbp[kj_index]*Z[i][k-1]*ZB[k][j];
//   }
// 
//   Z[i][j] = Zpf;
//  
//   return Zns;
// }
// 
// 
// double GetPartFunc(int *seq_int, int n){
//   
//   double **Z;
//   double **ZB;
//   double **ZM;
//   double **ZM1;
//   double Z1n,Zn1;
//   int i,j,d;
// 
// 
//   //Allocate MAtrices
//   Z = (double**)malloc((n+1)*sizeof(double*));
//   ZB = (double**)malloc((n+1)*sizeof(double*));
//   ZM = (double**)malloc((n+1)*sizeof(double*));
//   ZM1 = (double**)malloc((n+1)*sizeof(double*));
//   for(i=0;i<n+1;i++){
//     Z[i] = (double*) calloc(n+1,sizeof(double));
//     ZB[i] = (double*) calloc(n+1,sizeof(double));
//     ZM[i] = (double*) calloc(n+1,sizeof(double));
//     ZM1[i] = (double*) calloc(n+1,sizeof(double));
//   }
// 
//   //Initialize Z
//   for(d=0;d<=SIGMA;d++)
//     for(i=1;i<=n-d;i++){
//       j=i+d;
//       Z[i][j] = 1;
//       Z[j][i] = 1;
//     }
// 
//   //Start Recursions
//   for (d = SIGMA+1; d <= n; d++){
//     for(i=1; (i+d<=n); i++){
//       j=i+d;
//       
//       ZB[j][i] += calculateZBpfns(seq_int,i,j,ZB,ZM,ZM1,n);
// 
//       ZM1[j][i] = calculateZM1pfns(seq_int,i,j,ZB,ZM1);
// 
//       ZM[j][i] = calculateZMpfns(seq_int,i,j,ZM1,ZM);
//    
//       Z[j][i] = calculateZpfns(seq_int,i,j,ZB,Z);
//     }
//   }
// 
//   //Save results  
//   Z1n = Z[1][n];//partition function
//   Zn1 = Z[n][1];//number of structures
// 
//   /* printf("ZBBL:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",ZBBL[i][j]);
//     printf("\n");
//   }
// 
//    printf("ZBBR:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",ZBBR[i][j]);
//     printf("\n");
//   }
// 
//   printf("ZBB:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",ZBB[i][j]);
//     printf("\n");
//   }
//   */
// 
//   /*printf("ZB:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",ZB[i][j]);
//     printf("\n");
//   }
// 
//   printf("ZM1:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",ZM1[i][j]);
//     printf("\n");
//   }
//   
//   printf("ZM:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",ZM[i][j]);
//     printf("\n");
//   }
// 
//   printf("Z:\n");
//   for(i=1;i<=n;i++){
//     for(j=1;j<=n;j++)
//       printf("%f ",Z[i][j]);
//     printf("\n");
//   }*/
// 
//   //Free Matrices
//   for(i=0; i<n+1; i++){
//     free(Z[i]);
//     free(ZB[i]);
//     free(ZM[i]);
//     free(ZM1[i]);
//   }
//   free(Z);
//   free(ZB);
//   free(ZM);
//   free(ZM1);
// 
//   printf("Number of Structures: %f\n",Zn1);
// 
//   return Z1n;
// }
