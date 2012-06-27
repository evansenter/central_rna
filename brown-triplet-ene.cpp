/* brown-triplet.c
modified from brown-3pl.c
by Valya Ilyin
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NST 0     /* Energy for nonstandard stacked pairs */
#define PUBLIC
#define GASCONST 1.98717  /* in [cal/K] */
#define K0  273.15
#define INF 1000000
#define FORBIDDEN 9999
#define BONUS 10000
#define NBPAIRS 7


PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
/*           CG     GC     GU     UG     AU     UA      */
{ {   INF ,  INF ,  INF ,  INF ,  INF ,  INF ,  INF , INF} ,
  {   INF , -330 , -240 , -140 , -210 , -210 , -210 , NST} , /* CG */
  {   INF , -340 , -330 , -150 , -250 , -240 , -220,  NST} , /* GC */
  {   INF , -250 , -210 , -50 ,   130 , -130 , -140 , NST} , /* GU */
  {   INF , -150 , -140 ,  30 ,   -50 , -100 ,  -60 , NST} , /* UG */
  {   INF , -220 , -210 , -60 ,  -140 , -90 ,  -110,  NST} , /* AU */
  {   INF , -240 , -210 , -100 , -130 , -130 ,  -90 , NST} , /* UA */
  {   INF ,  NST ,  NST ,  NST ,  NST ,  NST ,  NST , NST}};
/*     0,     1,     2,     3,     4,     5,     6,   7  */


/* enthalpies (0.01*kcal/mol at 37 C) for stacked pairs */
/* different from mfold-2.3, which uses values from mfold-2.2 */
PUBLIC int enthalpies[NBPAIRS+1][NBPAIRS+1] = 
/*          CG     GC     GU     UG     AU     UA  */
{ {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF}, 
  {  INF, -1060, -1340, -1210,  -560, -1050, -1040, NST},
  {  INF, -1340, -1490, -1260,  -830, -1140, -1240, NST},
  {  INF, -1210, -1260, -1460, -1350,  -880, -1280, NST},
  {  INF,  -560,  -830, -1350,  -930,  -320,  -700, NST},
  {  INF, -1050, -1140,  -880,  -320,  -940,  -680, NST},
  {  INF, -1040, -1240, -1280,  -700,  -680,  -770, NST},
  {  INF,   NST,   NST,   NST,   NST,   NST,   NST, NST}};



#define M 3
#define K 216
    /*---------------------------------------------------------
    M-tuples of the 6 possible nucleotide pairs, making the probability 
    distribution p[] and target distribution q[] range over 6^M=K = 216
    values.
    array : p[7][7][7] (0 means nonstandard pair; 0 is not used here)
    example: p[1][1][1] = CG-CG-CG or 5'-CCC-3'/ 3'-GGG-5'
    
    marginals:
    ql[7][7], qr[7][7]; target marginal are computed from stack37 tables; 0 is not used
    ---------------------------------------------------------*/



#define BETA(x,i) (( (x) >> (2*(i)) )%NBPAIRS)
	// Goedel-style BETA function
	// giving i-th coordinate of sequence encoded by x


double Z;


void convertEenergyToProbabilities(double mp[NBPAIRS][NBPAIRS],
				   double mpr[NBPAIRS][NBPAIRS],
				   double mpm[NBPAIRS][NBPAIRS]) {

   double eij[NBPAIRS][NBPAIRS];
   double ejk[NBPAIRS][NBPAIRS];
   double eik[NBPAIRS][NBPAIRS];
   double Zstar = 0.0;
   double Zrstar = 0.0;
   double Zmstar = 0.0;

   double Esum=0.0;
   double Ersum=0.0;
   double Emsum=0.0;

   for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) 
       mp[i][j] = eij[i][j] = ejk[i][j] = eik[i][j] = 0;

   double den = GASCONST*(K0+37);  // 20 degrees C

   for (int i=1; i<NBPAIRS; i++) {
     for (int j=1; j<NBPAIRS; j++) {
        double e = 0.0;
	double er = 0.0; 
	double em=0.0;

	//printf ("%d %d %f\n", i,j, e);
	for (int k=1; k<NBPAIRS; k++) {
	  e = stack37[i][j] + stack37[j][k];
	  eij[i][j] += exp(-e/den);
	  er =  stack37[i][j] + stack37[k][i];
	  ejk[i][j] += exp(-er/den);
	  em = stack37[i][k] + stack37[k][j];
	  eik[i][j] += exp(-em/den);
	}
	Zstar  += eij[i][j];
	Zrstar += ejk[i][j];
	Zmstar += eik[i][j];
     }
   }

   for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) { 
       mp[i][j] = eij[i][j]/Zstar;
       mpr[i][j] = ejk[i][j]/Zrstar;
       mpm[i][j] = eik[i][j]/Zmstar;
     }
 
    
    Z=Zmstar;
   //printf("Zstar= %7.5e, Zrstar= %7.5e, Zm= %7.5e\n", Zstar, Zrstar, Zmstar);
   //printf("E= %7.5e, Er= %7.5e, Em= %7.5e\n", Esum, Ersum, Emsum);

}


void calcMarginals(double trp[NBPAIRS][NBPAIRS][NBPAIRS],
              double ql[NBPAIRS][NBPAIRS], 
	      double qr[NBPAIRS][NBPAIRS], 
	      double qm[NBPAIRS][NBPAIRS]) {

    double tsum=0.0;
    double msum=0.0;
    for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) {
       double lsum = 0;
       double msumk = 0;
       for (int k=1; k<NBPAIRS; k++) {
         lsum += trp[i][j][k];
	 msumk += trp[i][k][j];
       }
       ql[i][j]=lsum;
       tsum +=lsum;
       qm[i][j]=msumk;
       msum +=msumk;
     }
     printf("tsum ql= %7.5e, qm= %7.5e, ",tsum, msum);
    
    tsum=0.0;
    for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) {
	double rsum = 0;
       for (int k=1; k<NBPAIRS; k++) {
         
	 rsum += trp[k][i][j];
       }
       
       qr[i][j]=rsum;
       tsum += rsum;
     }
     printf(" qr: %7.5e\n",tsum);

     
     return;
}

void compareMarginals(double m2[NBPAIRS][NBPAIRS], double q2[NBPAIRS][NBPAIRS]) {

     double diff=0.0;
     double ratio=0.0;
     double max=0.0;
     int count=0;
     for (int i=1; i<NBPAIRS; i++) 
      for (int j=1; j<NBPAIRS; j++) {
        double d = abs(m2[i][j]-q2[i][j]);
	if (d > max) max =d;
        diff += d;
	ratio += m2[i][j]/q2[i][j];
	count++;
      }
     double adiff = diff/count;
     double aratio = ratio/count;
     printf("Compare results: \
             Total diff=%7.5e, Total ratio=%7.5e, adiff=%7.5e, aratio=%7.5e, max=%7.5e\n",
             diff, ratio, adiff,aratio, max);
}


void printMarginals(double m2[NBPAIRS][NBPAIRS]) {

     double sum=0.0;
    for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) {
       printf("%d %d %7.5f\n", i, j, m2[i][j]);
       sum += m2[i][j];
     }
    printf("sum=%7.5e\n", sum);
}

void printM3(double m3[NBPAIRS][NBPAIRS][NBPAIRS]) {

    for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) 
       for (int k=1; k<NBPAIRS; k++)  
         printf("%d %d %d %7.5e\n", i, j, k, m3[i][j][k]);

}


int main(int argc, char *argv[]) {

    if (argc != 2) {
       fprintf(stderr,"%s numberIterations\n",argv[0]);
       
       return 0;
    }

    int count=0;
    int numIter = atoi(argv[1]);
 

    double mp[NBPAIRS][NBPAIRS]; // duplet marginals
    double mpr[NBPAIRS][NBPAIRS]; // duplet marginals
    double mpm[NBPAIRS][NBPAIRS]; // duplet marginals

    double ql[NBPAIRS][NBPAIRS]; // left marginals from distribution 
    double qr[NBPAIRS][NBPAIRS]; // right marginals from distribution 
    double qm[NBPAIRS][NBPAIRS]; // right marginals from distribution 
 

    convertEenergyToProbabilities(mp,mpr,mpm);
    printMarginals(mp);
    printMarginals(mpr);
    printMarginals(mpm);

    double trp[NBPAIRS][NBPAIRS][NBPAIRS];

    // initialize triplet uniform probablity distribution
    for (int i=1; i<NBPAIRS; i++) 
      for (int j=1; j<NBPAIRS; j++) 
        for (int k=1; k<NBPAIRS; k++) trp[i][j][k] = (double) 1/K; 

     printf("Printing initial triple distribution\n");
     printM3(trp);


    // initialize marginals to 0    
    for (int i=1; i<NBPAIRS; i++) 
     for (int j=1; j<NBPAIRS; j++) {
       ql[i][j] = 0;
       qr[i][j] = 0;
     }

    calcMarginals(trp,ql,qr,qm);
    
    for (count=0;count<numIter;count++){
	// update according to Brown's algorithm
      for (int i=1; i<NBPAIRS; i++) 
        for (int j=1; j<NBPAIRS; j++) 
          for (int k=1; k<NBPAIRS; k++) {
	     printf("i=%d, j=%d, k=%d, ql=%7.5e, mpij=%7.5e, trp=%7.5e, ",
	            i, j, k, ql[i][j], mp[i][j], trp[i][j][k]);
	     trp[i][j][k] *= mp[i][j]/ql[i][j];
	     printf(" new trp=%7.5e\n", trp[i][j][k]);
	}

       calcMarginals(trp,ql,qr,qm);
       // update according to Brown's algorithm
      for (int i=1; i<NBPAIRS; i++) 
        for (int j=1; j<NBPAIRS; j++) 
          for (int k=1; k<NBPAIRS; k++) {
	     
	     trp[k][i][j] *= mpr[i][j]/qr[i][j];
	     
	  }
	calcMarginals(trp,ql,qr,qm);

      for (int i=1; i<NBPAIRS; i++) 
        for (int j=1; j<NBPAIRS; j++) 
          for (int k=1; k<NBPAIRS; k++) {
	     
	     trp[i][k][j] *= mpm[i][j]/qm[i][j];
	     
	  }



      calcMarginals(trp,ql,qr,qm);
      // check sum	  
      double onesum=0.0;
      for (int i=1; i<NBPAIRS; i++) 
        for (int j=1; j<NBPAIRS; j++) 
          for (int k=1; k<NBPAIRS; k++) onesum += trp[i][j][k];
	  
       printf("Cycle: %d (1sum=%5.3e), Printing distribution\n", count, onesum);
       compareMarginals(mp,ql);
       compareMarginals(mpr,qr);
       compareMarginals(mpm,qm);
       //printM3(trp);
     }	

     printM3(trp);
     double den = GASCONST*(K0+37.0);  // 20 degrees C
     //energy
     printf("Energy for triplets:\n");
      for (int i=1; i<NBPAIRS; i++) 
        for (int j=1; j<NBPAIRS; j++) 
          for (int k=1; k<NBPAIRS; k++) {
	     double eout=0.0;
	     eout = -den*log(Z*trp[i][j][k]);
	     printf("i=%d, j=%d, k=%d, eout=%7.1f, sum2=%7.1f\n", 
	             i,j,k,eout, (double) (stack37[i][j]+stack37[j][k]));
          }

     return 0;
}
