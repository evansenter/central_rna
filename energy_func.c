#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include"energy_par.h"
#include"energy_func.h"
#include"enthalpy_par.h"
#include"energy_constant.h"
DEBUG = 0;
/* Energy functions according to Matthews Rules located at
   http://rna.urmc.rochester.edu/NNDB/

*/

/* GetTripletEnergy

   Input:
   1) position of i
   2) position of j closing i
   3) position of k
   4) position of l closing k
   5) position of u
   6) position of v closing u
   7) sequence of integers representing nucleotides in original input sequence

   output:
   1) corresponding triplet energy


   triplet is modeled as follows:
   3'  5'
     uv
     kl
     ij
   5'  3' 
*/

double GetTripletEnergy(int i, int j, int k, int l, int u, int v, int * trans_seq)
{ 
  int ij_trans_index, kl_trans_index, uv_trans_index, ij_stack_index, kl_stack_index, ijkl_stack_index, uv_stack_index;
  double energy; 
  
  // ij_trans_index and kl_trans_index provide the indices to look up within d1 that correspond to ij and kl 
  
  ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// ij_trans_index corresponds to the index location of ij in d1 
  kl_trans_index = GetIndex(trans_seq[k],trans_seq[l]);// kl_trans_index corresponds to the index location of kl in d1 
  uv_trans_index = GetIndex(trans_seq[u],trans_seq[v]);// uv_trans_index corresponds to the index location of uv in d1
  
  // ij_stack_index and kl_stack_index provide the indices ij(row) and kl(column) in stack37
  
  ij_stack_index = d1[ij_trans_index];// ij_stack_index is the index of ij in energy tables located in energy_par.c
  kl_stack_index = d1[kl_trans_index];// kl_stack_index is the index of kl in energy tables located in energy_par.c
  ijkl_stack_index = ij_stack_index*16 + kl_stack_index; // ijkl_stack_index is the index of ijkl located in row of triplet37 in energy_par.c
  uv_stack_index = d1[uv_trans_index];// uv_stack_index is the index of uv located in column of triplet37 in energy_par.c

  energy = triplet[ijkl_stack_index][uv_stack_index];// energy value is assigned as the one at the intersect of ijkl_stack_index(row) and uv_stack_index(column) in triplet37 in energy_par.c

  if(DEBUG){
    printf("Triplet: (%d,%d) %c%c; (%d,%d) %c%c; (%d,%d) %c%c: %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),u,v,GetNucl(trans_seq[u]),GetNucl(trans_seq[v]),energy);
  }
  return energy;
}

/* GetStackEnergy

   Input:
   1) position of i
   2) position of j closing i
   3) position of k
   4) position of l closing k
   5) sequence of integers representing nucleotides in original input sequence

   output:
   1) corresponding stacking energy


   stacking is modeled as follows:
   3'  5'
     kl
     ij
   5'  3' 
*/
double GetStackEnergy(int i, int j, int k, int l, int * trans_seq)
{
  int ij_trans_index, kl_trans_index, ij_stack_index, kl_stack_index;
  double energy;

  // ij_trans_index and kl_trans_index provide the indices to look up within d1 that correspond to ij and kl 
  
  ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// ij_trans_index corresponds to the index location of ij in d1
  kl_trans_index = GetIndex(trans_seq[k],trans_seq[l]);// kl_trans_index corresponds to the index location of kl in d1
 
  // ij_stack_index and kl_stack_index provide the indices ij(row) and kl(column) in stack37

  ij_stack_index = d1[ij_trans_index];// ij_stack_index is the index of ij in energy tables located in energy_par.c
  kl_stack_index = d1[kl_trans_index];// kl_stack_index is the index of kl in energy tables located in energy_par.c

  energy = stack[ij_stack_index][kl_stack_index];// energy value is assigned as intersection of ij_stack_index(row) and kl_stack_index(column) in stack37 in energy_par.c


  if(DEBUG){
    printf("Stacking: (%d,%d) %c%c; (%d,%d) %c%c: %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),energy);
  }
  return energy;  
}

/*

  GetDangleEnergy

  Input:

  1) position of i
  2) position of j closing i
  3) position of k; 0 if 5' dangle
  4) position of l; 0 if 3' dangle

  Output:
  1) dangle energy based on input


  Dangles are interior not EXTERIOR:

 5' dangle

  3' 5'
   0l
   ij
  5' 3'

 3' dangle

  3' 5'
   k0
   ij
  5' 3'
 
*/

double GetDangleEnergy(int i, int j, int k, int l, int * trans_seq)
{
  int ij_trans_index, kl_trans_index, ij_dangle_index, kl_dangle_index;
  double energy;

  // dangle in table is interior, this is for when dangle is exterior  
  // input is reverse from way table is set up  
  if (i == 0 || j == 0){
    return GetDangleEnergy(l,k,j,i,trans_seq);
  }  
  else{
    ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// need to get index to look up ij in d1
    if(k==0){
      kl_trans_index = GetIndex(0,trans_seq[l]);// need to get index to look up 5' dangle in d1
    }
    else if(l==0){
      kl_trans_index = GetIndex(trans_seq[k],0);// need to get index to look up 3' dangle in d1
    }
    ij_dangle_index = d1[ij_trans_index];// ij_dangle_index corresponds to index of ij(row) in dangle37 located in energy_par.c
    kl_dangle_index = d1[kl_trans_index];// kl_dangle_index corresponds to index of dangle in dangle37 located in energy_par.c

    energy = dangle[ij_dangle_index][kl_dangle_index];// dangle energy corresponds to intersect of ij_dangle_index(row) and kl_dangle_index(column) in dangle37
    
    if (DEBUG){
      if(k == 0){      
        printf("3' Dangle: (%d,%d) %c%c; (%d,%d) 0%c: %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[l]),energy);
      }
      else if (l==0){
        printf("5' Dangle: (%d,%d) %c%c; (%d,%d) %c0: %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),energy);
      } 
    }

    return energy;
  }
}

/* GetTMMEnergy

   input:
   1) position of i in integer_seq
   2) position of j closing i
   3) position of k which is i+1 
   4) position of l which is j-1
   5) integer sequence correspond to original sequence

   output:
   1) Terminal mismatch energy

   TMM modeled as follows, interior:
 
   3'  5'
     kl
     ij
   5'  3'

   * ij is closing base pair, kl is not base paired but a mismatch. 

*/
double GetTMMEnergy(int i, int j, int k, int l, int * trans_seq)
{
  int ij_trans_index, kl_trans_index, ij_tmm_index, kl_tmm_index;
  double energy;

  ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// ij_trans_index corresponds to index of ij in d1
  kl_trans_index = GetIndex(trans_seq[k],trans_seq[l]);// kl_trans_index corresponds to index of kl in d1

  ij_tmm_index = d1[ij_trans_index];// ij_tmm_index corresponds to index of ij(row) necessary to look up in TMM located in energy_par.c
  kl_tmm_index = d1[kl_trans_index];// kl_tmm_index corresponds to index of kl(column) necessary to look up in TMM located in energy_par.c; kl is not a basepair but the Terminal mismatch

  energy = TMM[ij_tmm_index][kl_tmm_index];// energy corresponds to intersect of ij_tmm_index(row_) and kl_tmm_index(column) located in TMM in energy_par.c
  
  if (DEBUG){
    printf("TMM: (%d,%d) %c%c; (%d,%d) %c%c: %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),energy);
  }
  
  return energy;
}


/*

  GetBulgeEnergy
  
  Input:
  1) position of i in input sequence
  2) position of j closing i in input sequence
  3) position of k in input sequence
  4) position of l closing k in input sequence

  Output:
  1) corresponding bulge energy

  Bulge energy is modeled as follows:

  3' 5'        3' 5'
   kl           kl
 .        or      .
   ij           ij
  5' 3'        5' 3'

  AU/GU penalty in incorporated in energy such that if ij is AU or GU a penalty is added to the bulge cost

*/
double GetBulgeEnergy(int i, int j, int k, int l, int u, int v, int* trans_seq, int triplet_flag)
{
  int diffi, diffj;// diffi is the difference between k and i, diffj is the difference between j and l; a value of 1 means there is no nucleotide between the positions

  double bulge_energy,stack_energy;
  bulge_energy = 0.0; 
  bulge_energy += Initiation; // Intermolecular initiation 

  diffi = k-i;
  diffj = j-l;

  if(DEBUG){
   printf("\n In Bulge: %d %d %d %d %d %d\n",i,j,k,l,u,v);
  }

  if (diffi == 1){// there is no nucleotide in between i and k so bulge is between l and j
    
    // bulge_cost + special c bulge + stacking energy between i,j and k,l
    if(diffj-1 <31){
      bulge_energy+=bulge[diffj-1];// # of nucleotides in between l and j is diffj -1, so not inclusive of k or j; bulge37 is a 1D integer array where index corresponds to size of bulge and value stored at index is cost.
    }
    else{//if(diffj-1<31,size of bulge is greater than 30
      bulge_energy+=bulge[6]+1.75*R*T*log((diffj-1)/6.0);
    }
    if(DEBUG){
      printf("Bulge: size: (%d) bulge energy: %f\n",diffj-1,bulge_energy);
    }
    // if bulge just size of 1 base
    if(diffj == 2){// remember that diffj-1 corresponds to bulge of size 1; must check if nucleotide inbetween l and j is C
      if(trans_seq[l+1] == 7){ // check to see if its a C
        bulge_energy+= C_bulge;// special penatly for c bulge, C_bulge is located in energy_par.c
        if(DEBUG){
           printf("C penalty: %f\n, bulge energy %f\n",C_bulge,bulge_energy);
        }
      }
      // add stacking base pair energy
      if (triplet_flag && (u-k-1 + l-v-1) >= 0 && (u-k-1 +l-v-1) <= 1){
        stack_energy=GetTripletEnergy(i,j,k,l,u,v,trans_seq);// for bulge of size one must also add the stacking basepair of ij with k
        bulge_energy+=stack_energy;
        if(DEBUG){           
          printf("Triplet for Bulge: triplet value: %f, bulge energy: %f\n",stack_energy, bulge_energy);
        }
      }
      else{
        stack_energy=GetStackEnergy(i,j,k,l,trans_seq);// for bulge of size one must also add the stacking basepair of ij with kl
        bulge_energy+=stack_energy;
        if(DEBUG){
          printf("Stacking for Bulge: stack value: %f bulge energy: %f\n",stack_energy,bulge_energy);
        }
      }
    }
    // check for AU or GU penalty at end of helix
    else if(trans_seq[j] == 12){// i is U
     if (trans_seq[i] == 1 || trans_seq[i] == 3){// if i is A or G then there is a penalty
       bulge_energy+= GUAU_penalty;// addition of penalty; GUAU_penalty is in energy_par.c
       if(DEBUG){
         printf("GU/AU penalty: (%d,%d) %c%c penalty: %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),GUAU_penalty);
       }
     }
    }
 
  }

  else{// if(diffi == 1);  there is no nucleotide between l and j so bulge is between i and k, same procedure as before when bulge was between l and j 

    // bulge_cost + special c bulge + stacking energy between i,j and k,l
    if(diffi-1 <31){
      bulge_energy+=bulge[diffi-1];// # of nucleotides in between l and j is diffj -1, so not inclusive of k or j; bulge37 is a 1D integer array where index corresponds to size of bulge and value stored at index is cost.
    }
    else{//if(diffj-1<31,size of bulge is greater than 30
      bulge_energy+=bulge[6]+1.75*R*T*log((diffi-1)/6.0);
    }
    if(DEBUG){
      printf("Bulge: (size %d) bulge energy: %f\n",diffi-1,bulge_energy);
    }
    // if bulge just size of 1 base
    if(diffi == 2){
      if(trans_seq[i+1] == 7){// check to see if its a C
        bulge_energy+= C_bulge;
        if(DEBUG){
           printf("C penalty: %f,  bulge energy: %f\n",C_bulge,bulge_energy);
        }
      }
      // add stacking base pair energy
      if (triplet_flag && (u-k-1+l-v-1) >= 0 && (u-k-1+l-v-1) <= 1){
        stack_energy=GetTripletEnergy(i,j,k,l,u,v,trans_seq);// for bulge of size one must also add the stacking basepair of ij with k
        bulge_energy+=stack_energy;
        if(DEBUG){
          printf("Triplet for Bulge: triplet value %f: bulge energy: %f\n",stack_energy,bulge_energy);
        }
      }
      else{  
        stack_energy=GetStackEnergy(i,j,k,l,trans_seq);// for bulge of size one must also add the stacking basepair of ij with kl
        bulge_energy+=stack_energy;
        if(DEBUG){
          printf("In bulge\n");
          printf("Stacking for Bulge: stack value %f: bulge energy: %f\n",stack_energy,bulge_energy);
        }
      }
    }
    else if(trans_seq[j] == 12){// i is U
     if (trans_seq[i] == 1 || trans_seq[i] == 3){// if i is A or G then there is a penalty
       bulge_energy+= GUAU_penalty;// addition of penalty; GUAU_penalty is in energy_par.c
       if(DEBUG){
         printf("GU/AU penalty: (%d,%d) %c%c penalty %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),GUAU_penalty);
       }
     }
    }
   
  }

  if(DEBUG){
    printf("Total Bulge Energy: %f\n",bulge_energy);   
  }
  return bulge_energy;
}

/* GetRightBulgeEnergy

   input:
   1) position i
   2) position j closing i
   3) position v
   4) integer representation of sequence

   output:
   1) Right bulge energy

   Modeled as follows:

   5'   3'
      v
     ij
   3'   5'  

   v-j represents right bulge; i+1 is basepaired to v

*/

double GetRightBulgeEnergy(int i, int j, int l, int * trans_seq)
{
  int diffj;
  double bulge_energy, stack_energy;
  bulge_energy = 0.0;
  bulge_energy += Initiation;
  diffj = j-l;
  
  if(DEBUG){
    printf("\n");
  }

  // bulge_cost + special c bulge
  if(diffj-1 <31){
      bulge_energy+=bulge[diffj-1];// # of nucleotides in between l and j is diffj -1, so not inclusive of k or j; bulge37 is a 1D integer array where index corresponds to size of bulge and value stored at index is cost.
    }
  else{//if(diffj-1<31,size of bulge is greater than 30
      bulge_energy+=bulge[6]+1.75*R*T*log((diffj-1)/6.0);
  }
  if(DEBUG){
    printf("Bulge: (%d) %f\n",diffj-1,bulge_energy);
  }
  // if bulge just size of 1 base
  if(diffj == 2){
    if(trans_seq[l+1] == 7){ // check to see if its a C
      bulge_energy+= C_bulge;
      if(DEBUG){
         printf("C penalty: %f\n",C_bulge);
      }
    }
  }
  else if(trans_seq[j] == 12){// i is U
     if (trans_seq[i] == 1 || trans_seq[i] == 3){// if i is A or G then there is a penalty
       bulge_energy+= GUAU_penalty;// addition of penalty; GUAU_penalty is in energy_par.c
       if(DEBUG){
         printf("GU/AU penalty: (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),GUAU_penalty);
       }
     }
  }

  if(DEBUG){
    printf("\n");
  }

  return bulge_energy;
}

/*
  GetLeftBulgeEnergy

  Input:
  1) position of i in sequence
  2) position of j closing i
  3) position of k in sequence
  4) integer representation of sequence

  Output:
  1) Left Bulge energy

  Modeled as follows

  3'   5'
   k
   i j
  5'   3'

  i->k represents the left bulge, one knows k is base paired to j-1

*/

double GetLeftBulgeEnergy(int i, int j, int k, int * trans_seq)
{

    int diffi;

    double bulge_energy, stack_energy; 
    bulge_energy = 0.0;
    bulge_energy += Initiation; 
    diffi = k-i;


    if(DEBUG){
      printf("\n");
    }


    // bulge_cost + special c bulge 
    if(diffi-1 <31){
      bulge_energy+=bulge[diffi-1];// # of nucleotides in between l and j is diffj -1, so not inclusive of k or j; bulge37 is a 1D integer array where index corresponds to size of bulge and value stored at index is cost.
    }
    else{//if(diffi-1<31,size of bulge is greater than 30
      bulge_energy+=bulge[6]+1.75*R*T*log((diffi-1)/6.0);
    }
    if(DEBUG){
      printf("Bulge: (%d) %f\n",diffi-1,bulge_energy);
    }
    // if bulge just size of 1 base
    if(diffi == 2){
      if(trans_seq[i+1] == 7){// check to see if its a C
        bulge_energy+= C_bulge;
        if(DEBUG){
           printf("C penalty: %f\n",C_bulge);
        }
      }
    }

    else if(trans_seq[j] == 12){// i is U
     if (trans_seq[i] == 1 || trans_seq[i] == 3){// if i is A or G then there is a penalty
       bulge_energy+= GUAU_penalty;// addition of penalty; GUAU_penalty is in energy_par.c
       if(DEBUG){
         printf("GU/AU penalty: (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),GUAU_penalty);
       }
     }
    }

  
    return bulge_energy;
}

/*

  GetInternal11
  
  input:
  1) position of i in seqeunce
  2) position of j closing i
  3) position of k in sequence
  4) position of l closing k
  5) integer converted sequence
 
  output:
  1) Internal Loop Energy of 1x1 internal loop
  

  Internal Loop 1x1 modeled as follows:

 3' 5'  
  k l
 x   y
  i j
 5'  3'

 ij and kl are basepairs
 xy corresponds to the 1x1 internal loop, neither x or y are basepaired

*/
double GetInternal11(int i, int j, int k, int l, int * trans_seq){

  double energy;
  int xy_trans_index,xy_index,ij_trans_index,kl_trans_index,ij_index,kl_index,ijkl;

  // obtain id of 1x1 loop
  xy_trans_index = GetIndex(trans_seq[i+1],trans_seq[l+1]);// xy_trans_index corresponds to the index of xy in d1
  xy_index = d1[xy_trans_index]; // xy_index corresponds to the index of xy to look up in the internal11 table in energy_par.c
  
  // i,j,k,l
  ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// ij_trans_index corresponds to the index of ij in d1
  kl_trans_index = GetIndex(trans_seq[k],trans_seq[l]);// kl_trans_index corresponds to the index of kl in d1

  ij_index = d1[ij_trans_index];// ij_index corresponds to the index of ij necessary for lookup in energy tables
  kl_index = d1[kl_trans_index];// kl_index corresponds to the index of kl necessary for looup in energy tables

  ijkl = linearize2D(ij_index,kl_index);// *6 b/c there are 6 base pairs, linearize the 2D stack array using row columnar multiplication to get the row for internal11: For example to get an index for ijkl one would take the ij_index * 6 to get the row and kl_index to know how many columns to move over. This value provides the row to look up in internal11

  energy = internal11[ijkl][xy_index];// corresponds to the 1x1 internal loop energy of base pairs ij and kl represented by ij_index(row) with a 1x1 internal loop identified by xy with xy_index(column) 
  energy += Initiation;
  if(DEBUG){
    printf("\n");
    printf("Internal Loop 1x1: (%d,%d) %c%c to (%d,%d) %c%c, Loop (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),i,l,GetNucl(trans_seq[i+1]),GetNucl(trans_seq[l+1]),energy);
    printf("\n");
  }    

  return energy;
}

/* GetInternal12

   Input:
   1) position of i in sequence
   2) position of j closing i
   3) position of k in sequence
   4) position of l closing k in sequence
   5) integer representation of sequence

   output:
   1) Internal 1x2 energy

   Modeled as follows; important to realize that there are modifications taken into account if internal loop is 2x1 with 2 nucleotides between i and k and 1 between l and j:
  
  3'  5'
    kl  
   0  b
   x  y
    ij
   5' 3'

*/

double GetInternal12(int i, int j, int k, int l, int * trans_seq){


  int ijkl, ij_trans_index, kl_trans_index, ij_index, kl_index;
  int x,y,xy,xy_index, b, ab , ab_index, xy_ab;
  
  double energy;

  /*
    k l
   0   b 
   x   y 
    i j

  */

  if (k-i == 2){// between i and k there is one nucleotide so there are two nucleotides between k and l
  
    // get i component

    ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// corresponds to index for ij in d1
    kl_trans_index = GetIndex(trans_seq[k],trans_seq[l]);// corresponds to index for kl in d1

    ij_index = d1[ij_trans_index];// corresponds to index for ij in energy tables
    kl_index = d1[kl_trans_index];// corresponds to index for kl in energy tables

    ijkl = linearize2D(ij_index,kl_index);// row columnar multilplication used to obtain index of ijkl(row) used in GetInternal12

    x = trans_seq[i+1];// x corresponds to integer representation of nucleotide at i+1
    y = trans_seq[l+2];// y corresponds to integer representation of nucleotide at l+2
  
    xy = GetIndex(x,y);// xy corresponds to index of xy in d1
    xy_index = d1[xy];// xy_index corresponds to index of xy in energy tables
  
    b = trans_seq[l+1];// b corresponds to integer representation of nucleotide at l+1
    ab = GetIndex(0,b);// ab corresponds to index of ab in d1; here a is 0 b/c this a 1x2 internal loop where there is 1 un base paired nucleotide between i and k

    ab_index = d1[ab];// ab_index corresponds to index of ab in energy tables

    xy_ab = linearize2D(ab_index,xy_index);// xy_ab used to indentify 1x2 loop using row(ab_index) columnar(xy_index) multilplication. 
    energy = internal12[ijkl][xy_ab]; // energy is located in internal12 table in energy_par.c at intersection of ijkl(row) and xy_ab(column)
    energy += Initiation;
    if(DEBUG){
      printf("\n");
      printf("Internal Loop 1x2(J): (%d,%d) %c%c/(%d,%d) %c%c, J: (%d,%d), I: (%d) %f",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),l+1,l+2,i+1,energy);
      printf("\n"); 
    }

  }

  /*
    k l
   b   0
   x   y
    i j
  */

  else{ //if(k-i)==2;  have to switch around to look up in table when bulge of 1 is on j side, ij becomes ji, kl becomes lk
  
    // get i component

    ij_trans_index = GetIndex(trans_seq[j],trans_seq[i]); //ij to ji 
    kl_trans_index = GetIndex(trans_seq[l],trans_seq[k]); //kl to lk

    ij_index = d1[ij_trans_index];
    kl_index = d1[kl_trans_index];

    ijkl = linearize2D(ij_index,kl_index);

    // get j component, only 1 bp missin gon ik side

    x = trans_seq[l+2]; // change from i+1 to l+2
    y = trans_seq[i+1]; // change from l+2 to i+1

    xy = GetIndex(x,y);
    xy_index = d1[xy];
  
    b = trans_seq[i+2];
    ab = GetIndex(0,b); // although in reality this is b,0 it is switched here for the design of the table, the table only supports the model depicted above.

    ab_index = d1[ab];

    xy_ab = linearize2D(ab_index,xy_index);
    energy = internal12[ijkl][xy_ab];
    energy += Initiation;
    if(DEBUG){
      printf("\n");
      printf("Internal Loop 1x2(I): (%d,%d) %c%c/(%d,%d) %c%c, I: (%d,%d), J: (%d) %f",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),i+1,i+2,j-1,energy);
      printf("\n");
    } 
  }

  return energy;
}


// draw picture to represent 2x2
/*

 GetInternal22

 Input:
 1) positon of i in sequence
 2) position of j closing i
 3) position of k in sequence
 4) position of l closing k
 5) integer representation of sequence

 output:
 1) 2x2 internal loop energy corresponding to sequence 

  k l
a     b
x     y
  i j
 5' 3'

 ij and kl are base paired with xy and ab representing the 2x2 loop

*/

double GetInternal22(int i, int j, int k, int l, int * trans_seq){

  int ij_trans_index, kl_trans_index, ij_index, kl_index, ijkl;
  int x,y,a,b, xy_trans_index, ab_trans_index, xy_index, ab_index, xyab;

  double energy;

  // get row component

  ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// corresponds to index of ij in d1
  kl_trans_index = GetIndex(trans_seq[k],trans_seq[l]);// corresponds to index of kl in d1

  ij_index = d1[ij_trans_index];// corresponds to index of ij in energy tables
  kl_index = d1[kl_trans_index];// corresponds to index of kl in energy tables

  ijkl = linearize2D(ij_index,kl_index);// ijkl represents ij as opening bp of internal loop and kl as closing bp of internal loop. Row(ij) columnar(kl) multiplication is used to obtain ijkl which will be row value to look up in internal22
  
  // get column component

  // positions of x,y,a,b in sequence
  x = i+1;
  y = j-1;
  a = i+2;
  b = j-2;
  
  xy_trans_index = GetIndex(trans_seq[x], trans_seq[y]);// corresponds to index of xy in d1
  ab_trans_index = GetIndex(trans_seq[a], trans_seq[b]);// corresponds to index of ab in d1

  xy_index = d1[xy_trans_index];// corresponds to index of xy in energy tables
  ab_index = d1[ab_trans_index];// corresponds to index of ab in energy tables

  xyab = linearize2D(xy_index,ab_index);//uses row(xy) columnar(ab) multiplication to represent 2x2 internal loop

  energy = internal22[ijkl][xyab];// interest of ijkl(row) and column(xyab) in internal22
  energy += Initiation; 
  if(DEBUG){
    printf("\n");
    printf("Internal Loop 2x2: (%d,%d) %c%c/ (%d,%d) %c%c, Loop: (%d,%d) %c%c, (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),i+1,j-1,GetNucl(trans_seq[i+1]),GetNucl(trans_seq[j-1]),i+2,j-2,GetNucl(trans_seq[i+2]),GetNucl(trans_seq[j-2]),energy);
    printf("\n");
  }

  return energy;
}

/* GetInternalEnergy

   Input:
   1) position of i in sequence
   2) position of j closing i
   3) position of k in sequence
   4) position of l closing k
   5) integer representation of sequence 

   Output:
   1) Internal Loop energy * NOTE internal loop of 1x1, 1x2, and 2x2 are not handled here
  

   There are two mismatch penalties to add i+1,j-1 and k-1,l+1
*/ 

double GetILEnergy(int i, int j, int k, int l, int *trans_seq){
  int diffi,diffj,diff;
  diffi = k-i-1;// corresponds to # nucleotides between k and i
  diffj = j-l-1;// corresponds to # nucleotides between j and l
  diff = diffi*diffj;
  if((diffi>2)||(diffj>2))
    return GetInternalEnergy(i,j,k,l,trans_seq);
  switch(diff){
    case 1:
	return GetInternal11(i,j,k,l,trans_seq);
    case 2:
        return GetInternal12(i,j,k,l,trans_seq);
    case 4: 
	return GetInternal22(i,j,k,l,trans_seq);
  /*  default:
	return GetInternalEnergy(i,j,k,l,trans_seq);*/
  }
  return 10000000;
}

double GetILSizeCost(int size){

  if( size <= 30){
    return internal[size];
  }

   return  internal[6]+1.08*log((size)/6.0);
}

double GetInternalEnergy(int i, int j, int k, int l, int * trans_seq){


  int diffi,diffj,asym;
  int x,y, ij_trans_index, xy_trans_index, ij_index, xy_index; 
  double energy = 0;
  energy += Initiation; 
  diffi = k-i-1;// corresponds to # nucleotides between k and i
  diffj = j-l-1;// corresponds to # nucleotides between j and l
  asym = 0;  // part of internal loop penalty 

  if (diffj > diffi){// more unbased pairs on j-l than k-l
   asym = diffj - diffi;// assym is assigned as positive value difference between # unbased paired nucleotides on k-i side vs j-l side
  }
  else if(diffi > diffj){
   asym = diffi - diffj;
  } 

  // internal cost energy
  energy+=GetILSizeCost(diffi+diffj);// cost of internal loop; just the sum of diffi and diffj
  if(DEBUG){ 
    printf("\n");
    printf("Internal Loop cost:(%d) (%d,%d) %c%c and (%d,%d) %c%c %f\n",diffi+diffj,i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),energy);
  }

  // assymetry
  energy+=assymetry*(asym);// assymetry penalty, value is in energy_par.c
  if(DEBUG){
    printf("Assymetry: %d Energy: %f\n", asym, assymetry*(asym));
  }

  //mismatch 1
  x = i+1;
  y = j-1;

  ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// corresponds to ij index in d1
  xy_trans_index = GetIndex(trans_seq[x],trans_seq[y]);// corresponds to kl index in d1

  ij_index = d1[ij_trans_index];// corresponds to index of ij in energy tables
  xy_index = d1[xy_trans_index];// corresponds to index of xy in energy tables

  energy+=TMM[ij_index][xy_index];// corresponds to ij(row) and xy(column) in TMM table in energy_par.c
  if(DEBUG){
    printf("Internal Energy 1st TMM: (%d,%d) %c%c %f\n", x,y,GetNucl(trans_seq[x]), GetNucl(trans_seq[y]), TMM[ij_index][xy_index]);
  }


  //mismatch 2
  x = k-1;
  y = l+1;

  ij_trans_index = GetIndex(trans_seq[l],trans_seq[k]);// corresponds to ij index in d1
  xy_trans_index = GetIndex(trans_seq[y],trans_seq[x]);// corresponds to kl index in d1

  ij_index = d1[ij_trans_index];// corresponds to index of ij in energy tables
  xy_index = d1[xy_trans_index];// corresponds to index of xy in energy tables

  energy+=TMM[ij_index][xy_index];// corresponds to ij(row) and xy(column) in TMM table in energy_par.c; 2nd mismatch
  if(DEBUG){
    printf("Internal Energy 2nd TMM: (%d,%d) %c%c %f\n", x,y,GetNucl(trans_seq[x]), GetNucl(trans_seq[y]), TMM[ij_index][xy_index]);
  }

  // AU, GU closure

  if(trans_seq[j] == 12){// j is U
    if(trans_seq[i] == 1 || trans_seq[i] == 3){// i is either A or G
      energy+=AUGUinternal;// add AU/GU penalty
      if(DEBUG){
        printf("Internal Energy GU/AU closure: %f\n",GUAU_penalty);
      }
    }
  }

  if(trans_seq[l] == 12){// j is U
    if(trans_seq[k] == 1 || trans_seq[k] == 3){// i is either A or G
      energy+=AUGUinternal;// add AU/GU penalty
      if(DEBUG){
        printf("Internal Energy GU/AU opening: %f\n",GUAU_penalty);
      }
    }
  }

    //printf("Internal Loop: Total Energy: %f\n\n",energy);

  return energy;

}

double GetHairpinEnergy(int i, int j, int * trans_seq){

  int ij_trans_index, xy_trans_index, ij_index, xy_index, pos, x, y;
  double energy;
  double cost;   
  energy = 0.0; 
 
  int diff = j-i-1;//diff refers to the size of the loop
  int Cflag = 1; // var used to check for loop of all C's; there is a special penalty in this case; 1 is true for all C's, 0 is false 
  int firstunbp_trans_index;// represents i+1,i+2
  int middleunbp_trans_index;
  int lastunbp_trans_index;// represents i+3 and i+4 if necessary(loop size 4) 
  int firstunbp;// var to look up index in energy tables
  int lastunbp;// var to look up index in energy tables
  int middleunbp;
  int bpfirst2inloop;
  int last4inloop;
  int bp;
  int bp_trans_index;
  int loop; 

  // There are special hairpin costs for certain haiprin sequences of length 3,4 and 6.
  if(diff >= 3){
    
    // cost
    if(diff == 3){// want to represent 3 nucleotides as IJ0K where I is i+1, J is i+2, and K is i+3; this allows lookup in hairpin3 table located in energy_par.c
      
      bp_trans_index = GetIndex(trans_seq[i],trans_seq[j]);      
      firstunbp_trans_index= GetIndex(trans_seq[i+1],trans_seq[i+2]);// corresponds to i+1,i+2 index in d1
      lastunbp_trans_index = GetIndex(0,trans_seq[i+3]);// corresponds to i+3,0 index in d1
      bp = d1[bp_trans_index];
      firstunbp = d1[firstunbp_trans_index];// corresponds to index of i+1,i+2 in energy tables
      lastunbp =  d1[lastunbp_trans_index];// corresponds to index of i+3 in energy tables 
      
      loop = linearize2D(lastunbp,firstunbp); 
      
      cost =  hairpin3[bp][loop];
     
      if (firstunbp == 12 && lastunbp == 3){
        cost += GGGhairpin; // special bonus for GGG hairpin
      } 
      
    
    }else if(diff == 4){ //if(diff==3)
      
      bp_trans_index = GetIndex(trans_seq[i],trans_seq[j]); 
      firstunbp_trans_index= GetIndex(trans_seq[i+1],trans_seq[i+2]);// corresponds to i+1,i+2 index in d1
      lastunbp_trans_index = GetIndex(trans_seq[i+3],trans_seq[i+4]);// corresponds to i+3,i+4 index in d1

      bp = d1[bp_trans_index];
      firstunbp = d1[firstunbp_trans_index];// corresponds to index of i+1,i+2 in energy tables
      lastunbp = d1[lastunbp_trans_index];

      loop = linearize2D(firstunbp,lastunbp);

      cost = hairpin4[bp][loop];
    }
    
    else if(diff == 6){// if(diff==4)
      bp_trans_index = GetIndex(trans_seq[i],trans_seq[j]);
      firstunbp_trans_index= GetIndex(trans_seq[i+1],trans_seq[i+2]);// corresponds to i+1,i+2 index in d1
      middleunbp_trans_index = GetIndex(trans_seq[i+3],trans_seq[i+4]);// corresponds to i+3,i+4 index in d1
      lastunbp_trans_index = GetIndex(trans_seq[i+5],trans_seq[i+6]);

      bp = d1[bp_trans_index];
      firstunbp = d1[firstunbp_trans_index];// corresponds to index of i+1,i+2 in energy tables
      middleunbp = d1[middleunbp_trans_index];
      lastunbp = d1[lastunbp_trans_index];
  
      bpfirst2inloop = linearize2D(bp,firstunbp);
      last4inloop = linearize2D(middleunbp,lastunbp);

      cost = hairpin6[bpfirst2inloop][last4inloop];
    }

    // regular case
    else if(diff <= 30){// if (diff==6)
      cost=hairpin[diff];
    }
   
    else{// if diff > 30 use extrapolation below
      cost=hairpin[6]+1.75*R*T*log((diff)/6.0);
    }
    
    energy+=cost;
    if(DEBUG){
      printf("\n");
      printf("Hairpin Loop: Cost (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),cost);
    }
    if(diff > 3){
      // mismatch located at i+1, j-1, note that i+1,j-1 are not basepairs
      x = i+1;// x is position of i+1
      y = j-1;// y is position of j-1

      ij_trans_index = GetIndex(trans_seq[i],trans_seq[j]);// corresponds to index of i,j in d1
      xy_trans_index = GetIndex(trans_seq[x],trans_seq[y]);// corresponds to index of i+1, j-1 in d1

      ij_index = d1[ij_trans_index];// corresponds to index of i,j in energy table TMM in energy_par.c
      xy_index = d1[xy_trans_index];// corresponds to index of i+1,j-1 in energy table TMM in energy_par.c
      //printf("%d - %d : %d, %d - %d : %d\n", i,j,ij_index,x,y,xy_index);
      energy+=TMM[ij_index][xy_index];
      if(DEBUG){
        printf("Hairpin Loop: TMM (%d,%d) %c%c and (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),i+1,j-1,GetNucl(trans_seq[i+1]),GetNucl(trans_seq[j-1]),TMM[ij_index][xy_index]);
      }
    }
    else{
      if(DEBUG){
        printf("TMM not added because loop is too small\n");
      }
    } 
    // uu or ga mismatch
    if((xy_trans_index == 28 || xy_trans_index == 10) && diff > 3){//special mismatch penalty for uu(index 28 in d1) /ga(index 10 in d1) 
      energy+= UUGAmismatch;
      if(DEBUG){
        printf("Haiprin Loop: UU or GA mismatch: (%d,%d) %c%c  %f\n",i+1,j-1,GetNucl(trans_seq[i+1]),GetNucl(trans_seq[j-1]),UUGAmismatch);
      }
    } 
    
    // gg mismatch
    else if(xy_trans_index == 25 && diff > 3){// gg is at index 25 in d1
      energy+= GGmismatch;
      if(DEBUG){
        printf("Hairpin Loop, GG mismatch: (%d,%d) %c%c %f\n",i+1,j-1,GetNucl(trans_seq[i+1]),GetNucl(trans_seq[j-1]),GGmismatch); 
      }
    } 
   
    // au/gu closure
    if(ij_trans_index == 1 && diff > 3){//gu is at index 1 in d1 and au is at index 3 in d1
      energy+=GUclosure;
      if(DEBUG){
        printf("Hairpin Loop, Special GU closure: (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),GUclosure);
      }
    }

    // check for c loop
    for(pos = i+1; pos < j; pos++){
      if(trans_seq[pos] != 7){
        Cflag = 0;// Cflag becomes false if nucleotide other than C present in loop
        break;
      }
    }
      
    if(Cflag){// if Cflag == 1; so all C's in loop
      if (diff == 3){
        energy+= C3loop;// special penalty for loop size of 3 all C
        if(DEBUG){
          printf("Hairpin Loop, C3 loop: (%d,%d,%d) %c%c%c %f\n",i+1,i+2,i+3,GetNucl(trans_seq[i+1]),GetNucl(trans_seq[i+2]),GetNucl(trans_seq[i+3]),C3loop);
        }
      }
      else{// if(diff == 3)
        energy+=(CloopA*diff+CloopB);// formula given for loop size greater than three, all C a*diff+b where CloopA and CloopB can be found in energy_par.c 
        if(DEBUG){
          printf("Hairpin Loop, C loop: size(%d) %f\n", diff, CloopA*diff+CloopB);
        }
      }
    }
  }

  else{//if(diff >= 3)
    energy+= hairpin[diff];
    if(DEBUG){
      printf("Hairpin Loop Cost: (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),hairpin37[diff]);
    }
  }

  if(DEBUG){
    printf("Hairpin Loop Total Energy: %f\n\n",energy);
  }

  return energy;
}

/* GetMultiloopCoaxEnergy

   Input:
   1) Int array of Helix where index is helix number and value is BP structure that has startin bp(int i, int j) and closing bp (int i, int j)
   2) 2D int array of children where row is helix number and column is int representing children of helix
   3) var to give helix number
   4) bp to represent 
   5) integer representation of sequence

   Output:
   1) energy of multiloop including coaxial stacking

   Function is used when a multiloop is closed in GetStructureEnergy. Function calculates number of branches in ML, number base pairs inbetween branches and calculates penalty. 

*/

double GetMultiloopCoaxEnergy(H* Helix, int** children, int pos, int * bp, int * trans_seq, int CD_flag, int tmm_flag)
{
  double energy,assym;

  int num_branch,index,Helix_index,unbp,coax_flag;// coax_flag is used for when there are more than 2 helices coaxially stacked
  int i,j,k,l;
  double tmp_coax_energy, min_coax_energy; // necessary for coaxial stacking. In situation where ()()() only 2 of helices can be coaxially stacked so MFE coaxial stacking is calculated

  unbp = 0;
  coax_flag = 0;

  num_branch = children[pos][0];// 1st column holds number of branches for this parent
  i = Helix[pos].close.i; // pos is the parent, so this is the parent i
  j = Helix[pos].close.j; // parent j
 
  if(DEBUG){
    printf("Multiloop: (%d,%d) %c%c: Number of branches: %d\n", Helix[pos].close.i, Helix[pos].close.j, GetNucl(trans_seq[Helix[pos].close.i]),GetNucl(trans_seq[Helix[pos].close.j]),num_branch+1);
  }
       
  for(index = 1; index <= num_branch; index++){//start at 1 b/c 0 holds number of branches
    Helix_index = children[pos][index];// Helix_index is the index of the child to look up in Helix
    k = Helix[Helix_index].start.i;// k is the start bp i
    l = Helix[Helix_index].start.j;// l is the start bp j
    
    if (index == 1){
      unbp+= k-i-1;// situation where looking at last bp in parent helix(i,j) and first bp(k,l) in first child helix: (((..(((..)))...(((...)))...)))
                                                                                                                 //     i  k      l               j 

      if(DEBUG){
        printf("Adding unbase paired: i is %d, k is %d, diff is %d\n",i,k,k-i-1);
      }
    }
    else{
      unbp+= k-j-1;//situation as follows: (((..(((..)))...(((...)))...)))
                                          //    i      j   k       l
       if(DEBUG){
         printf("Adding unbase paired: %d %d %d %d  k is %d, j is %d, diff is %d\n",i,j,k,l,k,j,k-j-1);
        }
    }
    
    // coaxial stacking
    if (CD_flag){ 
      if (coax_flag == 1){
        if(k-i-1 <= 1){// coaxially stacked, want to choose minimum
          tmp_coax_energy = GetCoaxEnergy(i,j,k,l,bp,trans_seq);
          if (tmp_coax_energy < min_coax_energy){
          min_coax_energy = tmp_coax_energy;
          }
        }
        else{// not coaxially stacked, so take the lowest coaxially stacked energy
          energy+=min_coax_energy;
          coax_flag = 0;
        }
      }

      // check for coaxial stacking 
      else if(k-i-1 == 1 || k-i-1 == 0 || k-j-1 == 1 || k-j-1 == 0){//if(coax_flag ==1)
          coax_flag = 1;
          min_coax_energy = GetCoaxEnergy(i,j,k,l,bp,trans_seq);
      }
    
     
      else if(bp[l+1]==-1){
        energy+= tmm_flag*GetTMMEnergy(l,k,l+1,k-1,trans_seq);
      }
   
    }

    else if(bp[l+1]==-1 && bp[k-1]==-1){
        energy+= tmm_flag*GetTMMEnergy(l,k,l+1,k-1,trans_seq);
    }

    i = Helix[Helix_index].start.i;
    j = Helix[Helix_index].start.j;

  }// close  for(index = 1; index <= num_branch; index++)
  
  // last branch is i,j and opening branch is k,l: situation where you are at last branch in 
  // loop and much check if coaxially stacked with parent helix   (((..(((..)))...(((...)))...)))
  //                                                                k             i       j   l                              
  k = Helix[pos].close.i;
  l = Helix[pos].close.j;
  unbp+= l-j-1;
  if(DEBUG){
     printf("Adding unbase paired: l is %d, j is %d, diff is %d\n",l,j,l-j-1);
  }

  // coaxial stacking with last branch

  if (CD_flag){
    if (coax_flag == 1){
        if(l-j-1 == 1 || l-j-1 == 0){
          tmp_coax_energy = GetCoaxEnergy(i,j,k,l,bp,trans_seq);
          if (tmp_coax_energy < min_coax_energy){
            energy += tmp_coax_energy;
          }
        }
        else{
          energy+=min_coax_energy;
          coax_flag = 0;
        }
    }

    else{
        if(l-j-l <= 1){
          energy+=GetCoaxEnergy(i,j,k,l,bp,trans_seq);
        }
    } 
  }

  // Multiloop Assymetry Calculation
  assym = unbp/(num_branch+1);
  if (assym > 2){
    assym = 2;
  }
 
  if(DEBUG){
    printf("Multiloop: (%d,%d) %c%c, unbp %d\n", Helix[pos].close.i, Helix[pos].close.j,GetNucl(trans_seq[Helix[pos].close.i]),GetNucl(trans_seq[Helix[pos].close.j]), unbp);
  }

  // a term in multiloop energy calculation
  energy += MultiloopA;
  if(DEBUG){
    printf("Multiloop A value: %f\n",MultiloopA);
  }

  /* average assymetry
  energy += MultiloopAssym*assym;
  if(DEBUG){
    printf("Average Assymetry: %f %f\n",assym,MultiloopAssym*assym);
  }*/
    
  // num_branch 
  energy += MultiloopB*(num_branch+1);
  if(DEBUG){
    printf("Number of Branches: %d %f\n",num_branch+1,MultiloopB*(num_branch+1));
  }

  /* Multiloop Strain
  if (num_branch + 1 == 3 && unbp < 2){
    energy+= MultiloopStrain;
  }*/

  if(DEBUG){
    printf("Multiloop Total Energy: %f\n\n",energy);
  }

  return energy;
}

double GetExteriorLoopEnergy(H* Helix, int** children, int * bp, int * trans_seq, int CD_flag, int tmm_flag)
{
  double energy,assym;

  int index,Helix_index,coax_flag,num_branch;
  int i,j,k,l;
  double tmp_coax_energy, min_coax_energy;

  energy = 0.0;

  coax_flag = 0;

  num_branch = children[0][0];

  if(DEBUG){
    printf("Exterior loop: (%d,%d) %c%c: Number of branches: %d\n", Helix[children[0][1]].start.i, Helix[children[0][num_branch]].start.j, GetNucl(trans_seq[Helix[children[0][1]].start.i]), GetNucl(trans_seq[Helix[children[0][num_branch]].start.j]),num_branch);
  }

  for(index = 1; index < num_branch; index++){
    printf("Energy now is: %f\n",energy);
    Helix_index = children[0][index];
    i = Helix[Helix_index].start.i;
    j = Helix[Helix_index].start.j;
    Helix_index = children[0][index+1]; 
    k = Helix[Helix_index].start.i;
    l = Helix[Helix_index].start.j;

    // coaxial stacking

    if (CD_flag){
      if (coax_flag == 1){
        if(k-j-1 <= 1){// coaxially stacked, want to choose minimum
          tmp_coax_energy = CD_flag*GetCoaxEnergy(i,j,k,l,bp,trans_seq);
          if (tmp_coax_energy < min_coax_energy){
            min_coax_energy = tmp_coax_energy;
          }
        }
        else{// not coaxially stacked, so take the lowest coaxially stacked energy
          energy+=min_coax_energy;
          coax_flag = 0;
        }
      }

      else if(k-j-1 <= 1){// if(coax_flag == 1)
          coax_flag = 1;
          min_coax_energy = CD_flag*GetCoaxEnergy(i,j,k,l,bp,trans_seq);
          if(index = num_branch-1){
            energy+=min_coax_energy;
            coax_flag = 0;
          }
      }
    
      // add T.M.M
      else if (index == 1 && bp[j+1] == -1){// check for TMM on first helix
        energy+= tmm_flag*GetTMMEnergy(j,i,j+1,i-1,trans_seq);
        if(bp[l+1] == -1){// check for TMM on next helix
          energy+= tmm_flag*GetTMMEnergy(l,k,l+1,k-1,trans_seq);
        }  
      }

      else if(bp[l+1] == -1){// check for TMM on next helix
        energy+= tmm_flag*GetTMMEnergy(l,k,l+1,k-1,trans_seq);
      }

    }
    
    else if (index == 1 && bp[j+1] == -1){// check for TMM on first helix
        energy+= tmm_flag*GetTMMEnergy(j,i,j+1,i-1,trans_seq);
        if(bp[l+1] == -1){// check for TMM on next helix
          energy+= tmm_flag*GetTMMEnergy(l,k,l+1,k-1,trans_seq);
        }
    }

    else if(bp[l+1] == -1){// check for TMM on next helix
        energy+= tmm_flag*GetTMMEnergy(l,k,l+1,k-1,trans_seq);
    }

  }// close for(index = 1; index < num_branch; index++)

  if(DEBUG){ 
    if(tmm_flag == 0){
      printf("TMM are disregarded for Exterior Loops\n");
    }
    if(CD_flag == 0){
      printf("Coaxial Stacking are disregarded for Exterior Loops\n");
    }
    printf("Exterior Loop: (%d,%d) %c%c %f\n", Helix[children[0][1]].start.i, Helix[children[0][num_branch]].start.j,GetNucl(trans_seq[Helix[children[0][1]].start.i]),GetNucl(trans_seq[Helix[children[0][num_branch]].start.j]), energy);
  }

  return energy;
}

/* GetCoaxEnergy
 
   Inputs: 
   1) position i in integer sequence; position i should correspond to the left bp in the opening of a helix
   2) position j closes i; right bp in the opening of a helix
   3) position k in integer sequence; correspond to left bp in opening of a helix that stacks on ij
   4) position l closes k
   5) integer representation of sequence

   Output: Coaxial Stacking energy based on parameters

   Note: Must account for following instances:

   Opening helix in ML and first branch: ((((...(((...)))...(((...)))...))))
                                            i   k       l               j

   Adjacent Helices:         ((((...(((...)))...(((...)))...))))
                                    i       j   k       l

   Last branch in ML and opening helix for ML: ((((...(((...)))...(((...)))...))))
                                                  k               i       j   l


   In each case one must compare coaxial stacking with 5' and 3' dangles with individual helices to coaxial stacking of both helices.

*/

double GetCoaxEnergy(int i,int j,int k,int l,int *bp, int * trans_seq){

  // if basepair mediated coaxial stacking is TMM plus discontinuity penalty

  /* if flush coaxial stacking treat as stacking basepairs
     must compare with having just dangles */


  double energy = 0.0;
  double tmp_energy, dangle_energy, coax_energy;

  /* ij is the closing bp of stem that opens multiloop, kl is opening bp of first branch
     one helix followed by another ().()
  */
  if (k-i-1 == 1){// unpaired bp in between

    if(bp[l+1] == -1){// unpaired nuc in l+1 pos
      energy = GetTMMEnergy(l,k,l+1,k-1,trans_seq);
      energy+= DisconPenalty; 
    
      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d, %d: Unpaired nucleotide in l+1 position, energy is %f\n\n",i,j,k,l,energy);
      }
    }

    if(bp[j-1] == -1){// unpaired nuc in j-1 pos
      tmp_energy = GetTMMEnergy(i,j,i+1,j-1,trans_seq);
      if (tmp_energy < energy){
        energy = tmp_energy;
      }

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d, %d: Unpaired nucleotide in j-1 position, energy is %f\n\n",i,j,k,l,energy);
      } 
    }
  }// close if (k-i-1 == 1)
  /*
   ij is the closing bp of stem that opens multiloop, kl is opening bp of first branch
   unmediated stacking ()()
  */
  else if(k-i-1 == 0){ //if(k-i-1 == 1) need to compare coaxial energy with dangle energy
    coax_energy = GetStackEnergy(i,j,k,l,trans_seq);   
    if(bp[j-1] == -1 && bp[l+1] == -1){
      dangle_energy = GetDangleEnergy(i,j,0,j-1,trans_seq);  // 3 prime dangle
      dangle_energy += GetDangleEnergy(0,l+1,k,l,trans_seq); // 3 prime dangle 
    } 
  
    else if(bp[j-1] == -1){//  if(bp[j-1] == -1 && bp[l+1] == -1)
      dangle_energy = GetDangleEnergy(i,j,0,j-1,trans_seq);  // 3 prime dangle
    }
    
    else if(bp[l+1] == -1){// else if(bp[j-1] == -1){
      dangle_energy += GetDangleEnergy(0,l+1,k,l,trans_seq); // 3 prime dangle
    }

    if(coax_energy <= dangle_energy){
      energy = coax_energy;
    
      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d %d: Flush, energy: %f\n\n",i,j,k,l,energy);
      }
    }

    else{// if(coax_energy <= dangle_energy)
      energy = dangle_energy;
      if(DEBUG){
        printf("Calculating dangles between %d, %d and %d %d: energy:  %f\n\n",i,j,k,l,energy);
      }
    }

  }

  /* one helix followed by another ().() */
  else if (k-j-1 == 1){// else if(k-i-1 == 0); unpaired bp in between

    if(bp[l+1] == -1){// unpaired nuc in l+1 pos
      energy = GetTMMEnergy(l,k,l+1,k-1,trans_seq);
      energy+= DisconPenalty;

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d, %d: Unpaired nucleotide in l+1 position, energy is %f\n\n",i,j,k,l,energy);
      }
    }

    if(bp[i-1] == -1){// unpaired nuc in i-1 pos
      tmp_energy = GetTMMEnergy(j,i,j+1,i-1,trans_seq);
      if (tmp_energy < energy){
        energy = tmp_energy;
      }

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d, %d: Unpaired nucleotide in i-1 position, energy is %f\n\n",i,j,k,l,energy);
      }
    }
  }
  
  /* unmediated coaxial stacking ()() */
  else if(k-j-1 == 0){ //else if (k-j-1 == 1); need to compare coaxial energy with dangle energy
    coax_energy = GetStackEnergy(i,j,k,l,trans_seq);
    if(bp[j+1] == -1 && bp[l+1] == -1){
      dangle_energy = GetDangleEnergy(i-1,0,i,j,trans_seq);  // 5 prime dangle
      dangle_energy += GetDangleEnergy(k,l,0,l+1,trans_seq); // 3 prime dangle 
    }

    else if(bp[j-1] == -1){
      dangle_energy = GetDangleEnergy(i-1,0,i,j,trans_seq);  // 5 prime dangle
    }

    else if(bp[l+1] == -1){
      dangle_energy += GetDangleEnergy(k,l,0,l+1,trans_seq); // 3 prime dangle
    }

    if(coax_energy <= dangle_energy){
      energy = coax_energy;

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d %d: Flush, energy: %f\n\n",i,j,k,l,energy);
      }
    }

    else{//if(coax_energy <= dangle_energy)
      energy = dangle_energy;
      if(DEBUG){
        printf("Calculating dangles between %d, %d and %d %d: energy:  %f\n\n",i,j,k,l,energy);
      }
    }

  }
 
  /* last branch of multiloop coaxial stacks with opening stem of ML */
  else if (l-j-1 == 1){// else if (k-j-1) == 1; unpaired bp in between

    if(bp[k+1] == -1){// unpaired nuc in k+1 pos
      energy = GetTMMEnergy(k,l,k+1,l-1,trans_seq);
      energy+= DisconPenalty;

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d, %d: Unpaired nucleotide in k+1, energy is %f\n\n",i,j,k,l,energy);
      }
    }

    if(bp[i-1] == -1){// unpaired nuc in i-1 pos
      tmp_energy = GetTMMEnergy(j,i,j+1,i-1,trans_seq);
      if (tmp_energy < energy){
        energy = tmp_energy;
      }

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d, %d: Unpaired nucleotide in i-1 position, energy is %f\n\n",i,j,k,l,energy);
      }
    }
  }
  
  else if(l-j-1 == 0){ // need to compare coaxial energy with dangle energy
    coax_energy = GetStackEnergy(i,j,k,l,trans_seq);
    if(bp[i-1] == -1 && bp[k+1] == -1){
      dangle_energy = GetDangleEnergy(k,l,k+1,0,trans_seq);
      dangle_energy += GetDangleEnergy(i-1,0,i,j,trans_seq);
    }

    else if(bp[j-1] == -1){
      dangle_energy = GetDangleEnergy(k,l,k+1,0,trans_seq);
    }

    else if(bp[l+1] == -1){
      dangle_energy += GetDangleEnergy(i-1,0,i,j,trans_seq);
    }

    if(coax_energy <= dangle_energy){
      energy = coax_energy;

      if(DEBUG){
        printf("Calculating Coaxial Stacking between %d, %d and %d %d: Flush, energy: %f\n\n",i,j,k,l,energy);
      }
    }

    else{
      energy = dangle_energy;
      if(DEBUG){
        printf("Calculating dangles between %d, %d and %d %d: energy:  %f\n\n",i,j,k,l,energy);
      }
    }

  }

  if(DEBUG){
    printf("Coaxial Stacking: (%d,%d) %c%c/ (%d,%d) %c%c %f\n",i,j,GetNucl(trans_seq[i]),GetNucl(trans_seq[j]),k,l,GetNucl(trans_seq[k]),GetNucl(trans_seq[l]),energy);
  }

  return energy;
}

double GetStructureEnergy(int * trans_seq, char * secstr, int triplet_flag, int CD_flag, int tmm_flag){
  double energy;
  energy = 0.0; 
  int pos,i,j,k,l,k_index,index,len, par_index, par_j,par_ij,child_index,diffi,diffj,flag;
  int * bp;// ordered bp list
  int ** BooleanMatrix;// need to generate Boolean Matrix to check bp in between i,j
  int trip_pos;
  //Added for triple model
  int u;//i+2
  int v;//j-2
  int diffi2;//diff between u and k
  int diffj2;//diff betwenn v and l
  //----------------------

  flag = 0;  
  len = strlen(secstr);
  bp = GetBPList(secstr,len);

  BooleanMatrix = Allocate2Dmatrix(len,len);
  BooleanMatrix = BPBetweenBool(BooleanMatrix,bp,len);
  index = 1;

  // helix, stores start and closing BP
  H * Helix;  
  Helix = (H*)(malloc(MAXBP * sizeof(H)));
   
  // parent

  int par[MAXBP];
  
  // children
  
  int **children = Allocate2Dmatrix(MAXBP, MAXCHILDREN);
  for (pos = 0; pos < MAXBP; pos++){
    children[pos][0] = 0;// want to start 1 over from 0, 0 is where you store index
  }

  int numbp = 0;

  // initialize helix
  for(pos = 0; pos < len; pos++){
    if(bp[pos] != -1){
      i = pos;
      j = bp[pos];
      numbp++;
      break;
    }
  }
  Helix[0].start.i = 0;
  Helix[0].start.j = len;
  Helix[index].start.i = i;
  Helix[index].start.j = j;


  // initialize parent and child
  par[1] = 0;
  int lasti = i;
  int lastj = j;
  //iterate through bps

  for(pos = i+1; pos < len; pos++){
    //find i,j,k,l
   
    i = lasti;
    j = lastj;
 
    if(bp[pos] != -1 && bp[pos] > pos){
      k = pos;
      l = bp[pos];
      numbp++;
   
      lasti = k;
      lastj = l;

      for(trip_pos = k+1; trip_pos < len; trip_pos++){
        if (bp[trip_pos] != -1 && bp[trip_pos] > trip_pos){
          u = trip_pos;
          v = bp[trip_pos];
          break;
        }
      }
      // get space inbetween i and k and l and j, for loop energy calculation
      diffi = k-i;
      diffj = j-l;

      // check if i,j is closing bp of stem so no bases in between
      if(BooleanMatrix[i][j] == 0){
 
          // close helix
          Helix[index].close.i = i;
          Helix[index].close.j = j; 
    
          // open a new helix
          Helix[index+1].start.i = k;
          Helix[index+1].start.j = l;
    
          // get parent helix, need to get j of start of parent helix. 
          par_index = par[index];
          par_j = Helix[par_index].start.j;

          //parent of k,l is parent of parent of i,j - closing a multiloop
          if(k > par_j){
            if(DEBUG){
              printf("closing last branch of multiloop with %d %d %d %d\n",i,j,k,par_j);
            }
            
            // move up a level to get parent
            par_ij = par[index];
            par[index+1] = par[par_ij];
            energy+= GetHairpinEnergy(i,j,trans_seq);
            if(DEBUG){
              printf("ENERGY: %f\n", energy); 
            }

            //add child to parent
            children[par[index+1]][0]++;//increase number of branches by 1;  
            child_index = children[par[index+1]][0];
            children[par[index+1]][child_index] = index+1;
            if(DEBUG){
              printf(" Adding child %d to parent %d\n",index+1,par[index+1]);
            }
            energy+= GetMultiloopCoaxEnergy(Helix, children, par[index], bp, trans_seq, CD_flag, tmm_flag);
            index++;
          }
  
          else{
  
            // parent doesn't change
            par[index+1] = par[index];
            energy+= GetHairpinEnergy(i,j,trans_seq);
            printf("ENERGY: %f\n", energy); 

            //add child to parent 
            children[par[index+1]][0]++;// increase number of branches by 1 
            child_index = children[par[index+1]][0];
            children[par[index+1]][child_index] = index+1;
            printf("In ML: Adding child %d to parent %d\n",index+1,par[index+1]);
            index++;
          } 
      
      }

      // opening up a multiloop
      else if(BooleanMatrix[l][j] == 1){
      
          Helix[index].close.i = i;
          Helix[index].close.j = j;

          Helix[index+1].start.i = k;
          Helix[index+1].start.j = l;

          par[index+1] = index;
   

          children[index][0]++;// increase number of branches by 1   
          child_index = children[index][0];
          children[index][child_index] = index+1;
          printf("Opening ML: Adding child %d to parent %d\n",index+1,index);
    
          index++;
      }  
 
      else {
        
 
          // stacking and triplet energy 
          if(diffi == 1 && diffj == 1){
            if(pos<len-1){
              diffi2 = u-k;
              diffj2 = l-v;
              
              if((diffi2 + diffj2)>= 2 && (diffi2 +diffj2)<=3 && triplet_flag){
                //triplet_flag = triplet_flag + 1;
                energy+=GetTripletEnergy(i,j,k,l,u,v,trans_seq);
                //triplet_flag = 0;
              }
              else{
                energy+=GetStackEnergy(i,j,k,l,trans_seq);
              }
            }else{            
              energy+=GetStackEnergy(i,j,k,l,trans_seq);
            }
            printf("ENERGY: %f\n", energy);
          }
   
          else{// if(diffi == 1 && diffj == 1)
            // bulge
            if(diffi == 1 && diffj > 1){
              energy+= GetBulgeEnergy(i, j, k, l,u,v,trans_seq, triplet_flag);
              printf("ENERGY: %f\n", energy);
            }
            
            else if(diffi > 1 && diffj == 1){
              energy+= GetBulgeEnergy(i, j, k, l,u,v,trans_seq, triplet_flag);
              printf("ENERGY: %f\n", energy);
            }
            // internal 1x1
            else if(diffi == 2 && diffj == 2){
              energy+= GetInternal11(i,j,k,l,trans_seq);
              printf("ENERGY: %f\n", energy);
            }
            // internal 1x2
            else if(diffi == 2 && diffj == 3){
              energy+= GetInternal12(i,j,k,l,trans_seq);
              printf("ENERGY: %f\n", energy);
            }

            else if(diffi == 3 && diffj == 2){
              energy+= GetInternal12(i,j,k,l,trans_seq);
              printf("ENERGY: %f\n", energy);
            }
            //internal 2x2 
            else if(diffi == 3 && diffj == 3){
              energy+= GetInternal22(i,j,k,l,trans_seq);
              printf("ENERGY: %f\n", energy);
            }
            // internal energy
            else{// if(diffi == 1 && diffj > 1)
              energy+=GetInternalEnergy(i,j,k,l,trans_seq);
              printf("ENERGY: %f\n", energy);
            }

          }  
      }
    }
    else if(pos == len-1 && numbp == 1){
      energy+= GetHairpinEnergy(i,j,trans_seq);
      printf("ENERGY: %f\n", energy);
    }
    else if(pos == len-1 && numbp == 0){
      printf("Empty Structure: Energy is 0.0\n");
    }
    else{
      continue;
    }
  }


  // special case for the closing bp of very last branch, not accounted for in main loop

  if(BooleanMatrix[i][j] == 0){
    // close helix
    Helix[index].close.i = i;
    Helix[index].close.j = j;
    // parent doesn't change
    energy+= GetHairpinEnergy(i,j,trans_seq);
    if(par[index] != 0){
      printf("At Last branch and its in ML\n");
      energy+= GetMultiloopCoaxEnergy(Helix, children, par[index], bp, trans_seq,CD_flag,tmm_flag);      
    }
    else{
    //add child to parent 
    children[par[index]][0]++;// increase number of branches by 1 
    child_index = children[par[index]][0];
    children[par[index]][child_index] = index;

    printf("Adding child %d to parent %d\n",index,par[index]);
    }
  } 

  energy+=GetExteriorLoopEnergy(Helix, children, bp, trans_seq, CD_flag,tmm_flag);//exterior loop
 
  printf("FINAL ENERGY: %f\n", energy); 
  
  //free

  free(Helix);
  for(i = 0; i < len; i++){
    free(BooleanMatrix[i]);
  }
  free(BooleanMatrix);
  for(i = 0; i < MAXBP; i++){
    free(children[i]);
  }
  free(children);
  free(bp);
  return energy;
}

void CreateEnergyTable(double T, int tcflag){

 int i,j;
 double  T37 = 310.15;
 T+=273.15;
 //printf("T37 is %f and T is %f\n",T37,T);
 if(T != T37){
   for(i = 0; i < 16; i++){
     for(j = 0; j < 16; j++){
       if (stack[i][j] != INF){
         if(tcflag == 0){
            stack[i][j] = stack_enthalpy37[i][j] - T*((stack37[i][j]-stack_enthalpy37[i][j])/(-1 *T37));
         } 
         else{ 
            stack[i][j] = stack37[i][j] * (T/T37) - stack_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       }
       else if (TMM[i][j] != INF){
         if (tcflag == 0){ 
           TMM[i][j] = TMM_enthalpy37[i][j] - T*((TMM37[i][j]-TMM_enthalpy37[i][j])/(-1 *T37));
         }  
         else{
           TMM[i][j] = TMM37[i][j] * (T/T37) - TMM_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       }
     }     
     for(j = 0; j < 8; j++){
       if (dangle[i][j] != INF)
         if (tcflag == 0){
           dangle[i][j] = dangle_enthalpy37[i][j] -T*((dangle37[i][j]-dangle_enthalpy37[i][j])/(-1 *T37));
         }
         else{
           dangle[i][j] = dangle37[i][j] * (T/T37) - dangle_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }  
     }

     for(j = 0; j < 64; j++){
       if (hairpin3[i][j] != INF){
         if(tcflag == 0){
           hairpin3[i][j] = hairpin3_enthalpy37[i][j]- T*((hairpin3_37[i][j]-hairpin3_enthalpy37[i][j])/(-1 *T37));
         }
         else{
           hairpin3[i][j] = hairpin3_37[i][j] * (T/T37) - hairpin3_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       } 
     }
   
     for(j = 0; j < 256; j++){
       if (hairpin4[i][j] != INF){
         if(tcflag == 0){
           hairpin4[i][j] = hairpin4_enthalpy37[i][j]- T*((hairpin4_37[i][j]-hairpin4_enthalpy37[i][j])/(-1 *T37));
         }
         else{
           hairpin4[i][j] = hairpin4_37[i][j] * (T/T37) - hairpin4_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       }
     }
   }

   for(i = 0; i < 256; i++){
     for(j = 0; j < 16; j++){
       if (internal11[i][j] != INF){
         if(tcflag == 0){
           internal11[i][j] = internal11_enthalpy37[i][j] - T*((internal11_37[i][j]-internal11_enthalpy37[i][j])/(-1 *T37)); 
         }
         else{
           internal11[i][j] = internal11_37[i][j] * (T/T37) - internal11_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       }
       if (triplet[i][j] != INF){
         if(tcflag == 0){
           triplet[i][j] = triplet_enthalpy37[i][j] - T*((triplet37[i][j]-triplet_enthalpy37[i][j])/(-1 *T37));
         }
         else {
          triplet[i][j] = triplet37[i][j] * (T/T37) - triplet_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       }
     }
     for(j = 0; j < 64; j++){
      if (internal12[i][j] != INF){
        if (tcflag == 0){
          internal12[i][j] = internal12_enthalpy37[i][j] - T*((internal12_37[i][j]-internal12_enthalpy37[i][j])/(-1 *T37)) ;
        }
        else {
          internal12[i][j] = internal12_37[i][j] * (T/T37) - internal12_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
        }
      }
     }
     for(j = 0; j < 256; j++){
       if (internal22[i][j] != INF){
         if(tcflag == 0){
           internal22[i][j] = internal22_enthalpy37[i][j] - T*((internal22_37[i][j]-internal22_enthalpy37[i][j])/(-1 *T37)); 
         } 
         else{
           internal22[i][j] = internal22_37[i][j] * (T/T37) - internal22_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       }
       if (hairpin6[i][j] != INF){
         if(tcflag == 0){
           hairpin6[i][j] = hairpin6_enthalpy37[i][j] - T*((hairpin6_37[i][j]-hairpin6_enthalpy37[i][j])/(-1 *T37));
         } 
         else{
           hairpin6[i][j] = hairpin6_37[i][j] * (T/T37) - hairpin6_enthalpy37[i][j] * T * ((T-T37)/(T37*T37));
         }
       } 
     }
   }

   for(i = 0; i < 31; i++){
     if (hairpin[i] != INF){
       if(tcflag == 0){
         hairpin[i] = hairpin_enthalpy37[i] - T*((hairpin37[i]-hairpin_enthalpy37[i])/(-1 *T37));
       }
       else{
         hairpin[i] = hairpin37[i] * (T/T37) - hairpin_enthalpy37[i] * T * ((T-T37)/(T37*T37));
       }
     }
     if (internal[i] != INF){
       if(tcflag == 0){
         internal[i] = internal_enthalpy37[i] - T*((internal37[i]-internal_enthalpy37[i])/(-1 *T37));
       }
       else{
         internal[i] = internal37[i] * (T/T37) - internal_enthalpy37[i] * T * ((T-T37)/(T37*T37));
       }
     }
     if (bulge[i] != INF){
       if(tcflag == 0){
         bulge[i] = bulge_enthalpy37[i] - T*((bulge37[i]-bulge_enthalpy37[i])/(-1 *T37));
       }
       else{
         bulge[i] = bulge37[i] * (T/T37) - bulge_enthalpy37[i] * T * ((T-T37)/(T37*T37));
       }
     }
   }
     
if (tcflag == 0){
 GUAU_penalty = GUAU_penalty_enthalpy37 -T*((GUAU_penalty37-GUAU_penalty_enthalpy37)/(-1 *T37));
 GGmismatch = GGmismatch_enthalpy37 -T*((GGmismatch37-GGmismatch_enthalpy37)/(-1 *T37)); 
 UUGAmismatch = UUGAmismatch_enthalpy37 -T*((UUGAmismatch37-UUGAmismatch_enthalpy37)/(-1 *T37)); 
 GGGhairpin = GGGhairpin_enthalpy37 -T*((GGGhairpin37-GGGhairpin_enthalpy37)/(-1*T37));
 AUGUinternal = AUGUinternal_enthalpy37 -T*((AUGUinternal37-AUGUinternal_enthalpy37)/(-1 *T37)); 
 GUclosure = GUclosure_enthalpy37 - T*((GUclosure37-GUclosure_enthalpy37)/(-1*T37)); 
 C_bulge = C_bulge_enthalpy37 - T*((C_bulge37 - C_bulge_enthalpy37)/(-1 *T37)); 
 C3loop = C3loop_enthalpy37 -T*((C3loop37-C3loop_enthalpy37)/(-1 *T37));
 CloopA = CloopA_enthalpy37 -T*((CloopA37- CloopA_enthalpy37)/(-1 *T37));
 CloopB = CloopB_enthalpy37 -T*((CloopB37-CloopB_enthalpy37)/(-1 *T37));
 MultiloopStrain = MultiloopStrain_enthalpy37 -T*((MultiloopStrain37-MultiloopStrain_enthalpy37)/(-1 *T37));
 MultiloopA = MultiloopA_enthalpy37 -T*((MultiloopA37-MultiloopA_enthalpy37)/(-1 *T37));
 MultiloopB = MultiloopB_enthalpy37 -T*((MultiloopB37-MultiloopB_enthalpy37)/(-1 *T37));
 MultiloopC = MultiloopC_enthalpy37-T*((MultiloopC37-MultiloopC_enthalpy37)/(-1 *T37));
 MultiloopAssym = MultiloopAssym_enthalpy37-T*((MultiloopAssym37-MultiloopAssym_enthalpy37)/(-1 *T37));
 DisconPenalty = DisconPenalty_enthalpy37 -T*((DisconPenalty37-DisconPenalty_enthalpy37)/(-1 *T37));
 Initiation = 0.0;
}
 //Initiation = Initiation37 * (T/T37) - Initiation_enthalpy37 * T * ((T-T37)/(T37*T37));
  else{
    Initiation = 0.0;
    GUAU_penalty = GUAU_penalty37 * (T/T37) - GUAU_penalty_enthalpy37 * T * ((T-T37)/(T37*T37));
    GGmismatch = GGmismatch37 * (T/T37) - GGmismatch_enthalpy37 * T * ((T-T37)/(T37*T37));
    UUGAmismatch = UUGAmismatch37 * (T/T37) - UUGAmismatch_enthalpy37 * T * ((T-T37)/(T37*T37));
    GGGhairpin = GGGhairpin37 * (T/T37) - GGGhairpin_enthalpy37 * T * ((T-T37)/(T37*T37)); 
    AUGUinternal = AUGUinternal37 * (T/T37) - AUGUinternal_enthalpy37 * T * ((T-T37)/(T37*T37));
    GUclosure = GUclosure37 * (T/T37) - GUclosure_enthalpy37 * T * ((T-T37)/(T37*T37));
    C_bulge = C_bulge37 * (T/T37) - C_bulge_enthalpy37 * T * ((T-T37)/(T37*T37));
    C3loop = C3loop37 * (T/T37) - C3loop_enthalpy37 * T * ((T-T37)/(T37*T37));
    CloopA = CloopA37 * (T/T37) - CloopA_enthalpy37 * T * ((T-T37)/(T37*T37));
    CloopB = CloopB37 * (T/T37) - CloopB_enthalpy37 * T * ((T-T37)/(T37*T37));
   // printf("A %f B %f\n",CloopA,CloopB);
    MultiloopStrain = MultiloopStrain37 * (T/T37) - MultiloopStrain_enthalpy37 * T * ((T-T37)/(T37*T37));
    MultiloopA = MultiloopA37 * (T/T37) - MultiloopA_enthalpy37 * T * ((T-T37)/(T37*T37));
    MultiloopB = MultiloopB37 * (T/T37) - MultiloopB_enthalpy37 * T * ((T-T37)/(T37*T37));
    MultiloopC = MultiloopC37 * (T/T37) - MultiloopC_enthalpy37 * T * ((T-T37)/(T37*T37));
    MultiloopAssym = MultiloopAssym37 * (T/T37) - MultiloopAssym_enthalpy37 * T * ((T-T37)/(T37*T37));
    DisconPenalty = DisconPenalty37 * (T/T37) - DisconPenalty_enthalpy37 * T * ((T-T37)/(T37*T37));
  }
 }
}
