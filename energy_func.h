#include "misc.h"
#include "DataStruct.h"
/* declaration of functions for initializaing, transformation */

double GetTripletEnergy(int i, int j, int k, int l, int u, int v, int* trans_seq); // i, j, k, l, u, v are the values within the sequence to integer, transformation array
double GetILEnergy(int i, int j, int k, int l, int *trans_seq);
double GetStackEnergy(int i, int j, int k, int l, int* trans_seq); // i, j, k, l are the values within the sequence to integer, transformation array
double GetDangleEnergy(int i, int j, int k, int l, int* trans_seq); // i, j, k, l are the values within the sequence to integer, transformation array 
double GetTMMEnergy(int i, int j, int k, int l, int * trans_seq); // i, j, k, l are the values within the sequence to integer, transformation array 
double GetStructureEnergy(int * trans_seq, char * ss, int triplet_energy, int CD_flag, int tmm_flag); // seq and ss  
double GetBulgeEnergy(int i, int j, int k, int l, int u, int v, int* trans_seq, int triplet_flag);// i,j,k,l are indices 
double GetLeftBulgeEnergy(int i, int j, int k, int * trans_seq);
double GetRightBulgeEnergy(int i, int j, int l, int * trans_seq);
double GetInternal11(int i, int j, int k, int l, int * trans_seq);
double GetInternal12(int i, int j, int k, int l, int * trans_seq);
double GetInternal22(int i, int j, int k, int l, int * trans_seq);
double GetInternalEnergy(int i, int j, int k, int l, int * trans_seq);
double GetHairpinEnergy(int i, int j, int * trans_seq);
double GetMultiloopCoaxEnergy(H* Helix, int** children, int pos, int * bp, int * trans_seq, int CD_flag, int tmm_flag);
double GetExteriorLoopEnergy(H* Helix, int** children, int * bp, int * trans_seq, int CD_flag, int tmm_flag);
double GetCoaxEnergy(int i,int j,int k,int l,int *bp, int * trans_seq);//i,j,k,l are indices
void CreateEnergyTable(double T, int tcflag);
double GetILSizeCost(int size);
