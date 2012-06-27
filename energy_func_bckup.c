#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"energy_par.h"
#include"energy_func.h"
#define DEBUG 1

double GetStackEnergy(int i, int j, int k, int l)
{
  int ij_trans_index, kl_trans_index, ij_stack_index, kl_stack_index;
  double energy;

  ij_trans_index = GetIndex(i,j);
  kl_trans_index = GetIndex(k,l);
 
  ij_stack_index = d1[ij_trans_index];
  kl_stack_index = d1[kl_trans_index];

  printf("%d, %d, %d, %d    %d %d %d %d\n",i,j,k,l, ij_trans_index, kl_trans_index, ij_stack_index, kl_stack_index);
  energy = stack37[ij_stack_index][kl_stack_index];
  return energy;  
}

double GetDangleEnergy(int i, int j, int k, int l)
{
  int ij_trans_index, kl_trans_index, ij_dangle_index, kl_dangle_index;
  double energy;

  if (i == 0 || j == 0){
    return GetDangleEnergy(l,k,j,i);
  }  
  else{
    ij_trans_index = GetIndex(i,j);
    kl_trans_index = GetIndex(k,l);

    ij_dangle_index = d1[ij_trans_index];
    kl_dangle_index = d1[kl_trans_index];

    printf("%d, %d, %d, %d    %d %d   %d %d\n",i,j,k,l, ij_trans_index, kl_trans_index, ij_dangle_index, kl_dangle_index);
    energy = dangle37[ij_dangle_index][kl_dangle_index];
    return energy;
  }
}

double GetTMMEnergy(int i, int j, int k, int l)
{
  int ij_trans_index, kl_trans_index, ij_tmm_index, kl_tmm_index;
  double energy;

  ij_trans_index = GetIndex(i,j);
  kl_trans_index = GetIndex(k,l);

  ij_tmm_index = d1[ij_trans_index];
  kl_tmm_index = d1[kl_trans_index];

  printf("%d, %d, %d, %d    %d %d   %d %d\n",i,j,k,l, ij_trans_index, kl_trans_index, ij_tmm_index, kl_tmm_index);
  energy = TMM[ij_tmm_index][kl_tmm_index];
  return energy;
}

double GetEnergy(BP* bp, int len, int numbp ){

  double energy;
  int pos, i,j,k,l,diffi,diffj;

  
  // need to initailize structures to hold current parent, list of Parents(start of multiloop) and children(start of branches), and # of branches per multiloop
  BP * parent = (BP *)malloc(sizeof(BP));
  N * tree = (N *)malloc(len*sizeof(N));
  int ML_Branch[len];

  // initialize parent to (0,0)
  (*parent).i = 0;
  (*parent).j = 0;
 
  // initialize energy to 0.0
  energy = 0.0;

  for(pos = 0; pos < numbp; pos++){
    // looking at i,j and k,l
    i = bp[pos].i;
    j = bp[pos].j;
    k = bp[pos+1].i;
    l = bp[pos+1].j;

    diffi = k-i;
    diffj = j-l;
    
    if (BPbetweenBool(i,j,bp,numbp) == 0){
      if(DEBUG){
        printf("HP loop size %d from %d, %d\n",j-i-1,i,j);
        printf("TMM at %d +1, %d -1\n",i,j);
        printf("\n");
        }
     
      // alloc size for parent
      tree[i].parent = (BP*)malloc(sizeof(BP));
      tree[i].parent = parent;
      
      // closing multiloop
      if (k>(*parent).j){
        printf("Closing of Last Branch of multiloop of %d %d with %d %d\n",(*parent).i,(*parent).j,i,j);
        printf("\n");
        parent = tree[(*parent).i].parent; // go up a level
        
        tree[bp[pos+1].i].parent = (BP*)malloc(sizeof(BP));
        tree[bp[pos+1].i].parent = parent; // assign parent of first bp outside of multiloop to be on same level as start of multiloop
      }

      // within multiloop
      else if(ParentDictBool(parent,tree,len)){
        ML_Branch[(*parent).i]++;
        
        tree[(*parent).i].children = realloc(tree[(*parent).i].children,ML_Branch[(*parent).i]*sizeof(BP)); // realloc # of BP that children can hold to add 1 more branch
        tree[(*parent).i].children[ML_Branch[(*parent).i]-1] = bp[pos+1];//assign start of next helix to be child of start of multiloop

        tree[bp[pos+1].i].parent = (BP*)malloc(sizeof(BP));
        tree[bp[pos+1].i].parent = parent;// assign parent of next helix to be start of multiloop
        
        printf("Added branch starting with %d %d to multiloop %d %d\n",k,l,(*parent).i,(*parent).j); 
        printf("\n");
      } 

      else if(BPbetweenBool(l,j,bp,numbp) == 1){
        printf("Opening multiloop with %d %d\n",i,j);
        (*parent) = bp[pos];
        tree[bp[pos+1].i].parent = (BP*)malloc(sizeof(BP));
        tree[bp[pos+1].i].parent = parent;// assign parent of next bp to be current bp
        
        printf("Added %d %d to multiloop start list\n",i,j);
        printf("\n");
        ML_Branch[(*parent).i] = 1;
        tree[(*parent).i].children = (BP*)malloc(sizeof(BP));
        tree[(*parent).i].children[0] = bp[pos+1];// add opening branch to multiloop as child to parent   	 
      }
    }
    else{
       
        tree[bp[pos+1].i].parent = (BP*)malloc(sizeof(BP));
        tree[bp[pos+1].i].parent = tree[bp[pos].i].parent; // no bifurcation so everything on same level

        printf("%d %d Parent is: %d %d\n",i,j,(*tree[bp[pos].i].parent).i,(*tree[bp[pos].i].parent).j);
        if(diffi == 1 && diffj == 1){
          printf("Stacking base pair from %d %d to %d %d\n",i,j,k,l);
        }
        else if (diffi == 1 && diffj > 1){
          printf("Bulge on j side from %d %d to %d %d\n",i,j,k,l);
        }
        else if (diffi > 1 && diffj == 1){
          printf("Bulge on i side from %d %d to %d %d\n ",i,j,k,l);
        }
        else if (diffi == 2 && diffj == 2){
          printf("1x1 internal loop from %d %d to %d %d\n", i,j,k,l);
        }
        else if (diffi == 3 && diffj == 2){
          printf("1x2 internal loop from on i side from %d %d to %d %d\n", i,j,k,l);
        }
        else if (diffi == 2 && diffj == 3){
          printf("1x2 internal loop from on j side from %d %d to %d %d\n", i,j,k,l);
        }
        else if (diffi == 3 && diffj == 3){
          printf("2x2 internal looop from %d %d to %d %d\n",i,j,k,l);
        }
        else{
          printf("internal loop from %d %d to %d %d\n",i,j,k,l);
        }
        printf("\n");
    }     
  }
  return energy;
}
