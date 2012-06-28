#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "energy_func.h"

/* Vinodh Mechery
   Implementation of Partition Function calculation
   with coaxial stacking and dangles using energy
   parameters from Mattews 2004
*/

#define DEBUG 1

int main(int argc, char *argv[]){

  char sequence[MAXSIZE];
  char secstr[MAXSIZE];
  int i;
  int k;
  int * integer_seq; 
  int len;
  int numbp; 

  sequence[0] = '0';
  secstr[0] = '0';

  if (argc==4){
    i=1;
    
    // take in inputs
    while (i<argc){
      if(!strcmp(argv[i],"-s")){
        len = strlen(argv[i+1]);
        for(k = 0; k < len; k++){
          sequence[k+1] = argv[i+1][k];
          secstr[k+1] = argv[i+2][k];
        }
        //strncpy(sequence,argv[i+1],strlen(argv[i+1]));
        sequence[strlen(argv[i+1])+1] = '\0';//end sequence  
        //strncpy(secstr,argv[i+2],strlen(argv[i+2]));
        secstr[strlen(argv[i+2])+1] = '\0';
        i+=3;
        if (DEBUG){
          printf("input string: %s");
          PrintStr(sequence);
          printf("\n");
        }
      } else {
        printf("Error: Usage: %s -s SEQUENCE SECSTR\n",argv[0]);
        exit(1);
      }
    }

    integer_seq = TransformSeq(sequence);
    if (DEBUG){
      printf("sequence transformation to integer: \n"); 
      PrintInt(integer_seq,strlen(secstr));    
      PrintStr(secstr);
      printf("\n");
    }
  
    printf("%f\n",GetStructureEnergy(integer_seq,secstr)); 
    free(integer_seq);
 
  } else{
    printf("Error: Usage: %s -s SEQUENCE SECSTR\n",argv[0]);
    exit(1);
  }
  return 0;
}  
