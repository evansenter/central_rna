#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "energy_func.h"
#include "pfunc.h"
#include "mfe.h"
#include "tpfunc.h"
#include "tmfe.h"
//#include "tsampling.h"
//#include "tcalculateProbs1.h"
#include "misc.h"
#include "energy_par.h"
#include "enthalpy_par.h"
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

/* Ivan Dotu
   Implementation of Partition Function calculation
   with triplet model using energy
   parameters from Mattews 2004
*/


//------------Time functions
static struct rusage res;
static double virtual_time;
/*
 *  Start the stopwatch.
 */
void start_timer()
{
  //timed_out = 0;
  getrusage( RUSAGE_SELF, &res );
  virtual_time = (double) res.ru_utime.tv_sec +
    (double) res.ru_stime.tv_sec +
    (double) res.ru_utime.tv_usec / 1000000.0 +
    (double) res.ru_stime.tv_usec / 1000000.0;
  //if (initial_virtual_time == 0) initial_virtual_time = virtual_time;
}

/*
 *  Stop the stopwatch and return the CPU time used in seconds.
 */
double elapsed_time()
{
  getrusage( RUSAGE_SELF, &res );
  return( (double) res.ru_utime.tv_sec +
          (double) res.ru_stime.tv_sec +
          (double) res.ru_utime.tv_usec / 1000000.0 +
          (double) res.ru_stime.tv_usec / 1000000.0
          - virtual_time );
}


int main(int argc, char *argv[]){

  char sequence[MAXSIZE];
  char secstr[MAXSIZE];
  int i;
  int k;
  int * integer_seq; 
  int len;
  int numbp; 
  int mode;// 0 - RNAeval 1 - Partition Fx 2- MFE
  int triplet = 0; // to use triplet model(0) or not(1)
  double temp = 37;
  int pflag = 1;  
  int samples = 10; 
  int cdflag = 0;
  int tmmflag = 0;
  int tcflag = 0;
  int dflag = 0;

  sequence[0] = '0';
  secstr[0] = '0';
  if (argc>=4){// need -m mode -s sequence 
    i=1;
    
    // take in inputs
    while (i<argc){
      if(!strcmp(argv[i],"-d")){
        dflag = atoi(argv[i+1]);
        DEBUG = dflag;
      }
      else if(!strcmp(argv[i],"-m")){// mode
        mode = atoi(argv[i+1]);// store selection, its the next argument
      }      
      else if(!strcmp(argv[i],"-s")){// sequence
        len = strlen(argv[i+1]);
        for(k = 0; k < len; k++){
          sequence[k+1] = argv[i+1][k];// store sequence
        }
        sequence[strlen(argv[i+1])+1] = '\0';//end character for sequence  
        PrintStr(sequence);
      } 
      else if(!strcmp(argv[i],"-ss")){//secondary structure, necessary if mode is 0(RNAeval)
        len = strlen(argv[i+1]);
        for(k = 0; k < len; k++){
          secstr[k+1] = argv[i+1][k];// store secondary structure
        }
        secstr[strlen(argv[i+1])+1] = '\0';//end character for secondary structure  
        PrintStr(secstr);
        printf("\n\n");
      }
      else if(!strcmp(argv[i],"-tr")){// triplet model 0 for no, 1 for yes 
        triplet = atoi(argv[i+1]);// store triplet model selection
        //triplet = !triplet;
      }
      else if(!strcmp(argv[i],"-t")){
        temp = atoi(argv[i+1]);
        T = temp+273.15;
      }
      else if(!strcmp(argv[i],"-p")){
        pflag = atoi(argv[i+1]); 
      }
      else if(!strcmp(argv[i],"-cd")){
        cdflag = atoi(argv[i+1]);
      }
      else if(!strcmp(argv[i],"-ns")){
        samples = atoi(argv[i+1]);
      }
      else if(!strcmp(argv[i],"-tmm")){
        tmmflag = atoi(argv[i+1]);
      }
      else if(!strcmp(argv[i],"-tc")){
        tcflag = atoi(argv[i+1]);
      } 
      else { //if(!strcmp(argv[i],"s")) - error usage
        printf("Error: Usage: %s -m MODE -s SEQUENCE\n",argv[0]);
        printf("Options:\n");
        printf("-m : 0 (RNAeval), 1(Partition Functon), 2(MFE) 3(Sampling) 4(BP probability) 5(RNAeval and MFE)\n");
        printf("-tr: 0(Yes), 1(No): the default is 0\n");
        printf("-ss: secondary structure(needed for RNAeval)\n");
        printf("-t: temperature in C (int)\n");
        printf("-ns: number of samples\n");
        printf("-p: bp probability method 1,2,3\n");
        printf("-cd: Coaxial,Dangles 0(NO - default), 1(YES)\n");
        printf("-tmm: TMM 0(NO- default), 1(YES)\n"); 
        printf("-tc: Temperature Method 0(Traditional), 1(Matthews)\n");
        printf("-d: DEBUG 0 (defualt -no), 1 (yes)\n");
        exit(1);
      }
      i+=2;// want to go to next commandline argument flag
    }// end while of looping through commandline arguments
    
    integer_seq = TransformSeq(sequence);// transform sequence to integer array in order to use energy tables
    /*if (DEBUG){
      printf("sequence transformation to integer: \n"); 
      PrintInt(integer_seq,strlen(secstr));    
      PrintStr(secstr);
      printf("\n");
    }*/

    CreateEnergyTable(temp,tcflag); 

    start_timer();

    char* mfe_seq = (char*)malloc((len+2)*sizeof(char));
    mfe_seq[0] = '0';
    for(i=1;i<=len;i++)
      mfe_seq[i] = '.';
    mfe_seq[len+1] = '\0';
   
    switch(triplet){
      case 1: 
        switch(mode){// mode selected is RNAeval
          case 0:
            if(strlen(sequence)!=strlen(secstr)){
              printf("Error: Sequence and Secondary Structure are not same length\n");
              exit(1);
            } 
            else{//if(strlen(sequence)!=strlen(secstr))
                printf("%f\n",GetStructureEnergy(integer_seq,secstr,triplet,cdflag,tmmflag));  
            }      
            break;
          case 1://mode selected is Partition Function
            printf("%f\n",GetTPartFunc(integer_seq,len));
            break;
          case 2://MFE
            GetTMFE(integer_seq,len,mfe_seq);
            break;
          case 3:// sample
            tsample(integer_seq,len,samples);
            break;  
          case 4:// BP probability
            printf("BP probability\n");
            switch(pflag){
              case 1:
                tcalculateProbs1(integer_seq,len,samples);
                break;
              case 2:
                tcalculateProbs2(integer_seq,len,samples);
                break;
              /*case 3:
                calculateProbs3(integer_seq,len);
               break;*/
            default:
              printf("Please signinfy method using -p with option 1,2,3\n");
           }   
           break;
           case 5:// RNAeval and MFE
             if (dflag == 1){
               DEBUG = 0;
             }
             GetTMFE(integer_seq,len,mfe_seq);
             DEBUG = dflag; 
             printf("\nEval Output\n\n");
             printf("%f\n",GetStructureEnergy(integer_seq,mfe_seq,triplet,cdflag,tmmflag));
        }// switch mode  
        break; 
      case 0:
        switch(mode){// mode selected is RNAeval
          case 0:
            if(strlen(sequence)!=strlen(secstr)){
              printf("Error: Sequence and Secondary Structure are not same length\n");
              exit(1);
            }
            else{//if(strlen(sequence)!=strlen(secstr))
                printf("%f\n",GetStructureEnergy(integer_seq,secstr,triplet,cdflag,tmmflag));
            }
            break;
          case 1://mode selected is Partition Function
            printf("%f\n",GetPartFunc(integer_seq,len));
            break;
          case 2://MFE
            GetMFE(integer_seq,len,mfe_seq);
            break;
          case 3:// sample
            sample(integer_seq,len,samples);
            break;
          case 4:// BP probability
            printf("BP probability\n");
            switch(pflag){
              case 1:
                calculateProbs1(integer_seq,len,samples);
                break;
              case 2:
                calculateProbs2(integer_seq,len,samples);
                break;
              /*case 3:
                calculateProbs3(integer_seq,len);
               break;*/
            default:
              printf("Please signinfy method using -p with option 1,2,3\n");
           }
           break;
           case 5:// RNAeval and MFE
            if (dflag == 1){
               DEBUG = 0;
             }
             GetMFE(integer_seq,len,mfe_seq);
             DEBUG = dflag;
             printf("\nEval Output\n\n");
             printf("%f\n",GetStructureEnergy(integer_seq,mfe_seq,triplet,cdflag,tmmflag));
        }// switch mode   
   } 

    printf("Elapsed time = %f\n",elapsed_time());

    free(integer_seq);// free integer_seq after use
    free(mfe_seq);
  } else{// if (argc>=4); need atleast 4 arguments
        printf("Error: Usage: %s -m MODE -s SEQUENCE\n",argv[0]);
        printf("Options:\n");
        printf("-m : 0 (RNAeval), 1(Partition Functon), 2(MFE) 3(Sampling) 4(BP probability) 5(RNAeval and MFE)\n");
        printf("-tr: 0(Yes), 1(No): the default is 0\n");
        printf("-ss: secondary structure(needed for RNAeval)\n");
        printf("-t: temperature in C (int)\n");
        printf("-ns: number of samples\n");
        printf("-p: bp probability method 1,2,3\n");
        printf("-cd: Coaxial,Dangles 0(NO - default), 1(YES)\n");
        printf("-tmm: TMM 0(NO- default), 1(YES)\n");
        printf("-tc: Temperature Method 0(Traditional), 1(Matthews)\n");
        printf("-d: DEBUG 0 (defualt -no), 1 (yes)\n");
        exit(1);
  }
  return 0;
}  
