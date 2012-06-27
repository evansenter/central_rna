#include <time.h>
#include <string.h>
//AntonioXXX Elimino el siguiente include para compilar en EWindows
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include "random.h"
#include <stdio.h>

long get_local_seed()
{
  struct timeval t;
  unsigned short sd[3];
  long s;
  
  gettimeofday(&t,NULL);
  sd[0] = (t.tv_usec) & 0177777;
  sd[1] = getpid();
  sd[2] = getppid();
////AntonioXXX A�ado lo siguiente para no desvirtuar demasido y evitar que s sea siempre 0
  sd[2] = 1;
  s=sd[0]*sd[1]*sd[2];
  if(s>0)
    s=-s;
  return s;
//AntonioXXX Redefino esto para compilar en Windows!!
  //return 1;  
}

/*********************************************************************
  3. This random number generator is from William H. Press, et al.,
     _Numerical Recipes in C_, Second Ed. with corrections (1994), 
     p. 282.  This excellent book is available through the
     WWW at http://nr.harvard.edu/nr/bookc.html.
     The specific section concerning ran2, Section 7.1, is in
     http://cfatab.harvard.edu/nr/bookc/c7-1.ps
*********************************************************************/



/* ran2() - Return a random floating point value between 0.0 and
   1.0 (not included).  If idum is negative, a new series starts (and
   idum is made positive so that subsequent calls using an unchanged
   idum will continue in the same sequence). */

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }
  
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}

