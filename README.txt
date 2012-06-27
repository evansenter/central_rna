Energy Tables are based on integer representations of nucleotides based on following:

G-> 1
A-> 3
C-> 7
U-> 12

index number      // 0 1  2 3  4 5 6  7 8 9 10 11  12 13 14 15 16  17 18 19 20 21  22  23 24  25  26  27  28
extern int d1[29] = {1,4,-1,0,-1,2,2, 9,6,0,10, 3, -1, 7, 7, 4, 8, 11, 3, 6,-1, 1, -1,  5, 5, 12, 13, 14, 15};

int d2[13] = {0,13,0,14,0,0,0,15,0,0,0,0,16};


  00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 
  AU UA GC CG GU UG AC AG CA CU GA UC GG AA CC UU

  00 01 02 03 04 05 06 07
  0A 0U 0C 0G A0 U0 C0 G0


To represent basepairs use the following:

index = i-j+12+((i==j)*d2[i]);

For example (CG):

C = 7
G = 1

index = 7-1+12+(5==1)*d2[5]
index = 7-1+12 = 16

Now look in d1:
d1[18] = 3

In the energy tables in energy_par.c to get to values corresponding to
basepair CG you will need to look at index 3. For example the stacking energy
for CG on CG one would look at stack37[3][3] to get the value.


**************************************

To represent dangles the following is done:

A0 represents a 5' dangle:

A = 3
0 = 0 ( the j value is 0)

index = 3-0+12+(3==0)*d2[0]
index = 15

d1[15] = 4 and one can see above that 4 represents A0

*****************************************

To represent a dinucleotide of the same composition such as AA,CC,GG,UU (used
for terminal mismatches) the following is done:

For example CC:

C = 7
C = 7

index = 7-7+12+(7==7)*d2[7]
index = 12+1*15 = 27

d1[27] = 14

From above one can see that 14 corresponds to CC; and in energy tables an
index of 14 will correspond to CC

******************************************

For instances such as 2x2 Internal Loop of 1x2 Internal Loop row columnar
multiplication

For example: 

3' 5'
G--C
 AA
 UG
A--U
5' 3'

In this 2x2 loop AU is a basepair and GC is a basepair. The 2x2 loop is made
of UG(not base paired) with AA completing the internal loop. The energy tables
are constructed so that one can represent the AU and GC base pairs as row
components and the UG and AA internal bases as columnar components in a 2D
array.


In this 2D array the row will be the basepairs. To represent AU opening with a
GC closing we shall do the following:

AU index = 3-12+12+(3==12)*d2[3]
         = 3
         d1[3] = 0 


GC index = 1-7+12+(1==7)*d2[1]
         = 6
         d1[6] = 2


To represnt AU and GC one would do the following:

AU index *6 + GC inde:
0*6+2 = 2

In the table for Internal22 the index 2 so the third row applies to values for
AU opening the internal loop and GC closing the internal loop

Similarly to represent in the internal bases:

UG index = 12-1+12+(12==1)*d2[12]
         = 23
         d1[23] = 5 

AA index = 3-3+12+(3==3)*d2[3]
         = 12+14 = 26
         d1[26] = 13 
To represent UG and AA one would do the following:

UG index * 16 + AA index:
5*16+13 = 93

This provides the column component to look up in internal22 in enery_par.c:

Therefore if one looks at internal22[2][93] one will get the energy for the
internal loop displayed above
  
