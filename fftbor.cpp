#include "fftbor.h"
#include "misc.h"
#include "energy_func.h"
#include "energy_par.h"
#include <stdio.h>
#include <iostream>
#include <fftw3.h>
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define PRECISION 4
#define FFTW_REAL 0
#define FFTW_IMAG 1
#define STRUCTURE_COUNT 0
#define FFTBOR_DEBUG 0
#define WORKING_COUNTER 1

extern int DEBUG;

void FFTbor(int *intSequence, char *structure, int length, double temperature) {
  // Variable declarations.
  int i, j, d, k, l, delta, root;
  int *basePairs, **bpCount;
  double RT, energy, scalingFactor;
  dcomplex x;
  
  dcomplex **Z            = new dcomplex*[length + 1];
  dcomplex **ZB           = new dcomplex*[length + 1];
  dcomplex **ZM           = new dcomplex*[length + 1];
  dcomplex **rootsOfUnity = new dcomplex*[length + 1];
  double    *coefficients = new double[length + 1];
  
  // RT        = 0.0019872370936902486 * (temperature + 273.15) * 100; // 0.01 * (kcal K)/mol
  RT        = R * T;
  basePairs = GetBPList(structure, length + 1);
  bpCount   = numBasePairs(basePairs, length);
  
  // Matrix allocation.
  for (i = 0; i <= length; ++i) {
    Z[i]               = new dcomplex[length + 1];
    ZB[i]              = new dcomplex[length + 1];
    ZM[i]              = new dcomplex[length + 1];
    rootsOfUnity[i]    = new dcomplex[2];
    rootsOfUnity[i][0] = dcomplex(cos(2 * M_PI * i / (length + 1)), sin(2 * M_PI * i / (length + 1)));
  }
  
  // Start main recursions (root <= round(length / 2.0) is an optimization for roots of unity).
  for (root = 0; root <= round(length / 2.0); ++root) {
    // Flush the matrices.
    for (i = 0; i <= length; ++i) {
      for (j = 0; j <= length; ++j) {
        if (i > 0 && j > 0 && abs(j - i) <= MIN_PAIR_DIST) {
       Z[i][j] = ONE_C;
     } else {
       Z[i][j] = ZERO_C;
     }
     
        ZB[i][j] = ZERO_C;
        ZM[i][j] = ZERO_C;
      }
    }
    
    x = rootsOfUnity[root][0];
    
    // ****************************************************************************
    // Main recursions
    // ****************************************************************************
    for (d = MIN_PAIR_DIST + 1; d < length; ++d) {
      for (i = 1; i <= length - d; ++i) {
        j = i + d;
      
        if (canPair(intSequence[i], intSequence[j])) {
          // ****************************************************************************
          // Solve ZB 
          // ****************************************************************************
          // In a hairpin, [i + 1, j - 1] unpaired.
          energy    = GetHairpinEnergy(i, j, intSequence);
          delta     = bpCount[i][j] + jPairedTo(i, j, basePairs);
          ZB[i][j] += pow(x, delta) * exp(-energy / RT);
  
          if (STRUCTURE_COUNT) {
            ZB[j][i] += 1;
          }
  
          // Interior loop / bulge / stack / multiloop.
          for (k = i + 1; k < j - MIN_PAIR_DIST; ++k) {
            for (l = MAX(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
              if (canPair(intSequence[k], intSequence[l])) {
                 // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
                 energy    = GetInteriorStackingAndBulgeEnergy(i, j, k, l, intSequence);
                 delta     = bpCount[i][j] - bpCount[k][l] + jPairedTo(i, j, basePairs);
                 ZB[i][j] += (ZB[k][l] * pow(x, delta) * exp(-energy / RT));
  
                 if (STRUCTURE_COUNT) {
                   ZB[j][i] += ZB[l][k];
                 }
  
                if (k > i + MIN_PAIR_DIST + 2) {
                  // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1).
                  energy    = MultiloopA + 2 * MultiloopB + MultiloopC * (j - l - 1);
                  delta     = bpCount[i][j] - bpCount[i + 1][k - 1] - bpCount[k][l] + jPairedTo(i, j, basePairs);
                  ZB[i][j] += ZM[i + 1][k - 1] * ZB[k][l] * pow(x, delta) * exp(-energy / RT);
  
                  if (STRUCTURE_COUNT) {
                    ZB[j][i] += ZM[k - 1][i + 1] * ZB[l][k];
                  }
                }
              }
            }
          }
        }
        
        // ****************************************************************************
        // Solve ZM
        // ****************************************************************************
        energy    = 0;// P->MLbase;
        delta     = jPairedIn(i, j, basePairs);
        ZM[i][j] += ZM[i][j - 1] * pow(x, delta) * exp(-energy / RT);
  
        if (STRUCTURE_COUNT) {
          ZM[j][i] += ZM[j - 1][i];
        }
  
        for (k = i; k < j - MIN_PAIR_DIST; ++k) {
          if (canPair(intSequence[k], intSequence[j])) {
            // Only one stem.
            energy    = MultiloopB + MultiloopC * (k - i);
            delta     = bpCount[i][j] - bpCount[k][j];
            ZM[i][j] += ZB[k][j] * pow(x, delta) * exp(-energy / RT);
  
            if (STRUCTURE_COUNT) {
              ZM[j][i] += ZB[j][k];
            }
  
            // More than one stem.
            if (k > i + MIN_PAIR_DIST + 2) {
              // I believe this needs a MultiloopC penalty as well, but since it's set to 0 it's not a big deal.
              energy    = MultiloopB;
              delta     = bpCount[i][j] - bpCount[i][k - 1] - bpCount[k][j];
              ZM[i][j] += ZM[i][k - 1] * ZB[k][j] * pow(x, delta) * exp(-energy / RT);
  
              if (STRUCTURE_COUNT) {
                ZM[j][i] += ZM[k - 1][i] * ZB[j][k];
              }
            }
          }
        }
        
        // **************************************************************************
        // Solve Z
        // **************************************************************************
        delta    = jPairedIn(i, j, basePairs);
        Z[i][j] += Z[i][j - 1] * pow(x, delta);
  
        if (STRUCTURE_COUNT) {
          Z[j][i] += Z[j - 1][i];
        }
  
        for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
          // (k, j) is the rightmost base pair in (i, j).
          if (canPair(intSequence[k], intSequence[j])) {
            // if d1[k, j] == 2 or 3, it's not an AU or GU pair.
            energy = d1[GetIndex(intSequence[k], intSequence[j])] == 2 || d1[GetIndex(intSequence[k], intSequence[j])] == 3 ? 0 : GUAU_penalty;
            
            if (k == i) {
              delta    = bpCount[i][j] - bpCount[k][j];
              Z[i][j] += ZB[k][j] * pow(x, delta) * exp(-energy / RT);
  
              if (STRUCTURE_COUNT) {
                Z[j][i] += ZB[j][k];
              }
            } else {
              delta    = bpCount[i][j] - bpCount[i][k - 1] - bpCount[k][j];
              Z[i][j] += Z[i][k - 1] * ZB[k][j] * pow(x, delta) * exp(-energy / RT);
  
              if (STRUCTURE_COUNT) {
                Z[j][i] += Z[k - 1][i] * ZB[j][k];
              }
            }
          }
        }
      }
    }
    
    rootsOfUnity[root][1] = Z[1][length];
    
    if (!root) {
      scalingFactor = Z[1][length].real();
    }
  
    if (WORKING_COUNTER) {
      std::cout << "." << std::flush;
    }
  }
  
  if (WORKING_COUNTER) {
    printf("\n");
  }
  
  // Optimization leveraging complementarity of roots of unity.
  if (length % 2) {
    i = root - 2;
  } else {
    i = root - 1;
  }
  
  for (; root <= length && i > 0; --i, ++root) {
    rootsOfUnity[root][1] = dcomplex(rootsOfUnity[i][1].real(), -rootsOfUnity[i][1].imag());
  }
  
  if (STRUCTURE_COUNT) {
    printf("Number of structures: %.0f\n", Z[length][1].real());
  }
  
  solveSystem(length, rootsOfUnity, coefficients, scalingFactor);
}

void solveSystem(int length, dcomplex **rootsOfUnity, double *coefficients, double scalingFactor) {
  int i, j;
  dcomplex sum = ZERO_C;
  
	if (FFTBOR_DEBUG) {
		std::cout << "START ROOTS AND SOLUTIONS" << std::endl << std::endl;
		
		for (i = 0; i <= length; ++i) {
      for (j = 0; j <= 1; ++j) {
        printf("%+.15f, %-+25.15f", rootsOfUnity[i][j].real(), rootsOfUnity[i][j].imag());
      }
      std::cout << std::endl;
    }
    
	  std::cout << "END ROOTS AND SOLUTIONS" << std::endl << std::endl;
	}

  fftw_complex signal[length + 1];
  fftw_complex result[length + 1];
  
  for (i = 0; i <= length; i++) {
    signal[i][FFTW_REAL] = (pow(10, PRECISION) * rootsOfUnity[i][1].real()) / scalingFactor;
    signal[i][FFTW_IMAG] = (pow(10, PRECISION) * rootsOfUnity[i][1].imag()) / scalingFactor;
  }
  
  fftw_plan plan = fftw_plan_dft_1d(length + 1, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  for (i = 0; i <= length; i++) {
    coefficients[i] = PRECISION == 0 ? result[i][FFTW_REAL] / (length + 1) : pow(10.0, -PRECISION) * static_cast<int>(result[i][FFTW_REAL] / (length + 1));
    sum            += coefficients[i];
    
    std::cout << i << "\t" << coefficients[i] << std::endl;
  }
  
	printf("Scaling factor (Z{1, n}): %.15f\n", scalingFactor);
  std::cout << "Sum: " << sum << std::endl;
}

inline int canPair(int a, int b) {
  // Takes the integer representation of two nucleotides.
  if (DEBUG) {
    if (!((a == 1 || a == 3 || a == 7 || a == 12) && (b == 1 || b == 3 || b == 7 || b == 12))) {
      printf("WARNING: calling canPair(%d, %d) with arguments that aren't the integer encoding of nucleotides. This is probably an error\n", a, b);
    }
  }
  
  return d1[GetIndex(a, b)] < 6;
}

inline int jPairedTo(int i, int j, int *basePairs) {
  return basePairs[i] == j ? -1 : 1;
}

inline int jPairedIn(int i, int j, int *basePairs) {
  return basePairs[j] >= i && basePairs[j] < j ? 1 : 0;
}

int **numBasePairs(int *basePairs, int length) {
  int d, i, j;
  
  int **bpCount = Allocate2Dmatrix(length + 1, length + 1);
  
  for (d = MIN_PAIR_DIST + 1; d < length; d++) {
    for (i = 1; i <= length - d; i++) {
      j = i + d;
      bpCount[i][j] = bpBetween(i, j, basePairs);
    }
  }
  
  return bpCount;
}

int bpBetween(int i, int j, int *basePairs) {
  int k, n = 0;
  for (k = i; k <= j; k++) {
    if (k < basePairs[k] && basePairs[k] <= j) {
      n++;
    }
  }
  
  return n;
}

inline double GetInteriorStackingAndBulgeEnergy(int i, int j, int k, int l, int *intSequence) {
  if (k == i + 1 && l == j - 1) {
    // Stacking.
    return GetStackEnergy(i, j, k, l, intSequence);
  } else if (l == j - 1) {
    // Left bulge.
    return GetLeftBulgeEnergy(i, j, k, intSequence);
  } else if (k == i + 1) {
    // Right bulge.
    return GetRightBulgeEnergy(i, j, l, intSequence);
  } else {
    // Interior loop.
    return GetInternalEnergy(i, j, k, l, intSequence);
  }
}