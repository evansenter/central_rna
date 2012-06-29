#include <complex>

typedef std::complex<double> dcomplex;

void FFTbor(int *intSequence, char *structure, int length, double temperature);
void solveSystem(int sequenceLength, dcomplex **rootsOfUnity, double *coefficients, double scalingFactor);
int canPair(int a, int b);
int jPairedTo(int i, int j, int *basePairs);
int jPairedIn(int i, int j, int *basePairs);
int **numBasePairs(int *basePairs, int length);
int bpBetween(int i, int j, int *basePairs);
double GetInteriorStackingAndBulgeEnergy(int i, int j, int k, int l, int *intSequence);