#include <complex>

typedef std::complex<double> dcomplex;

void FFTbor(int *intSequence, char *structure, int length, double temperature);
void solveSystem(int sequenceLength, dcomplex **rootsOfUnity, double *coefficients, double scalingFactor);
inline int canPair(int a, int b);
inline int jPairedTo(int i, int j, int *basePairs);
inline int jPairedIn(int i, int j, int *basePairs);
int **numBasePairs(int *basePairs, int length);
int bpBetween(int i, int j, int *basePairs);
inline double GetInteriorStackingAndBulgeEnergy(int i, int j, int k, int l, int *intSequence);