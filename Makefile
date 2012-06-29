CFLAGS  = -Wall -c -O3
LDFLAGS = -lfftw3
CC      = g++

all: TripletPF

TripletPF: TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o fftbor.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o
	$(CC) TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o fftbor.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o $(LDFLAGS) -o TripletPF

TripletPF.o: TripletPF.cpp
	$(CC) $(CFLAGS) TripletPF.cpp

tpfunc.o: tpfunc.cpp
	$(CC) $(CFLAGS) tpfunc.cpp

tmfe.o: tmfe.cpp
	$(CC) $(CFLAGS) tmfe.cpp

tsampling.o: tsampling.cpp
	$(CC) $(CFLAGS) tsampling.cpp

tcalculateProbs1.o: tcalculateProbs1.cpp
	$(CC) $(CFLAGS) tcalculateProbs1.cpp

tcalculateProbs2.o: tcalculateProbs2.cpp
	$(CC) $(CFLAGS) tcalculateProbs2.cpp

pfunc.o: pfunc.cpp
	$(CC) $(CFLAGS) pfunc.cpp

fftbor.o: fftbor.cpp misc.h misc.cpp
	$(CC) $(CFLAGS) fftbor.cpp

mfe.o: mfe.cpp
	$(CC) $(CFLAGS) mfe.cpp

sampling.o: sampling.cpp
	$(CC) $(CFLAGS) sampling.cpp

calculateProbs1.o: calculateProbs1.cpp
	$(CC) $(CFLAGS) calculateProbs1.cpp

calculateProbs2.o: calculateProbs2.cpp
	$(CC) $(CFLAGS) calculateProbs2.cpp

energy_func.o: energy_func.cpp
	$(CC) $(CFLAGS) energy_func.cpp

energy_par.o: energy_par.cpp
	$(CC) $(CFLAGS) energy_par.cpp

enthalpy_par.o: enthalpy_par.cpp
	$(CC) $(CFLAGS) enthalpy_par.cpp

misc.o: misc.cpp
	$(CC) $(CFLAGS) misc.cpp

random.o: random.cpp
	$(CC) $(CFLAGS) random.cpp

clean:
	rm -rf *o TripletPF
# DO NOT DELETE
