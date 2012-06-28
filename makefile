CFLAGS  = -c -O3
LDFLAGS = -lfftw3 -L.
BINDIR  = /usr/local/bin # Change this to the BINDIR
CC      = g++

all: TripletPF

TripletPF: TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o
	$(CC) -o TripletPF TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o -lm -g

TripletPF.o: TripletPF.cpp
	$(CC) -c TripletPF.cpp -g

tpfunc.o: tpfunc.cpp
	$(CC) -c tpfunc.cpp -g

tmfe.o: tmfe.cpp
	$(CC) -c tmfe.cpp -g

tsampling.o: tsampling.cpp
	$(CC) -c tsampling.cpp -g

tcalculateProbs1.o: tcalculateProbs1.cpp
	$(CC) -c tcalculateProbs1.cpp -g

tcalculateProbs2.o: tcalculateProbs2.cpp
	$(CC) -c tcalculateProbs2.cpp -g

pfunc.o: pfunc.cpp
	$(CC) -c pfunc.cpp -g

mfe.o: mfe.cpp
	$(CC) -c mfe.cpp -g

sampling.o: sampling.cpp
	$(CC) -c sampling.cpp -g

calculateProbs1.o: calculateProbs1.cpp
	$(CC) -c calculateProbs1.cpp -g

calculateProbs2.o: calculateProbs2.cpp
	$(CC) -c calculateProbs2.cpp -g

energy_func.o: energy_func.cpp
	$(CC) -c energy_func.cpp -g

energy_par.o: energy_par.cpp
	$(CC) -c energy_par.cpp

enthalpy_par.o: enthalpy_par.cpp
	$(CC) -c enthalpy_par.cpp

misc.o: misc.cpp
	$(CC) -c misc.cpp

random.o: random.cpp
	$(CC) -c random.cpp

clean:
	rm -rf *o TripletPF
# DO NOT DELETE
