all: TripletPF

#MatthewsRNAeval: MatthewsRNAeval.o energy_func.o energy_par.o misc.o
#	gcc -o MatthewsRNAeval MatthewsRNAeval.o energy_func.o energy_par.o misc.o -lm 

TripletPF: TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o
	gcc -o TripletPF TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o -lm -g

#MatthewsRNAeval.o: MatthewsRNAeval.c
#	gcc -c MatthewsRNAeval.c

TripletPF.o: TripletPF.c
	gcc -c TripletPF.c -lm -g

tpfunc.o: tpfunc.c
	gcc -c tpfunc.c -lm -g

tmfe.o: tmfe.c
	gcc -c tmfe.c -lm -g

tsampling.o: tsampling.c
	gcc -c tsampling.c -lm -g

tcalculateProbs1.o: tcalculateProbs1.c
	gcc -c tcalculateProbs1.c -lm -g

tcalculateProbs2.o: tcalculateProbs2.c
	gcc -c tcalculateProbs2.c -lm -g

pfunc.o: pfunc.c
	gcc -c pfunc.c -lm -g

mfe.o: mfe.c
	gcc -c mfe.c -lm -g

sampling.o: sampling.c
	gcc -c sampling.c -lm -g

calculateProbs1.o: calculateProbs1.c
	gcc -c calculateProbs1.c -lm -g

calculateProbs2.o: calculateProbs2.c
	gcc -c calculateProbs2.c -lm -g

#tcalculateProbs3.o: tcalculateProbs3.c
#	gcc -c tcalculateProbs3.c -lm -g

energy_func.o: energy_func.c
	gcc -c energy_func.c -g -lm

energy_par.o: energy_par.c
	gcc -c energy_par.c

enthalpy_par.o: enthalpy_par.c
	gcc -c enthalpy_par.c

misc.o: misc.c
	gcc -c misc.c

random.o: random.c
	gcc -c random.c -lm

clean:
	rm -rf *o 
