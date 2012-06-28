all: TripletPF

#MatthewsRNAeval: MatthewsRNAeval.o energy_func.o energy_par.o misc.o
#	gcc -o MatthewsRNAeval MatthewsRNAeval.o energy_func.o energy_par.o misc.o -lm 

TripletPF: TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o
	g++ -o TripletPF TripletPF.o tpfunc.o tmfe.o tsampling.o tcalculateProbs1.o tcalculateProbs2.o pfunc.o mfe.o sampling.o calculateProbs1.o calculateProbs2.o energy_func.o energy_par.o enthalpy_par.o misc.o random.o -lm -g

#MatthewsRNAeval.o: MatthewsRNAeval.cpp
#	gcc -c MatthewsRNAeval.cpp

TripletPF.o: TripletPF.cpp
	g++ -c TripletPF.cpp -lm -g

tpfunc.o: tpfunc.cpp
	g++ -c tpfunc.cpp -lm -g

tmfe.o: tmfe.cpp
	g++ -c tmfe.cpp -lm -g

tsampling.o: tsampling.cpp
	g++ -c tsampling.cpp -lm -g

tcalculateProbs1.o: tcalculateProbs1.cpp
	g++ -c tcalculateProbs1.cpp -lm -g

tcalculateProbs2.o: tcalculateProbs2.cpp
	g++ -c tcalculateProbs2.cpp -lm -g

pfunc.o: pfunc.cpp
	g++ -c pfunc.cpp -lm -g

mfe.o: mfe.cpp
	g++ -c mfe.cpp -lm -g

sampling.o: sampling.cpp
	g++ -c sampling.cpp -lm -g

calculateProbs1.o: calculateProbs1.cpp
	g++ -c calculateProbs1.cpp -lm -g

calculateProbs2.o: calculateProbs2.cpp
	g++ -c calculateProbs2.cpp -lm -g

#tcalculateProbs3.o: tcalculateProbs3.cpp
#	gcc -c tcalculateProbs3.cpp -lm -g

energy_func.o: energy_func.cpp
	g++ -c energy_func.cpp -g -lm

energy_par.o: energy_par.cpp
	g++ -c energy_par.cpp

enthalpy_par.o: enthalpy_par.cpp
	g++ -c enthalpy_par.cpp

misc.o: misc.cpp
	g++ -c misc.cpp

random.o: random.cpp
	g++ -c random.cpp -lm

clean:
	rm -rf *o TripletPF
