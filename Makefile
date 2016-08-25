CXXFLAGS=-O3 -D _RANDOM -D _FORCEONWALLS
LDFLAGS=-lpthread

current: main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o forceOnWalls.o write_vtk.o generateInitialState.o
	g++ -o run main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o forceOnWalls.o write_vtk.o generateInitialState.o $(LDFLAGS)
prog: main_prog.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_prog main_prog.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o $(LDFLAGS)
benchmark: main_benchmark.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o 
	g++ -o parallel_benchmark main_benchmark.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o $(LDFLAGS)
main.o: main.cpp
	g++ -o main.o -c main.cpp $(CXXFLAGS)
main_prog.o: main_prog.cpp
	g++ -o main_prog.o -c main_prog.cpp $(CXXFLAGS)
main_benchmark.o:
	g++ -o main_benchmark.o -c main_benchmark.cpp $(CXXFLAGS)
initialize_lattice_arrays.o: initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c initialize_lattice_arrays.cpp $(CXXFLAGS)
streamCollCompute.o: streamCollCompute.cpp
	g++ -o streamCollCompute.o -c streamCollCompute.cpp $(CXXFLAGS)
domain_noSlipWalls.o: domain_noSlipWalls.cpp
	g++ -o domain_noSlipWalls.o -c domain_noSlipWalls.cpp $(CXXFLAGS)
square.o: square.cpp
	g++ -o square.o -c square.cpp $(CXXFLAGS)
force.o: force.cpp
	g++ -o force.o -c force.cpp $(CXXFLAGS)
forceOnWalls.o: forceOnWalls.cpp
	g++ -o forceOnWalls.o -c forceOnWalls.cpp $(CXXFLAGS)
write_vtk.o: write_vtk.cpp
	g++ -o write_vtk.o -c write_vtk.cpp $(CXXFLAGS)
generateInitialState.o: generateInitialState.cpp
	g++ -o generateInitialState.o -c generateInitialState.cpp $(CXXFLAGS)
clean:
	rm -rf *.o
mrproper: clean
	rm -rf run
	rm -rf benchmark
	rm -rf prog

