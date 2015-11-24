CXXFLAGS=-O3
LDFLAGS=-lpthread

current: main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o write_vtk.o
	g++ -o run main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o write_vtk.o $(LDFLAGS)
benchmark: main_benchmark.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o write_vtk.o
	g++ -o parallel_benchmark main_benchmark.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o write_vtk.o $(LDFLAGS)
main.o: main.cpp
	g++ -o main.o -c main.cpp $(CXXFLAGS)
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
write_vtk.o: write_vtk.cpp
	g++ -o write_vtk.o -c write_vtk.cpp $(CXXFLAGS)
clean:
	rm -rf *.o
mrproper: clean
	rm -rf run
	rm -rf benchmark

