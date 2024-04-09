# Build tools
#NVCC = nvcc $(ARCH) -ccbin clang++
#CXX = clang++
ARCH = -arch=sm_52

NVCC = nvcc $(ARCH)
CXX = g++ -std=c++20

#QBS_REGION = 4
#D = -D QBS_REGION=$(QBS_REGION)
OPS_BLOCK=300

# here are all the objects
GPUOBJS = kernel.o
OBJS = dgm.o common.o gates.o lib_general.o lib_shor.o lib_grover.o
OBJS_ONLY_CPU = dgm_cpu.o common.o gates.o lib_general.o lib_shor.o lib_grover.o

# make and compile

# executables

shor: shor.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o shor.out shor.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

shor_cpu: shor.o $(OBJS_ONLY_CPU)
	$(CXX) $(CXXCPUFLAGS) -o shor_cpu.out shor.o $(OBJS) -fopenmp

grover: grover.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o grover.out grover.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

general: general.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o general.out general.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

general_cpu: general.o $(OBJS_ONLY_CPU)
	$(CXX) $(CXXCPUFLAGS) -o general_cpu.out general.o $(OBJS) -fopenmp
 
# objects

dgm.o: dgm.cu
	$(NVCC) -c dgm.cu -Xcompiler "-fopenmp -O3 -fcx-limited-range"

dgm_cpu.o: dgm.cu
	$(CXX) -c -x c++ dgm.cu -O3 -fopenmp -fcx-limited-range -D ONLY_CPU

kernel.o: kernel.cu
	$(NVCC) -c -D OPS_BLOCK=$(OPS_BLOCK) kernel.cu

gates.o: gates.cpp
	$(CXX) -c gates.cpp

common.o: common.cpp
	$(CXX) -c common.cpp

lib_general.o: lib_general.cpp
	$(CXX) -c lib_general.cpp

lib_shor.o: lib_shor.cpp
	$(CXX) -c lib_shor.cpp

lib_grover.o: lib_grover.cpp
	$(CXX) -c lib_grover.cpp
 
grover.o: grover.cpp
	$(CXX) -c grover.cpp -fopenmp

shor.o: shor.cpp
	$(CXX) -c shor.cpp -fopenmp
 
general.o: general.cpp
	$(CXX) -c general.cpp -fopenmp

# test

test_cpu: test.o lib_test.o
	$(CXX) -o test_cpu.out test.o lib_test.o

test.o: test.cpp
	$(CXX) -c test.cpp

lib_test.o: lib_test.cu
	$(CXX) -c -x c++ lib_test.cu

# clean up

clean:
	rm *.o *.out
