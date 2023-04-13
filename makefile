# Build tools
#NVCC = nvcc $(ARCH) -ccbin clang++
#CXX = clang++
ARCH = -arch=sm_52

NVCC = nvcc $(ARCH)
CXX = g++
GCC = gcc

BIN=./bin
SRC=./src
INCLUDE=./include

#QBS_REGION = 4
#D = -D QBS_REGION=$(QBS_REGION)
OPS_BLOCK=300

# here are all the objects
GPUOBJS = $(BIN)/kernel.o
OBJS = $(BIN)/dgm.o $(BIN)/common.o $(BIN)/gates.o $(BIN)/lib_hadamard.o $(BIN)/lib_shor.o $(BIN)/lib_grover.o


# make and compile

all: shor grover hadamard

# executables

shor: $(BIN)/shor.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o $(BIN)/shor.out $(BIN)/shor.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp" -I$(INCLUDE)

grover: $(BIN)/grover.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o $(BIN)/grover.out $(BIN)/grover.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp" -I$(INCLUDE)

hadamard: $(BIN)/hadamard.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o $(BIN)/hadamard.out $(BIN)/hadamard.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp" -I$(INCLUDE)
 
# objects

$(BIN)/dgm.o: $(SRC)/dgm.cu
	$(NVCC) -c $< -Xcompiler "-fopenmp -O3 -fcx-limited-range" -I$(INCLUDE) -o $@

$(BIN)/kernel.o: $(SRC)/kernel.cu
	$(NVCC) -c -D OPS_BLOCK=$(OPS_BLOCK) $< -I$(INCLUDE) -o $@

$(BIN)/gates.o: $(SRC)/gates.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@

$(BIN)/common.o: $(SRC)/common.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@

$(BIN)/lib_hadamard.o: $(SRC)/lib_hadamard.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@

$(BIN)/lib_shor.o: $(SRC)/lib_shor.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@

$(BIN)/lib_grover.o: $(SRC)/lib_grover.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@
 
$(BIN)/grover.o: $(SRC)/grover.cpp
	$(CXX) -c $< -fopenmp -I$(INCLUDE) -o $@

$(BIN)/shor.o: $(SRC)/shor.cpp
	$(CXX) -c $< -fopenmp -I$(INCLUDE) -o $@
 
$(BIN)/hadamard.o: $(SRC)/hadamard.cpp
	$(CXX) -c $< -fopenmp -I$(INCLUDE) -o $@

clean:
	rm $(BIN)/*
