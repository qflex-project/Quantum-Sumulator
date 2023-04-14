# Build tools

CXX = clang++
ARCH = sm_70
NVCC = nvcc -arch=$(ARCH) #-ccbin $(CXX)

CXX_ARGS = -fopenmp -Ofast
NVCC_ARGS = -Xcompiler "$(CXX_ARGS)"
 
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
	$(NVCC) $< $(OBJS) $(GPUOBJS) $(NVCC_ARGS) -I$(INCLUDE) -o $(BIN)/shor.out

grover: $(BIN)/grover.o $(OBJS) $(GPUOBJS)
	$(NVCC) $< $(OBJS) $(GPUOBJS) $(NVCC_ARGS) -I$(INCLUDE) -o $(BIN)/grover.out

hadamard: $(BIN)/hadamard.o $(OBJS) $(GPUOBJS)
	$(NVCC) $< $(OBJS) $(GPUOBJS) $(NVCC_ARGS) -I$(INCLUDE) -o $(BIN)/hadamard.out
 
# objects

$(BIN)/common.o: $(SRC)/common.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/dgm.o: $(SRC)/dgm.cu
	$(NVCC) -c $< $(NVCC_ARGS) -I$(INCLUDE) -o $@

$(BIN)/gates.o: $(SRC)/gates.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/kernel.o: $(SRC)/kernel.cu
	$(NVCC) -c $< $(NVCC_ARGS) -I$(INCLUDE) -D OPS_BLOCK=$(OPS_BLOCK) -o $@

$(BIN)/lib_grover.o: $(SRC)/lib_grover.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/lib_hadamard.o: $(SRC)/lib_hadamard.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/lib_shor.o: $(SRC)/lib_shor.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/grover.o: $(SRC)/grover.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/hadamard.o: $(SRC)/hadamard.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

$(BIN)/shor.o: $(SRC)/shor.cpp
	$(CXX) -c $< $(CXX_ARGS) -I$(INCLUDE) -o $@

clean:
	rm $(BIN)/*
