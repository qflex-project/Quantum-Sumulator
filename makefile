# Build tools

CXX = mpic++
ARCH = sm_70
NVCC = nvcc -arch=$(ARCH) -std=c++03 -ccbin $(CXX)

CXX_ARGS = -fopenmp -std=c++03
NVCC_ARGS = -Xcompiler "$(CXX_ARGS)"
 
BIN=./bin
SRC=./src
INCLUDE=./include

#QBS_REGION = 4
#D = -D QBS_REGION=$(QBS_REGION)
OPS_BLOCK=300

# here are all the objects
GPUOBJS = $(BIN)/gpu.o
OBJS =  $(BIN)/common.o $(BIN)/pt.o $(BIN)/cpu.o $(BIN)/pcpu.o $(BIN)/hybrid.o $(BIN)/dgm.o $(BIN)/gates.o $(BIN)/lib_hadamard.o $(BIN)/lib_shor.o $(BIN)/lib_grover.o

# make and compile
all: shor grover hadamard

# debug
debug: CXX_ARGS += -DDEBUG -g -O0
debug: all

# release
release: CXX_ARGS += -O3
release: all

# executables

shor: $(BIN)/shor.o $(OBJS) $(GPUOBJS)
	$(NVCC) $< $(OBJS) $(GPUOBJS) -I$(INCLUDE) -o $(BIN)/shor.out $(NVCC_ARGS) 

grover: $(BIN)/grover.o $(OBJS) $(GPUOBJS)
	$(NVCC) $< $(OBJS) $(GPUOBJS) -I$(INCLUDE) -o $(BIN)/grover.out $(NVCC_ARGS) 

hadamard: $(BIN)/hadamard.o $(OBJS) $(GPUOBJS)
	$(NVCC) $< $(OBJS) $(GPUOBJS) -I$(INCLUDE) -o $(BIN)/hadamard.out $(NVCC_ARGS) 
 
# objects

$(BIN)/pt.o: $(SRC)/pt.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/common.o: $(SRC)/common.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/cpu.o: $(SRC)/cpu.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS)

$(BIN)/pcpu.o: $(SRC)/pcpu.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/dcpu.o: $(SRC)/dcpu.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/hybrid.o: $(SRC)/hybrid.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/dgm.o: $(SRC)/dgm.cu
	$(NVCC) -c $< -I$(INCLUDE) -o $@ $(NVCC_ARGS) 

$(BIN)/gates.o: $(SRC)/gates.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/gpu.o: $(SRC)/gpu.cu
	$(NVCC) -c $< -I$(INCLUDE) -D OPS_BLOCK=$(OPS_BLOCK) -o $@ $(NVCC_ARGS) 

$(BIN)/lib_grover.o: $(SRC)/lib_grover.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/lib_hadamard.o: $(SRC)/lib_hadamard.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/lib_shor.o: $(SRC)/lib_shor.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/grover.o: $(SRC)/grover.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/hadamard.o: $(SRC)/hadamard.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

$(BIN)/shor.o: $(SRC)/shor.cpp
	$(CXX) -c $< -I$(INCLUDE) -o $@ $(CXX_ARGS) 

clean:
	rm $(BIN)/*