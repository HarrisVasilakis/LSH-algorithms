COMP = g++
CFLAGS = -Wall -std=c++11
ODIR = ./obj

# Dependencies
_DEPS = datasetarray ghashfunction class clusteringfunctions cubehash hypercubefunctions
DEPS := $(foreach D,$(_DEPS),./$D.h)
INCLUDE_HEADERS = $(patsubst %,-I ./%,$(_DEPS))

# Object Files
_LSH_OBJS = lsh datasetarray ghashfunction class
LSH_OBJS = $(patsubst %,$(ODIR)/%.o,$(_LSH_OBJS))

_CUBE_OBJS = cube datasetarray ghashfunction class cubehash hypercubefunctions
CUBE_OBJS = $(patsubst %,$(ODIR)/%.o,$(_CUBE_OBJS))

_CLUSTER_OBJS = clustering datasetarray ghashfunction class clusteringfunctions
CLUSTER_OBJS = $(patsubst %,$(ODIR)/%.o,$(_CLUSTER_OBJS))

# Executables
EXEC = lsh cluster cube

# Compilation rules for executables
all: $(EXEC)

lsh: $(LSH_OBJS)
	$(COMP) -o $@ $^ $(CFLAGS)
	@echo "Executable file \"$@\" created";

cube: $(CUBE_OBJS)
	$(COMP) -o $@ $^ $(CFLAGS)
	@echo "Executable file \"$@\" created";

clustering: $(CLUSTER_OBJS)
	$(COMP) -o $@ $^ $(CFLAGS)
	@echo "Executable file \"$@\" created";

$(ODIR)/lsh.o: ./lsh.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/cube.o: ./cube.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/clustering.o: ./clustering.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)	

# Compilation Rules for all other files
$(ODIR)/datasetarray.o: ./datasetarray.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/ghashfunction.o: ./ghashfunction.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/class.o: ./class.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/clusteringfunctions.o: ./clusteringfunctions.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/cubehash.o: ./cubehash.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

$(ODIR)/hypercubefunctions.o: ./hypercubefunctions.cpp $(DEPS) | $(ODIR)
	$(COMP) -c -o $@ $< $(INCLUDE_HEADERS) $(CFLAGS)

# Object files dir
$(ODIR):
	mkdir -p $@

# make clean
.PHONY: clean

clean:
	rm -rf $(ODIR)/
	rm -f $(EXEC)