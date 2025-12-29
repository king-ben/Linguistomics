
# define the C compiler to use
CC = clang++

# define any compile-time flags
CFLAGS = -Wall -g -std=c++17 -Ofast

# define any directories containing header files other than /usr/include
INCLUDES = -I/home/newhall/include  -I../include

# define library paths in addition to /usr/lib
LFLAGS =

# define any libraries to link into executable
LIBS = -lm

# define the C source files
SRCS  =   main.cpp \
Alignment.cpp \
AlignmentProposal.cpp \
IntVector.cpp \
JphMatrix.cpp \
LikelihoodCalculator.cpp \
MatrixMath.cpp \
Mcmc.cpp \
Model.cpp \
Msg.cpp \
Node.cpp \
NodeFactory.cpp \
NodeSet.cpp \
Parameter.cpp \
ParameterAlignment.cpp \
ParameterEquilibirumFrequencies.cpp \
ParameterExchangabilityRates.cpp \
ParameterIndelGammaShape.cpp \
ParameterIndelRates.cpp \
ParameterRatesGammaShape.cpp \
ParameterTree.cpp \
Partition.cpp \
Probability.cpp \
RandomVariable.cpp \
RateMatrix.cpp \
RateMatrixHelper.cpp \
RbBitSet.cpp \
SiteLikelihood.cpp \
Subset.cpp \
TaxonPair.cpp \
Threads.cpp \
TransitionProbabilities.cpp \
Tree.cpp \
UpdateInfo.cpp \
UserSettings.cpp

# define the object files 
OBJS = $(SRCS:.c=.o)

# define the executable file 
MAIN = twister

# generic code
.PHONY: depend clean

all: $(MAIN)
	@echo	The program $(MAIN) has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^


