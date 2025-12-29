LIBS     =  -stdlib=libc++

CFLAGS   =  -Ofast -std=c++17

CC       =   clang++

OBJECTS  =  main.o Alignment.o AlignmentProposal.o AlnMatrix.o Container.o IndelMatrix.o IntVector.o LikelihoodCalculator.o MatrixMath.o Mcmc.o McmcPhase.o Model.o Msg.o Node.o NodeFactory.o NodeSet.o Parameter.o ParameterAlignment.o ParameterEquilibirumFrequencies.o ParameterExchangabilityRates.o ParameterIndelGammaShape.o ParameterIndelRates.o ParameterRatesGammaShape.o ParameterTree.o Partition.o Probability.o RandomVariable.o RateMatrix.o RateMatrixHelper.o RbBitSet.o SteppingStones.o Stone.o Subset.o TaxonPair.o Threads.o TransitionProbabilities.o Tree.o UpdateInfo.o UserSettings.o

PROGS    = tonguetwister

all:		$(PROGS)

tonguetwister:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o tonguetwister
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Alignment.o:	Alignment.cpp
		$(CC) $(CFLAGS) -c Alignment.cpp

AlignmentProposal.o:	AlignmentProposal.cpp
		$(CC) $(CFLAGS) -c AlignmentProposal.cpp
		
AlnMatrix.o:	AlnMatrix.cpp
		$(CC) $(CFLAGS) -c AlnMatrix.cpp
		
Container.o:	Container.cpp
		$(CC) $(CFLAGS) -c Container.cpp

IndelMatrix.o:	IndelMatrix.cpp
		$(CC) $(CFLAGS) -c IndelMatrix.cpp

IntVector.o:	IntVector.cpp
		$(CC) $(CFLAGS) -c IntVector.cpp

LikelihoodCalculator.o:	LikelihoodCalculator.cpp
		$(CC) $(CFLAGS) -c LikelihoodCalculator.cpp

MatrixMath.o:	MatrixMath.cpp
		$(CC) $(CFLAGS) -c MatrixMath.cpp

Mcmc.o:	Mcmc.cpp
		$(CC) $(CFLAGS) -c Mcmc.cpp

McmcPhase.o:	McmcPhase.cpp
		$(CC) $(CFLAGS) -c McmcPhase.cpp

Model.o:	Model.cpp
		$(CC) $(CFLAGS) -c Model.cpp

Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Node.o:	Node.cpp
		$(CC) $(CFLAGS) -c Node.cpp

NodeFactory.o:	NodeFactory.cpp
		$(CC) $(CFLAGS) -c NodeFactory.cpp

NodeSet.o:	NodeSet.cpp
		$(CC) $(CFLAGS) -c NodeSet.cpp

Parameter.o:	Parameter.cpp
		$(CC) $(CFLAGS) -c Parameter.cpp

ParameterAlignment.o:	ParameterAlignment.cpp
		$(CC) $(CFLAGS) -c ParameterAlignment.cpp

ParameterEquilibirumFrequencies.o:	ParameterEquilibirumFrequencies.cpp
		$(CC) $(CFLAGS) -c ParameterEquilibirumFrequencies.cpp

ParameterExchangabilityRates.o:	ParameterExchangabilityRates.cpp
		$(CC) $(CFLAGS) -c ParameterExchangabilityRates.cpp

ParameterRatesGammaShape.o:	ParameterRatesGammaShape.cpp
		$(CC) $(CFLAGS) -c ParameterRatesGammaShape.cpp

ParameterIndelGammaShape.o:	ParameterIndelGammaShape.cpp
		$(CC) $(CFLAGS) -c ParameterIndelGammaShape.cpp
		
ParameterIndelRates.o:	ParameterIndelRates.cpp
		$(CC) $(CFLAGS) -c ParameterIndelRates.cpp

ParameterTree.o:	ParameterTree.cpp
		$(CC) $(CFLAGS) -c ParameterTree.cpp

Partition.o:	Partition.cpp
		$(CC) $(CFLAGS) -c Partition.cpp

Probability.o:	Probability.cpp
		$(CC) $(CFLAGS) -c Probability.cpp

RandomVariable.o:	RandomVariable.cpp
		$(CC) $(CFLAGS) -c RandomVariable.cpp

RateMatrix.o:	RateMatrix.cpp
		$(CC) $(CFLAGS) -c RateMatrix.cpp

RateMatrixHelper.o:	RateMatrixHelper.cpp
		$(CC) $(CFLAGS) -c RateMatrixHelper.cpp

RbBitSet.o:	RbBitSet.cpp
		$(CC) $(CFLAGS) -c RbBitSet.cpp

SteppingStones.o:	SteppingStones.cpp
		$(CC) $(CFLAGS) -c SteppingStones.cpp

Stone.o:	Stone.cpp
		$(CC) $(CFLAGS) -c Stone.cpp

Subset.o:	Subset.cpp
		$(CC) $(CFLAGS) -c Subset.cpp

TaxonPair.o:	TaxonPair.cpp
		$(CC) $(CFLAGS) -c TaxonPair.cpp

Threads.o:	Threads.cpp
		$(CC) $(CFLAGS) -c Threads.cpp

TransitionProbabilities.o:	TransitionProbabilities.cpp
		$(CC) $(CFLAGS) -c TransitionProbabilities.cpp

Tree.o:	Tree.cpp
		$(CC) $(CFLAGS) -c Tree.cpp

UpdateInfo.o:	UpdateInfo.cpp
		$(CC) $(CFLAGS) -c UpdateInfo.cpp

UserSettings.o:	UserSettings.cpp
		$(CC) $(CFLAGS) -c UserSettings.cpp

clean:		
		rm -f *.o
