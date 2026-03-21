LIBS     =  -stdlib=libc++

CFLAGS   =  -std=c++20 -O3 -ffast-math -fno-finite-math-only

CC       =   clang++

OBJECTS  =  Alignment.o AlignmentDistribution.o Analysis.o AnalysisComparison.o BitSet.o Container.o Exchangeabilities.o FileManager.o IndelRates.o main.o McmcSummary.o Msg.o Node.o ParameterStatistics.o Partition.o RandomVariable.o RateMatrix.o RateMatrixAverage.o RateMatrixF81.o RateMatrixGTR.o RateMatrixHelper.o RateMatrixJC69.o RateMatrixNaturalClass.o Sequence.o Simplex.o StateFrequencies.o Statistics.o Subset.o Threads.o Tree.o TreeSamples.o UserSettings.o

PROGS    = ttread

all:		$(PROGS)

ttread:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o ttread
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Alignment.o:	Alignment.cpp
		$(CC) $(CFLAGS) -c Alignment.cpp

AlignmentDistribution.o:	AlignmentDistribution.cpp
		$(CC) $(CFLAGS) -c AlignmentDistribution.cpp
		
Container.o:	Container.cpp
		$(CC) $(CFLAGS) -c Container.cpp
		
McmcSummary.o:	McmcSummary.cpp
		$(CC) $(CFLAGS) -c McmcSummary.cpp
		
Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Node.o:	Node.cpp
		$(CC) $(CFLAGS) -c Node.cpp

ParameterStatistics.o:	ParameterStatistics.cpp
		$(CC) $(CFLAGS) -c ParameterStatistics.cpp

Partition.o:	Partition.cpp
		$(CC) $(CFLAGS) -c Partition.cpp

Analysis.o:	Analysis.cpp
		$(CC) $(CFLAGS) -c Analysis.cpp

AnalysisComparison.o:	AnalysisComparison.cpp
		$(CC) $(CFLAGS) -c AnalysisComparison.cpp

BitSet.o:	BitSet.cpp
		$(CC) $(CFLAGS) -c BitSet.cpp

Sequence.o:	Sequence.cpp
		$(CC) $(CFLAGS) -c Sequence.cpp

Subset.o:	Subset.cpp
		$(CC) $(CFLAGS) -c Subset.cpp

Tree.o:	Tree.cpp
		$(CC) $(CFLAGS) -c Tree.cpp

UserSettings.o:	UserSettings.cpp
		$(CC) $(CFLAGS) -c UserSettings.cpp

Exchangeabilities.o:	Exchangeabilities.cpp
		$(CC) $(CFLAGS) -c Exchangeabilities.cpp

FileManager.o:	FileManager.cpp
		$(CC) $(CFLAGS) -c FileManager.cpp

IndelRates.o:	IndelRates.cpp
		$(CC) $(CFLAGS) -c IndelRates.cpp

RateMatrix.o:	RateMatrix.cpp
		$(CC) $(CFLAGS) -c RateMatrix.cpp

RateMatrixAverage.o:	RateMatrixAverage.cpp
		$(CC) $(CFLAGS) -c RateMatrixAverage.cpp

RateMatrixF81.o:	RateMatrixF81.cpp
		$(CC) $(CFLAGS) -c RateMatrixF81.cpp

RateMatrixGTR.o:	RateMatrixGTR.cpp
		$(CC) $(CFLAGS) -c RateMatrixGTR.cpp

RateMatrixHelper.o:	RateMatrixHelper.cpp
		$(CC) $(CFLAGS) -c RateMatrixHelper.cpp

RateMatrixJC69.o:	RateMatrixJC69.cpp
		$(CC) $(CFLAGS) -c RateMatrixJC69.cpp

RateMatrixNaturalClass.o:	RateMatrixNaturalClass.cpp
		$(CC) $(CFLAGS) -c RateMatrixNaturalClass.cpp

Simplex.o:	Simplex.cpp
		$(CC) $(CFLAGS) -c Simplex.cpp

StateFrequencies.o:	StateFrequencies.cpp
		$(CC) $(CFLAGS) -c StateFrequencies.cpp

Statistics.o:	Statistics.cpp
		$(CC) $(CFLAGS) -c Statistics.cpp

Threads.o:	Threads.cpp
		$(CC) $(CFLAGS) -c Threads.cpp

TreeSamples.o:	TreeSamples.cpp
		$(CC) $(CFLAGS) -c TreeSamples.cpp

clean:		
		rm -f *.o


