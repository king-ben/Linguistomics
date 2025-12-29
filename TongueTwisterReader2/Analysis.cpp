#include "AlignmentDistribution.hpp"
#include "Analysis.hpp"
#include "Exchangeabilities.hpp"
#include "McmcSummary.hpp"
#include "Partition.hpp"
#include "RateMatrixF81.hpp"
#include "RateMatrixGTR.hpp"
#include "RateMatrixJC69.hpp"
#include "RateMatrixNaturalClass.hpp"
#include "StateFrequencies.hpp"
#include "Threads.hpp"



Analysis::Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction) : 
    rng(r), pool(tp), rates(nullptr), part(nullptr), Q(nullptr), freqs(nullptr) {

    McmcSummary summary(pool, directoryName);
        
    if (summary.hasPartition() == true)
        part = new Partition(*summary.getStatePartition());
        
    if (summary.hasFrequencies() == true)
        freqs = new StateFrequencies(summary, burnFraction);
        
    if (summary.hasExchangeabilities() == true)
        rates = new Exchangeabilities(summary, burnFraction);
    
    if (freqs != nullptr && rates == nullptr)
        {
        std::cout << "   * Constructing F81 rate matrix" << std::endl;
        Q = new RateMatrixF81(*freqs);
        }
    else if (freqs != nullptr && rates != nullptr)
        {
        size_t numStates = freqs->getNumStates();
        size_t numRates = rates->getNumRates();
        if (numStates*(numStates-1)/2 == numRates)
            {
            std::cout << "   * Constructing GTR rate matrix" << std::endl;
            Q = new RateMatrixGTR(*freqs, *rates);
            }
        else
            {
            std::cout << "   * Constructing Natural Class rate matrix" << std::endl;
            Q = new RateMatrixNaturalClass(*freqs, *rates);
            }
        }
    else 
        {
        std::cout << "   * Constructing JC69 rate matrix" << std::endl;
        Q = new RateMatrixJC69(summary.getNumStates());
        }
    
    for (size_t n=0; n<summary.numAlignments(); n++)
        {
        AlignmentDistribution* alnDist = new AlignmentDistribution(rng, summary.getAlignments(n), summary.getAlignmentName(n));
        alignments.push_back(alnDist);
        }
}

Analysis::~Analysis(void) {

    if (freqs != nullptr)
        delete freqs;
    if (rates != nullptr)
        delete rates;
    if (part != nullptr)
        delete part;
    if (Q != nullptr)
        delete Q;
    for (AlignmentDistribution* ad : alignments)
        delete ad;
}

void Analysis::print(void) {

    if (part)
        part->print();
        
    if (freqs)
        freqs->print(part);
        
    if (rates)
        rates->print();
        
    if (Q != nullptr)
        {
        DoubleMatrix& m = Q->getMean();
        m.print();
        }

    std::cout << "alignments.size()=" << alignments.size() << std::endl;
    for (AlignmentDistribution* ad : alignments)
        ad->print(false, part);
}

void Analysis::printSorted(void) {

    if (freqs != nullptr)
        {
        std::map<float,StatInfo> sortedFreqs = freqs->sortByKL();
        for (auto [key,val] : sortedFreqs)
            {
            std::string range = std::to_string(val.index1);
            std::cout << std::setw(6) << range << " ";
            std::cout << key << " ";
            std::cout << val.priorMean << " -> ";
            std::cout << val.mean << " (";
            std::cout << val.lowerCI << ",";
            std::cout << val.upperCI << ")";
            std::cout << std::endl;
            }
        }

    if (rates != nullptr)
        {
        std::map<float,StatInfo> sortedRates = rates->sortByKL();
        for (auto [key,val] : sortedRates)
            {
            std::string range = std::to_string(val.index1) + "-" + std::to_string(val.index2);
            std::cout << std::setw(6) << range << " ";
            std::cout << key << " ";
            std::cout << val.priorMean << " -> ";
            std::cout << val.mean << " (";
            std::cout << val.lowerCI << ",";
            std::cout << val.upperCI << ")";
            std::cout << std::endl;
            }
        }
}

void Analysis::writeNytril(std::string fn) {

}
