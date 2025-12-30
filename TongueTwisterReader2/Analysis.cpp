#include <sstream>
#include "AlignmentDistribution.hpp"
#include "Analysis.hpp"
#include "Exchangeabilities.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "RateMatrixF81.hpp"
#include "RateMatrixGTR.hpp"
#include "RateMatrixJC69.hpp"
#include "RateMatrixNaturalClass.hpp"
#include "StateFrequencies.hpp"
#include "Threads.hpp"
#include "Tree.hpp"
#include "TreeSamples.hpp"



Analysis::Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction) : 
    rng(r), pool(tp), rates(nullptr), part(nullptr), Q(nullptr), freqs(nullptr), trees(nullptr) {

    McmcSummary summary(pool, directoryName);
    
    numStates = summary.getNumStates();
        
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
            if (part == nullptr)
                Msg::error("Must have partition to construct the Natural Class model");
            Q = new RateMatrixNaturalClass(*freqs, *rates, part);
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
        
    if (summary.hasTrees() == true)
        trees = new TreeSamples(summary, burnFraction);
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
    if (trees != nullptr)
        delete trees;
}

std::string Analysis::modelName(void) {

    if (freqs != nullptr && rates == nullptr)
        {
        return "F81";
        }
    else if (freqs != nullptr && rates != nullptr)
        {
        size_t numStates = freqs->getNumStates();
        size_t numRates = rates->getNumRates();
        if (numStates*(numStates-1)/2 == numRates)
            {
            return "GTR";
            }
        else
            {
            if (part == nullptr)
                Msg::error("Must have partition to construct the Natural Class model");
            return "Natural Class";
            }
        }
    return "JC69";
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
        
    if (trees != nullptr)
        trees->print();
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

void Analysis::randomlyChooseFreqs(std::vector<float>& f) {

    if (freqs == nullptr)
        {
        double x = static_cast<float>(1.0) / f.size();
        for (size_t i=0; i<f.size(); i++)
            f[i] = x;
        }
    else 
        {
        freqs->valuesAtIndex((size_t)(rng->uniformRv()*freqs->numSamples()), f);
        }
}

DoubleMatrix* Analysis::randomlyChooseRateMatrix(void) {

    std::vector<DoubleMatrix*>& rateMatrices = Q->getRateMatrices();
    if (rateMatrices.size() == 0)
        return nullptr;
    return rateMatrices[(int)(rng->uniformRv()*rateMatrices.size())];
}

DoubleMatrix* Analysis::randomlyChooseRateMatrixAndFreqs(std::vector<float>& f) {

    f.resize(numStates); 

    std::vector<DoubleMatrix*>& rateMatrices = Q->getRateMatrices();
    if (rateMatrices.size() == 0)
        return nullptr;
    size_t numSamples = rateMatrices.size();
    
    int idx = (int)(rng->uniformRv() * numSamples);
    if (freqs == nullptr)
        {
        double x = static_cast<float>(1.0) / f.size();
        for (size_t i=0; i<numStates; i++)
            f[i] = x;
        }
    else 
        {
        if (freqs->numSamples() != numSamples)
            Msg::error("Cannot match up rate matrices and stationary frequencies");
        freqs->valuesAtIndex(idx, f);
        }
    
    return rateMatrices[idx];
}

void Analysis::writeNytril(std::string pathName) {

    // output the full set of alignments (and analyses) to the nytril file
    auto file = new std::ofstream(pathName + "/consensus.tre", std::ios::out);
    ConsensusTree* ct = trees->getConsensusTree();
    std::string newickStr = ct->getNewick(5);
    *file << newickStr << std::endl;
    file->close();
    delete file;

    file = new std::ofstream(pathName + "/alignments_full.nytril", std::ios::out);
    nytrilOutput(*file, 1000); // Sanity limit
    file->close();
    delete file;
 
    file = new std::ofstream(pathName + "/alignments.nytril", std::ios::out);
    nytrilOutput(*file, 20);   // Display limit for papers and presentations
    file->close();
    delete file;
}

void Analysis::nytrilOutput(std::ofstream& file, int maxAlignment) {

    double cutoff = 0.95;

    file << "// This file is auto-generated by TongueTwisterReader. Do not edit\n\n";
    file << "AlignmentCutoff = " << cutoff << ";\n\n";

    // create a json object
    auto json = nlohmann::json::object();
    
    // output state partition information
    if (part != nullptr)
        json["state_part"] = part->toJson();

    // output the consensus tree
    ConsensusTree* conTree = trees->getConsensusTree();
    json["consensus"]["tree"]   = conTree->getNewick(4);
    json["consensus"]["n_taxa"] = conTree->getNumTaxa();
    
#   if 0
    // output information on mean and credible interval for all real-valued parameters
    auto jStats = nlohmann::json::array();
    for (int i=0; i<stats.size(); i++)
        {
        std::cout << stats[i] << std::endl;
        nlohmann::json cogStats = stats[i]->toJson();
        jStats.push_back(cogStats);
        }
    json["stats"] = jStats;
    
    for (int i=0; i< jStats.size(); i++)
        {
        std::string name = jStats[i]["cognate"];
        //std::cout << "Name = " << name << std::endl;
        if (name == "Insertion" || name == "Deletion")
            {
            file << name << "RateMean = " << jStats[i]["mean"] << ";\n";
            file << name << "RateLow = " << jStats[i]["lower"] << ";\n";
            file << name << "RateHigh = " << jStats[i]["upper"] << ";\n\n";
            }
        }

    int numRateClasses = inferNumberOfRates();
    std::cout << "numRateClasses = " << numRateClasses << std::endl;
    if (numRateClasses == 0)
        {
        // have the JC69 model
        }
    else if (numStates * (numStates-1) / 2 == numRateClasses)
        {
        // have the GTR model
        }
    else
        {
        // have the custom model
        int numSubsets = statePartitions->numSubsets();

        DoubleMatrix lo(numSubsets,numSubsets);
        DoubleMatrix md(numSubsets,numSubsets);
        DoubleMatrix hi(numSubsets,numSubsets);
        for (int i=0; i<numSubsets; i++)
            {
            for (int j=0; j<numSubsets; j++)
                {
                lo(i,j) = -1.0;
                md(i,j) = -1.0;
                hi(i,j) = -1.0;
                }
            }
        for (int i=0; i<stats.size(); i++)
            {
            ParameterStatistics* s = stats[i];
            if (s->getName()[0] == 'R')
                {
                int r1 = 0, r2 = 0;
                getRateElements(s->getName(), r1, r2);
                CredibleInterval interval = s->getCredibleInterval();
                double mdVal = s->getMean();
                double loVal = interval.lower;
                double hiVal = interval.upper;
                r1--;
                r2--;
                lo(r1,r2) = loVal;
                lo(r2,r1) = loVal;
                md(r1,r2) = mdVal;
                md(r2,r1) = mdVal;
                hi(r1,r2) = hiVal;
                hi(r2,r1) = hiVal;
                std::cout << s->getName() << " " << r1 << " " << r2 << std::endl;
                }
            }
        for (int i=0; i<numSubsets; i++)
            {
            for (int j=0; j<numSubsets; j++)
                {
                if (md(i,j) < 0.0)
                    Msg::error("Failed to initialize all partition rates");
                }
            }
            
        file << "StatClass[][] TransitionStats = [\n";
        for (int i = 0; i < numSubsets; ++i)
            {
            if (i)
                file << ",\n";
            file << "  [";
            for (int j = 0; j < numSubsets; ++j)
                {
                if (j)
                    file << ",";
                file << "new(" << lo(i,j) << "," << md(i,j) << "," << hi(i,j) << ")";
                }
            file << "]";
            }
        file << "\n];\n\n";
        
        
        // print out the information for the subsets
        std::vector<CredibleInterval> ncq = partitionRates();
        file << "StatClass[][] NCQRates = [\n";
        for (int i = 0; i < numSubsets; ++i)
            {
            if (i)
                file << ",\n";
            file << "  [";
            for (int j = 0; j < numSubsets; ++j)
                {
                if (j)
                    file << ",";
                CredibleInterval& ci = ncq[i * numSubsets + j];
                file << "new(" << ci.lower << "," << ci.median << "," << ci.upper << ")";
                }
            file << "]";
            }
        file << "\n];\n\n";
        
        std::vector<CredibleInterval> ncf = partitionFrequencies();
        file << "StatClass[] NCQFreqs = [";
        for (int i = 0; i < numSubsets; ++i)
            {
            if (i)
                file << ",\n";
            CredibleInterval& ci = ncf[i];
            file << "new(" << ci.lower << "," << ci.median << "," << ci.upper << ")";
            }
        file << "];\n\n";
        }
    
    // output average rates of change
    // No credible intervals on this information.
    int ns = inferNumberOfStates();
    if (numStates != ns && ns != 0)
        Msg::error("Inconsistency in the number of states");
    DoubleMatrix aveRates(numStates,numStates);
    calculateAverageRates(aveRates);
    writeMatrix(file, aveRates, "AverageRates");


    DoubleMatrix q(numStates,numStates);
    calculateRates(q);
    writeMatrix(file, q, "QRates");
        
    int numSubsets = statePartitions->numSubsets();
    DoubleMatrix sq(numSubsets,numSubsets);
    double maxValue = 0.0;
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                int jss = statePartitions->indexOfSubsetWithValue(j) - 1;
                std::cout << iss << " -> " << jss << std::endl;
                //double rate = q(i,j);
                double rate = aveRates(i,j);
                sq(iss,jss) += rate;
                if (sq(iss,iss) > maxValue)
                    maxValue = sq(iss,iss);
                }
            }
        }
    // temp, for printing circles
    for (int i=0; i<numSubsets; i++)
        {
        for (int j=0; j<numSubsets; j++)
            {
            double r = sq(i,j)/maxValue;
            double d = sqrt(5184.0 * r);
            std::cout << std::fixed << std::setprecision(0) << d << " ";
            }
        std::cout << std::endl;
        }
    std::cout << std::setprecision(6);

    sq.print();
    
    std::vector<double> f;
    calculateFrequencies(f);
    std::vector<double> catFreqs(numSubsets);
    for (int i=0; i<numStates; i++)
        {
        int iss = statePartitions->indexOfSubsetWithValue(i) - 1;
        catFreqs[iss] += f[i];
        }
    for (int i=0; i<numSubsets; i++)
        std::cout << i << " -- " << catFreqs[i] << std::endl;
    
    // output information on the gappiness of the alignment
    int nt1 = alignments[0]->getMapAlignment()->getNumTaxa()+1;
    std::vector<std::vector<double>> gapSpectrumCnt;
    gapSpectrumCnt.resize(3);
    for (int i=0; i<3; i++)
        {
        gapSpectrumCnt[i].resize(nt1);
        for (int j=0; j<nt1; j++)
            gapSpectrumCnt[i][j] = 0.0;
        }
    for (int n=0; n<alignments.size(); n++)
        {
        AlignmentDistribution* cogAlignments = alignments[n];
        Alignment* aln = cogAlignments->getMapAlignment();
        std::vector<double> spectrum = aln->getGapSpectrum();
        if (spectrum.size() == 2)
            {
            gapSpectrumCnt[0][(int)spectrum[0]]++;
            gapSpectrumCnt[1][(int)spectrum[1]]++;
            }
        else if (spectrum.size() == 3)
            {
            gapSpectrumCnt[0][(int)spectrum[0]]++;
            gapSpectrumCnt[1][(int)(round(spectrum[1]))]++;
            gapSpectrumCnt[2][(int)spectrum[2]]++;
            }
        }
    file << "GapSpectrum = [\n";
    for (int i=0; i<3; i++)
        {
        file << "[";
        for (int j=0; j<nt1; j++)
            {
            if (j != 0)
                file << ",";
            file << gapSpectrumCnt[i][j];
            }
        file << "]";
        if (i != 2)
            file << ",";
        file << "\n";
        }
    file << "];\n\n";

    file << "AlignIndexClass[] Alignments = [\n";

        
    // output the credible set for all alignments
    auto sampledAlignments = nlohmann::json::array();
    for (int i=0; i<alignments.size(); i++)
        {
        auto& a = alignments[i];
        
        if (a->size() > 0)
            {
            auto cogAlns = a->toJson(cutoff, maxAlignment, file);
            sampledAlignments.push_back(cogAlns);
            }
        }
    json["cog_alns"] = sampledAlignments;

    file << "];\n\n";

    // output the pooled partition subset frequencies (if a subset is available)
    if (statePartitions != NULL && ns != 0)
        {
        std::map<int,double> subsetFreqs;
        
        for (int i=0; i<stats.size(); i++)
            {
            std::string n = stats[i]->getName();
            if (n[0] == 'F' && n[1] == '[')
                {
                int x = parseNumberFromFreqHeader(n);
                double val = stats[i]->getMean();
                
                //std::cout << n << " " << x << std::endl;
                Subset* s = statePartitions->findSubsetWithValue(x);
                int subsetIndex = s->getIndex();
                if (s == NULL)
                    Msg::error("Could not find subset with value " + std::to_string(x));
                
                std::map<int,double>::iterator it = subsetFreqs.find(subsetIndex);
                if (it == subsetFreqs.end())
                    {
                    subsetFreqs.insert( std::make_pair(subsetIndex, val) );
                    }
                else
                    {
                    it->second += val;
                    }
                
                }
            }
            double sum = 0.0;
            for (std::map<int,double>::iterator it = subsetFreqs.begin(); it != subsetFreqs.end(); it++)
                {
                std::cout << it->first << " " << it->second << std::endl;
                sum += it->second;
                }
            std::cout << "sum = " << sum << std::endl;
        json["part_freqs"] = statePartitions->toFile(subsetFreqs, file);
        }
#   endif
}

