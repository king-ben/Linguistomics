#include <sstream>
#include "AlignmentDistribution.hpp"
#include "Analysis.hpp"
#include "Exchangeabilities.hpp"
#include "FileManager.hpp"
#include "IndelRates.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "RateMatrixAverage.hpp"
#include "RateMatrixF81.hpp"
#include "RateMatrixGTR.hpp"
#include "RateMatrixJC69.hpp"
#include "RateMatrixNaturalClass.hpp"
#include "StateFrequencies.hpp"
#include "Statistics.hpp"
#include "Threads.hpp"
#include "Tree.hpp"
#include "TreeSamples.hpp"



Analysis::Analysis(RandomVariable* r, ThreadPool* tp, std::string directoryName, double burnFraction) : 
    rng(r), pool(tp), indelRates(nullptr), rates(nullptr), part(nullptr), Q(nullptr), 
    aveQ(nullptr), freqs(nullptr), ncFreqs(nullptr), trees(nullptr) {

    // read the MCMC output files in the directory directoryName
    McmcSummary summary(pool, directoryName);
    
    // set the name
    analysisName = summary.getName();
    
    // initialize the number of model states
    numStates = summary.getNumStates();
    
    // post-burn insertion and deletion rate samples
    if (summary.hasIndelRates() == true)
        indelRates = new IndelRates(summary, burnFraction);
        
    // instantiate the partition object, if the definitiion exists in the config file
    if (summary.hasPartition() == true)
        part = new Partition(*summary.getStatePartition());
        
    // post-burn stationary frequencies of states (not present for JC69)
    if (summary.hasFrequencies() == true)
        freqs = new StateFrequencies(summary, burnFraction);
        
    // post-burn stationary frequencies of subset categories of partition
    if (freqs != nullptr && part != nullptr)
        ncFreqs = new StateFrequencies(*freqs, part);
        
    // post-burn rate parameters for GTR and Natural Class models
    if (summary.hasExchangeabilities() == true)
        rates = new Exchangeabilities(summary, burnFraction);
    
    // post-burn rate matrices (q_ij) and average rate (f_i X q_ij)
    if (freqs != nullptr && rates == nullptr)
        {
        std::cout << "   * Constructing F81 rate matrix" << std::endl;
        Q = new RateMatrixF81(*freqs);
        aveQ = new RateMatrixAverage(*Q, *freqs);
        }
    else if (freqs != nullptr && rates != nullptr)
        {
        size_t numStates = freqs->getNumStates();
        size_t numRates = rates->getNumRates();
        if (numStates*(numStates-1)/2 == numRates)
            {
            std::cout << "   * Constructing GTR rate matrix" << std::endl;
            Q = new RateMatrixGTR(*freqs, *rates);
            aveQ = new RateMatrixAverage(*Q, *freqs);
           }
        else
            {
            std::cout << "   * Constructing Natural Class rate matrix" << std::endl;
            if (part == nullptr)
                Msg::error("Must have partition to construct the Natural Class model");
            Q = new RateMatrixNaturalClass(*freqs, *rates, part);
            aveQ = new RateMatrixAverage(*Q, *freqs);
            }
        }
    else 
        {
        std::cout << "   * Constructing JC69 rate matrix" << std::endl;
        Q = new RateMatrixJC69(summary.getNumStates());
        aveQ = new RateMatrixAverage(summary.getNumStates());
        }
    
    // post-burn alignment samples
    for (size_t n=0; n<summary.numAlignments(); n++)
        {
        AlignmentDistribution* alnDist = new AlignmentDistribution(rng, summary.getAlignments(n), summary.getAlignmentName(n));
        alignments.push_back(alnDist);
        }
        
    //post-burn tree samples
    if (summary.hasTrees() == true)
        trees = new TreeSamples(summary, burnFraction);
}

Analysis::~Analysis(void) {

    if (indelRates != nullptr)
        delete indelRates;
    if (freqs != nullptr)
        delete freqs;
    if (ncFreqs != nullptr)
        delete ncFreqs;
    if (rates != nullptr)
        delete rates;
    if (part != nullptr)
        delete part;
    if (Q != nullptr)
        delete Q;
    if (aveQ != nullptr)
        delete aveQ;
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
        freqs->print();
        
    if (ncFreqs)
        ncFreqs->print();
        
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

    if (FileManager::ensureDirectoryExists(pathName) == true)
        std::cout << "   Creating directory for Nytril output" << std::endl;

    std::string fullNytrilFn = pathName + "/alignments_full.nytril";
    file = new std::ofstream(fullNytrilFn, std::ios::out);
    nytrilOutput(*file, 1000); // Sanity limit
    file->close();
    delete file;
 
    std::string smallNytrilFn = pathName + "/alignments.nytril";
    file = new std::ofstream(smallNytrilFn, std::ios::out);
    nytrilOutput(*file, 20);   // Display limit for papers and presentations
    file->close();
    delete file;
}

void Analysis::writeMatrix(std::ofstream& file, DoubleMatrix &m, std::string name) {

    file << name << " = [\n";
    auto rows = m.getNumRows();
    auto cols = m.getNumCols();
    double max = 0;
    for (int i = 0; i < rows; ++i)
        {
        if (i)
            file << ",\n";
        file << "  [";
        for (int j = 0; j < cols; ++j)
            {
            double r = m(i, j);
            if (r > max)
                max = r;
            if (j)
                file << ",";
            file << r;
            }
        file << "]";
        }
    file << "\n];\n";
    file << name << "Max = " << max << ";\n\n";
}

void Analysis::nytrilOutput(std::ofstream& file, int maxAlignment) {

    // Note: this is a horrible function
    
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
    
    // output information on mean and credible interval for all real-valued parameters
    auto jStats = nlohmann::json::array();
    if (indelRates)
        {
        nlohmann::json parmStats = indelRates->toJson();
        jStats.push_back(parmStats);
        }
    if (freqs)
        {
        nlohmann::json parmStats = freqs->toJson();
        jStats.push_back(parmStats);
        }
    if (rates)
        {
        nlohmann::json parmStats = rates->toJson();
        jStats.push_back(parmStats);
        }
    json["stats"] = jStats;
    
    file << "Insertion" << "RateMean = " << indelRates->getMean().first << ";\n";
    file << "Insertion" << "RateLow = " << indelRates->getLowerCI().first << ";\n";
    file << "Insertion" << "RateHigh = " << indelRates->getUpperCI().first << ";\n\n";

    file << "Deletion" << "RateMean = " << indelRates->getMean().second << ";\n";
    file << "Deletion" << "RateLow = " << indelRates->getLowerCI().second << ";\n";
    file << "Deletion" << "RateHigh = " << indelRates->getUpperCI().second << ";\n\n";

    if (auto* nc = dynamic_cast<RateMatrixNaturalClass*>(rates))
        {
        // have the custom/natural class model
        int numSubsets = part->numSubsets();
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

        for (size_t i=0; i<rates->getNumRates(); i++)
            {
            size_t r1 = rates->getFirstIndex(i);
            size_t r2 = rates->getSecondIndex(i);
            lo(r1,r2) = rates->getLowerCI(i);
            lo(r2,r1) = rates->getLowerCI(i);
            md(r1,r2) = rates->getMean(i);
            md(r2,r1) = rates->getMean(i);
            hi(r1,r2) = rates->getUpperCI(i);
            hi(r2,r1) = rates->getUpperCI(i);
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
        file << "StatClass[][] NCQRates = [\n";
        for (int i = 0; i < rates->size(); ++i)
            {
            if (i)
                file << ",\n";
            file << "  [";
            for (int j = 0; j < numSubsets; ++j)
                {
                if (j)
                    file << ",";
                file << "new(" << rates->getLowerCI(i) << "," << rates->getMean(i) << "," << rates->getUpperCI(i) << ")";
                }
            file << "]";
            }
        file << "\n];\n\n";
        
        if (ncFreqs)
            {
            file << "StatClass[] NCQFreqs = [";
            for (int i = 0; i < ncFreqs->size(); ++i)
                {
                if (i)
                    file << ",\n";
                file << "new(" << ncFreqs->getLowerCI(i) << "," << ncFreqs->getMean(i) << "," << ncFreqs->getUpperCI(i) << ")";
                }
            file << "];\n\n";
            }
        }
   
    // output average rates of change
    // No credible intervals on this information.
    // Shawn: We now can calculate credible intervals for the average rates. These
    //        values are obtained using aveQ->getLowerCI() and aveQ->getUpperCI()
    //        which both return a DoubleMatrix& where each element contains the value
    writeMatrix(file, aveQ->getMean(), "AverageRates");
    writeMatrix(file, Q->getMean(), "QRates");
    
    if (part && freqs)
        {
        // Note: I'm commenting this section out because it only appears to print
        //       to the screen. Moreover, these quantities are not calculated on 
        //       construction of this object.
#       if 0
        DoubleMatrix& averageQ = aveQ->getMean();
        int numSubsets = part->numSubsets();
        DoubleMatrix sq(numSubsets,numSubsets);
        double maxValue = 0.0;
        for (int i=0; i<numStates; i++)
            {
            int iss = part->indexOfSubsetWithValue(i) - 1;
            for (int j=0; j<numStates; j++)
                {
                if (i != j)
                    {
                    int jss = part->indexOfSubsetWithValue(j) - 1;
                    std::cout << iss << " -> " << jss << std::endl;
                    //double rate = q(i,j);
                    double rate = averageQ(i,j);
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
    
        std::vector<double> catFreqs(numSubsets);
        for (int i=0; i<numSubsets; i++)
            catFreqs[i] = ncFreqs->getMean(i);
        for (int i=0; i<numSubsets; i++)
            std::cout << i << " -- " << catFreqs[i] << std::endl;
#       endif
        }
    
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
    if (part != nullptr && freqs != nullptr)
        {
        std::map<int,double> subsetFreqs;
        for (int i=0; i<ncFreqs->size(); i++)
            {
            int subsetIndex = i+1;
            double val = ncFreqs->getMean(i);
            subsetFreqs.insert( std::make_pair(subsetIndex, val) );
            }
        json["part_freqs"] = part->toFile(subsetFreqs, file);
        }
}

void Analysis::writeR(std::string pathName, std::string fileName, std::string modelName, int idx) {

    // output the average rate matrix
    if (aveQ != nullptr)
        {
        std::string fn = pathName + "/" + analysisName + "_w_matrix_" + modelName + "_" + std::to_string(idx) + ".csv";
        std::ofstream file(fn, std::ios::out);

        DoubleMatrix& W = aveQ->getMean();
        std::string matrixId = modelName + "_" + std::to_string(idx);
        size_t nStates = W.getNumRows();
        file << "matrix_id,from_state,to_state,rate" << std::endl;
        for (size_t i=0; i<nStates; i++)
            {
            for (size_t j=0; j<nStates; j++)
                {
                file << "W_" << matrixId << "," << std::to_string(i) << "," << std::to_string(j) << "," << W(i,j) << std::endl;
                }
            }
    
        file.close();
        }

    // output the stationary frequencies
    if (freqs != nullptr)
        {
        std::string fn = pathName + "/" + analysisName + "_freqs_" + modelName + "_" + std::to_string(idx) + ".csv";
        std::ofstream file(fn, std::ios::out);
        
        file << "param_id,param_group,param_index,value_mean,kl_to_prior,q025,q975" << std::endl;
        size_t nStates = freqs->getNumStates();
        for (size_t i=0; i<nStates; i++)
            {
            file << "pi[" << i << "],";
            file << "pi,";
            file << i << ",";
            file << freqs->getMean(i) << ",";
            file << freqs->getKullbackLiebler(i) << ",";
            file << freqs->getLowerCI(i) << ",";
            file << freqs->getUpperCI(i) << std::endl;
            }
        
        file.close();
        }

    // output th exchangeability parameters
    if (rates != nullptr)
        {
        std::string fn = pathName + "/" + analysisName + "_rates_" + modelName + "_" + std::to_string(idx) + ".csv";
        std::ofstream file(fn, std::ios::out);
        
        file << "param_id,param_group,param_index,param_index_from,param_index_to,value_mean,kl_to_prior,q025,q975" << std::endl;
        size_t nRates = rates->getNumRates();
        for (size_t i=0; i<nRates; i++)
            {
            size_t from = rates->getFirstIndex(i);
            size_t to = rates->getSecondIndex(i);
            file << "R[" << i << "],";
            file << "R,";
            file << i+1 << "," << from << "," << to << ",";
            file << rates->getMean(i) << ",";
            file << rates->getKullbackLiebler(i) << ",";
            file << rates->getLowerCI(i) << ",";
            file << rates->getUpperCI(i) << std::endl;
            }
        
        file.close();
        }
        
    // output alignment distributions
    {
    std::string fn = pathName + "/" + analysisName + "_alns_" + modelName + "_" + std::to_string(idx) + ".csv";
    std::ofstream file(fn, std::ios::out);

    file << "dist_id,align_id,rank,prob" << std::endl;
    HashAlignment hasher;
    for (size_t i=0; i<alignments.size(); i++)
        {
        AlignmentDistribution* alnDist = alignments[i];
        std::string alnName = alnDist->getName();
        std::replace(alnName.begin(), alnName.end(), ' ', '_');
        std::replace(alnName.begin(), alnName.end(), '-', '_');
        if (alnName.starts_with("Alignment_for_cognate_"))
            alnName.erase(0, std::string("Alignment_for_cognate_").size());
        std::map<int,std::pair<Alignment*,float>> sortedAlns = alnDist->getSortedAlignments();
        
        int alnNum = 0;
        for (auto& [key,val] : sortedAlns)
            {
            alnNum++;
            size_t h = hasher(val.first);
            file << alnName << ",";
            file << h << ",";
            file << key << ",";
            file << val.second << std::endl;
            }
        }
 
    file.close();
    }

}
