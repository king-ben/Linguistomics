#include <unordered_set>
#include "AlignmentDistribution.hpp"
#include "Analysis.hpp"
#include "AnalysisComparison.hpp"
#include "Container.hpp"
#include "Msg.hpp"
#include "RateMatrix.hpp"



Distances AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples) {

    // compare Q and W
    std::pair<double,double> d = compareQ(a1, a2, numStates, nSamples);
    double dW = d.first;
    double dQ = d.second;
        
    // compare alignments
    double dA = compareAlignments(a1, a2);

    std::cout << "Model Comparison:           " << a1->modelName() << "-" << a2->modelName() << std::endl;
    std::cout << "Average W distance:         " << dW << std::endl;
    std::cout << "Average W distance:         " << sqrt(dW) << std::endl;
    std::cout << "Average Q distance:         " << dQ << std::endl;
    std::cout << "Average Q distance:         " << sqrt(dQ) << std::endl;
    std::cout << "Average alignment distance: " << dA << std::endl;
    std::cout << "Average alignment distance: " << sqrt(dA) << std::endl;
    std::cout << std::endl;
    
    Distances dist;
    dist.qD = dQ;
    dist.wD = dW;
    dist.aD = dA;
    return dist;
}

double AnalysisComparison::compareAlignments(Analysis* a1, Analysis* a2) {

    const std::vector<AlignmentDistribution*>& alns1 = a1->getAlignments();
    const std::vector<AlignmentDistribution*>& alns2 = a2->getAlignments();
    if (alns1.size() != alns2.size())
        Msg::error("Comparing alingment distributions with different number of aligments");
    size_t numAlignments = alns1.size();
    
    double sumSquares = 0.0;
    for (size_t n=0; n<numAlignments; n++)
        {
        AlignmentDistribution* alnDist1 = alns1[n];
        AlignmentDistribution* alnDist2 = alns2[n];
//        if (alnDist1->getName() != alnDist2->getName())
//            {
//            std::cout << alnDist1->getName() << " != " << alnDist2->getName() << std::endl;
//            Msg::error("Comparing different cognate distributions");
//            }
        const std::map<Alignment*,int>& post1 = alnDist1->getAlignmentSamples();
        const std::map<Alignment*,int>& post2 = alnDist2->getAlignmentSamples();
        int numSamples1 = alnDist1->numSamples();
        int numSamples2 = alnDist2->numSamples();

        HashAlignment hasher;
        std::unordered_map<size_t, double> alnList1;
        for (const auto& [aln, count] : post1)
            {
            size_t h = hasher(aln);
            double prob = static_cast<double>(count) / numSamples1;
            alnList1[h] = prob;
            }
        std::unordered_map<size_t, double> alnList2;
        for (const auto& [aln, count] : post2)
            {
            size_t h = hasher(aln);
            double prob = static_cast<double>(count) / numSamples2;
            alnList2[h] = prob;
            }

        std::unordered_set<size_t> uniqueAlignments;
        double ss = 0.0;
        for (const auto& [aln,prob] : alnList1)
            {
            if (uniqueAlignments.contains(aln) == false)
                {
                // have not compared this alignment
                double prob1 = prob;
                double prob2 = 0.0;
                if (alnList2.contains(aln) == true)
                    prob2 = alnList2[aln];
                ss += (prob1 - prob2) * (prob1 - prob2);
                uniqueAlignments.insert(aln);
                }
            }
         for (const auto& [aln,prob] : alnList2)
            {
            if (uniqueAlignments.contains(aln) == false)
                {
                // have not compared this alignment
                double prob2 = prob;
                double prob1 = 0.0;
                if (alnList1.contains(aln) == true)
                    prob1 = alnList1[aln];
                ss += (prob1 - prob2) * (prob1 - prob2);
                uniqueAlignments.insert(aln);
                }
            }
        sumSquares += ss;
        }
        
    return sumSquares;
}

void AnalysisComparison::compareAnalyses(std::string fileName, std::vector<Analysis*>& analyses) {

    DoubleMatrix dRateMatrix(analyses.size(),analyses.size());
    DoubleMatrix dWeightMatrix(analyses.size(),analyses.size());
    DoubleMatrix dAlignment(analyses.size(),analyses.size());
    
    if (analyses.size() > 1)
        {
        // the JC69 model has idential rate matrices, so sensible to look
        // at only one of that model
        int numJC69 = 0;
        int excludeIdx = -1;
        for (int i=0; i<analyses.size(); i++)
            {
            if (analyses[i]->modelName() == "JC69")
                {
                numJC69++;
                if (numJC69 == 2)
                    {
                    excludeIdx = i;
                    break;
                    }
                }
            }
            
        // compare analyses
        for (size_t i=0; i<analyses.size(); i++)
            {
            Analysis* a1 = analyses[i];
            for (size_t j=i+1; j<analyses.size(); j++)
                {
                Analysis* a2 = analyses[j];
                std::cout << i+1 << " - " << j+1 << std::endl;
                Distances d = AnalysisComparison::compare(a1, a2, a2->getNumStates(), 100);
                dRateMatrix(i,j)   = std::sqrt(d.qD);
                dRateMatrix(j,i)   = std::sqrt(d.qD);
                dWeightMatrix(i,j) = std::sqrt(d.wD);
                dWeightMatrix(j,i) = std::sqrt(d.wD);
                dAlignment(i,j)    = std::sqrt(d.aD);
                dAlignment(j,i)    = std::sqrt(d.aD);
                }
            }
            
        // print W distance matrix
        std::string fn = fileName + "distances_w.csv";
        std::ofstream file(fn, std::ios::out);
        bool first = true;
        for (size_t i=0; i<analyses.size(); i++)
            {
            if (i == excludeIdx)
                continue;
            Analysis* a1 = analyses[i];
            if (first == false)
                file << ", ";
            first = false;
            file << a1->getName();
            }
        file << std::endl;
        
        for (size_t i=0; i<analyses.size(); i++)
            {
            if (i == excludeIdx)
                continue;
            Analysis* a1 = analyses[i];
            file << a1->getName() << ", ";
            first = true;
            for (size_t j=0; j<analyses.size(); j++)
                {
                if (j == excludeIdx)
                    continue;
                if (first == false)
                    file << ", ";
                first = false;
                if (i == j)
                    file << 0.0;
                else
                    file << dWeightMatrix(i,j);
                }
            file << std::endl;
            }
        file << std::endl;
        file.close();
           
        // print alignment distance matrix
        fn = fileName + "distances_aln.csv";
        file.open(fn, std::ios::out);

        for (size_t i=0; i<analyses.size(); i++)
            {
            Analysis* a1 = analyses[i];
            if (i != 0)
                file << ", ";
            file << a1->getName();
            }
        file << std::endl;
        for (size_t i=0; i<analyses.size(); i++)
            {
            Analysis* a1 = analyses[i];
            file << a1->getName() << ", ";
            for (size_t j=0; j<analyses.size(); j++)
                {
                if (j != 0)
                    file << ", ";
                if (i == j)
                    file<< 0.0;
                else
                    file << dAlignment(i,j);
                }
            file << std::endl;
            }

        file.close();
        }
}


std::pair<double,double> AnalysisComparison::compareQ(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples) {

    RateMatrix* Q1 = a1->getRateMatrix();
    RateMatrix* Q2 = a2->getRateMatrix();
    RateMatrix* W1 = a1->getWeightMatrix();
    RateMatrix* W2 = a2->getWeightMatrix();

    DoubleMatrix& meanQ1 = Q1->getMean();
    DoubleMatrix& meanQ2 = Q2->getMean();
    DoubleMatrix& meanW1 = W1->getMean();
    DoubleMatrix& meanW2 = W2->getMean();
    
    double sumW = 0.0;
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            double x = meanW1(i,j) - meanW2(i,j);
            double y = meanW1(j,i) - meanW2(j,i);
            sumW += (x*x + y*y);
            }
        }
        
    double sumQ = 0.0;
    for (size_t i=0; i<numStates; i++)
        {
        for (size_t j=i+1; j<numStates; j++)
            {
            double x = meanQ1(i,j) - meanQ2(i,j);
            double y = meanQ1(j,i) - meanQ2(j,i);
            sumQ += (x*x + y*y);
            }
        }

    std::pair<double,double> d = std::make_pair(sumW,sumQ);
    return d;
}
