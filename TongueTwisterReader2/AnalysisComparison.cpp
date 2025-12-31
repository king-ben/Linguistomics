#include <unordered_set>
#include "AlignmentDistribution.hpp"
#include "Analysis.hpp"
#include "AnalysisComparison.hpp"
#include "Container.hpp"
#include "Msg.hpp"



void AnalysisComparison::compare(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples) {

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

std::pair<double,double> AnalysisComparison::compareQ(Analysis* a1, Analysis* a2, size_t numStates, size_t nSamples) {

    std::vector<float> f1;
    std::vector<float> f2;
    double sumW = 0.0, sumQ = 0.0;
    
    for (size_t n = 0; n < nSamples; n++)
        {
        DoubleMatrix* m1 = a1->randomlyChooseRateMatrixAndFreqs(f1);
        DoubleMatrix* m2 = a2->randomlyChooseRateMatrixAndFreqs(f2);
        numStates = m1->getNumRows();
        
        double ss = 0.0;
        for (size_t i=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                double x = f1[i] * (*m1)(i,j) - f2[i] * (*m2)(i,j);
                double y = f1[j] * (*m1)(j,i) - f2[j] * (*m2)(j,i);
                ss += (x*x + y*y);
                }
            }
        sumW += ss;
        
        ss = 0.0;
        for (size_t i=0; i<numStates; i++)
            {
            for (size_t j=i+1; j<numStates; j++)
                {
                double x = (*m1)(i,j) - (*m2)(i,j);
                double y = (*m1)(j,i) - (*m2)(j,i);
                ss += (x*x + y*y);
                }
            }
        sumQ += ss;
        }
    double dW = sumW / nSamples;
    double dQ = sumQ / nSamples;
    std::pair<double,double> d = std::make_pair(dW,dQ);
    return d;
}
