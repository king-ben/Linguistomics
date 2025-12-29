#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "IntVector.hpp"
#include "Mcmc.hpp"
#include "McmcPhase.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Parameter.hpp"
#include "ParameterAlignment.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "SteppingStones.h"
#include "Tree.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"




Mcmc::Mcmc(Model** m, RandomVariable* r) {

    modelPtr      = m;
    rv            = r;
    maxLikePrint  = 0;
    maxPriorPrint = 0;
    
    UserSettings& settings = UserSettings::userSettings();
    numChains              = settings.getNumChains();
    temperature            = settings.getTemperature();
    numMcmcCycles          = settings.getNumMcmcCycles();
    printFrequency         = settings.getPrintFrequency();
    sampleFrequency        = settings.getSampleFrequency();
    preburninLength        = settings.getPreburninLength();
    tuneLength             = settings.getTuneLength();
    burninLength           = settings.getBurnLength();
    sampleLength           = settings.getSampleLength();
    maxGenPrint            = numDigits(settings.getNumMcmcCycles());
}

std::vector<double> Mcmc::calculatePowers(int numStones, double alpha, double beta) {

    int ns = numStones - 1;
    std::vector<double> pwrs;
    double intervalProb = (double)1.0 / ns;
    pwrs.push_back(1.0);
    for (int i=ns-1; i>0; i--)
        pwrs.push_back( Probability::Beta::quantile(alpha, beta, i * intervalProb) );
    pwrs.push_back(0.0);
#   if 0
    for (int i=0; i<pwrs.size(); i++)
        std::cout << i+1 << " -- " << std::fixed << std::setprecision(15) << pwrs[i] << std::endl;
#   endif
    return pwrs;
}

void Mcmc::chooseModelsToSwap(int& idx1, int& idx2) {

    if (numChains == 1)
        Msg::error("Cannot swap when there is only one chain");
        
    idx1 = (int)(rv->uniformRv() * numChains);
    do {
        idx2 = (int)(rv->uniformRv() * numChains);
        } while (idx2 == idx1);
}

void Mcmc::closeOutputFiles(void) {

    parmStrm.close();
    treeStrm.close();
}

std::string Mcmc::formattedTime(Timer& t1, Timer& t2) {

    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    int s = (int)durationSecs.count();
    int m = s / 60;
    int h = s / 3600;
        
    std::string tStr = "";
    if (h > 0)
        {
        tStr += std::to_string(h) + "h:";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        tStr += std::to_string(m) + "m:";
        s -= m * 60;
        }
    tStr += std::to_string(s) + "s";
    
    return tStr;
}

Model* Mcmc::getColdModel(void) {

    for (int chain=0; chain<numChains; chain++)
        {
        if (modelPtr[chain]->getIndex() == 0)
            return modelPtr[chain];
        }
    return NULL;
}

void Mcmc::initialize(void) {

    openOutputFiles();
    numParmValues = modelPtr[0]->getNumParameterValues();
    parmValues = new double[numParmValues];
    for (int i=0; i<numParmValues; i++)
        parmValues[i] = 0.0;
}

int Mcmc::numDigits(double x) {

    if (x < 0.0)
        x = -x;
    return (int)(log(x) / log(10.0)) + 1;
}

void Mcmc::openOutputFiles(void) {

    // open files for samples
    UserSettings& settings = UserSettings::userSettings();
    std::string outPath = settings.getOutFile();
    std::string parmFileName = outPath + ".tsv";
    std::string treeFileName = outPath + ".tre";

    parmStrm.open( parmFileName.c_str(), std::ios::out );
    if (!parmStrm)
        Msg::error("Cannot open file \"" + parmFileName + "\"");
    treeStrm.open( treeFileName.c_str(), std::ios::out );
    if (!treeStrm)
        Msg::error("Cannot open file \"" + treeFileName + "\"");
}

double Mcmc::power(int idx) {

    return 1.0 / (1.0 + temperature * modelPtr[idx]->getIndex());
}

void Mcmc::print(int powIdx, int numPowers, std::string phase, int gen, double* curLnL, double* curLnP) {

    for (int chain=0; chain<numChains; chain++)
        {
        if (numDigits(curLnL[chain]) > maxLikePrint)
            maxLikePrint = numDigits(curLnL[chain]);
        if (numDigits(curLnP[chain]) > maxPriorPrint)
            maxPriorPrint = numDigits(curLnP[chain]);
        }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "   * ";
    std::cout << powIdx << "/" << numPowers << " " << phase << " ";
    std::cout << std::setw(maxGenPrint) << gen << " --   ";
    for (int chain=0; chain<numChains; chain++)
        {
        bool isCold = false;
        if (modelPtr[chain] == getColdModel())
            isCold = true;
        if (isCold == true)
            std::cout << "[";
        else
            std::cout << "(";
        std::cout << std::setw(maxLikePrint + 5) << curLnL[chain];
        if (isCold == true)
            std::cout << "] ";
        else
            std::cout << ") ";
        }

    std::cout << std::endl;

}

void Mcmc::print(int gen, double* curLnL, double* curLnP, bool /*accept*/, Timer& t1, Timer& t2) {

    for (int chain=0; chain<numChains; chain++)
        {
        if (numDigits(curLnL[chain]) > maxLikePrint)
            maxLikePrint = numDigits(curLnL[chain]);
        if (numDigits(curLnP[chain]) > maxPriorPrint)
            maxPriorPrint = numDigits(curLnP[chain]);
        }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "   * ";
    std::cout << std::setw(maxGenPrint) << gen << " --   ";
    for (int chain=0; chain<numChains; chain++)
        {
        bool isCold = false;
        if (modelPtr[chain] == getColdModel())
            isCold = true;
        if (isCold == true)
            std::cout << "[";
        else
            std::cout << "(";
        std::cout << std::setw(maxLikePrint + 5) << curLnL[chain];
        if (isCold == true)
            std::cout << "] ";
        else
            std::cout << ") ";
        }
    //std::cout << std::setw(maxLikePrint + 5) << curLnL << " -> " << std::setw(maxLikePrint + 5) << newLnL << "   ";

    // estimate time remaining based on time to this point
    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t1 - t2);
    double timePerCycle = (double)durationSecs.count() / gen;
    if (timePerCycle == 0)
        timePerCycle = 1.0 / printFrequency;
    int s = (int)((numMcmcCycles - gen) * timePerCycle);
    int m = s / 60;
    int h = s / 3600;

    if (h > 0)
        {
        std::cout << h;
        std::cout << ":";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        if (m < 10)
            std::cout << "0";
        std::cout << m;
        std::cout << ":";
        s -= m * 60;
        }
    if (s < 10) 
        std::cout << "0";
    std::cout << s;
    std::cout << " remaining  ";

    //std::cout << std::setw(maxPriorPrint + 5) << curLnP << " -> " << std::setw(maxPriorPrint + 5) << newLnP << "   ";
#   if 0
    if (accept == true)
        std::cout << "Accepted ";
    else
        std::cout << "Rejected ";
    std::cout << getColdModel()->getUpdatedParameterName();
    std::cout << " parameter";
    std::string& updateType = getColdModel()->getLastUpdate();
    if (updateType != "")
        {
        std::cout << " (";
        std::cout << updateType;
        std::cout << " update)";
        }
#   endif
    std::cout << std::endl;
}

void Mcmc::run(void) {

    UserSettings& settings = UserSettings::userSettings();
    if (settings.getCalculateMarginalLikelihood() == false)
        runPosterior();
    else
        runPathSampling();
}

void Mcmc::runPathSampling(void) {

    // path sampling only implemented for one chain
    if (numChains > 1)
        Msg::error("Cannot run path sampling with more than one chain");
        
    // get a reference to the user settings object
    UserSettings& settings = UserSettings::userSettings();

    // initialize the chain
    initialize();

    // get powers and initialize vector with phases
    int numStones = 127;
    double alpha  = 0.3;
    double beta   = 1.0;
    std::vector<double> powers = calculatePowers(numStones, alpha, beta);
    SteppingStones samples(powers);
    McmcPhase mcmcPhases(settings.getFirstBurnLength(), settings.getPreburninLength(), settings.getTuneLength(), settings.getBurnLength(), settings.getSampleLength(), settings.getSampleFrequency());
    std::vector<std::string>& phases = mcmcPhases.getPhases();
    int firstBurnLength = mcmcPhases.getFirstBurnLength();

    // get initial likelihoods
    for (int chain=0; chain<numChains; chain++)
        modelPtr[chain]->setUpdateLikelihood();
    double* curLnL = new double[numChains];
    double* curLnP = new double[numChains];
    for (int chain=0; chain<numChains; chain++)
        {
        curLnL[chain] = modelPtr[chain]->lnLikelihood();
        curLnP[chain] = modelPtr[chain]->lnPriorProbability();
        }
        
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();

    // Metropolis-Hastings algorithm
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int powIdx=0; powIdx<powers.size(); powIdx++)
        {
        double power = powers[powIdx];
        for (std::string phase : phases)
            {
            int extra = 0;
            if (powIdx == 0 && phase == "preburn")
                extra = firstBurnLength;
            for (int n=1; n<=mcmcPhases.getLengthForPhase(phase)+extra; n++)
                {
                // propose a new value for the chain
                double lnProposalRatio = modelPtr[0]->update(n);
                
                // calculate the likelihood and prior ratios (natural log scale)
                double newLnL = modelPtr[0]->lnLikelihood();
                double lnLikelihoodRatio = newLnL - curLnL[0];
                double newLnP = modelPtr[0]->lnPriorProbability();
                double lnPriorRatio = newLnP - curLnP[0];
                lnLikelihoodRatio *= power;

                // accept or reject the state
                bool accept = false;
                if ( log(rv->uniformRv()) < lnLikelihoodRatio + lnPriorRatio + lnProposalRatio )
                    accept = true;
                
                // print to the screen to give the user some clue where the chain is
                if (n == 1 || n == numMcmcCycles || n % printFrequency == 0)
                    print(powIdx, (int)powers.size(), phase, n, curLnL, curLnP);
                
                // update the state of the chain
                if (accept == false)
                    {
                    modelPtr[0]->reject();
                    updateInfo.reject(modelPtr[0]->getLastUpdate());
                    }
                else
                    {
                    modelPtr[0]->accept();
                    updateInfo.accept(modelPtr[0]->getLastUpdate());
                    curLnL[0] = newLnL;
                    curLnP[0] = newLnP;
                    }
                
                // sample the chain
                 if ( (n % mcmcPhases.getSampleFrequency() == 0  || n == sampleLength) && phase == "sample" )
                    samples.addSample(powIdx, curLnL[0]);

                }
            }
        }
    auto stop = std::chrono::high_resolution_clock::now();
        
    // print summary information
    std::cout << std::endl;
    std::cout << "   Markov chain Monte Carlo run summary" << std::endl;
    std::cout << "   * Run time = " << formattedTime(start, stop) << std::endl;
    updateInfo.print();
    std::cout << std::endl;

    double marginalLnL = samples.marginalLikelihood();
    std::cout << "   * Marginal likelihood = " << marginalLnL << std::endl;

    // clean up
    closeOutputFiles();
    delete [] parmValues;
    delete [] curLnL;
    delete [] curLnP;
}

void Mcmc::runPosterior(void) {

    std::cout << "   Markov chain Monte Carlo" << std::endl;
    
    // initialize the chain
    initialize();
    
    for (int chain=0; chain<numChains; chain++)
        modelPtr[chain]->setUpdateLikelihood();
    
    double* curLnL = new double[numChains];
    double* curLnP = new double[numChains];
    for (int chain=0; chain<numChains; chain++)
        {
        curLnL[chain] = modelPtr[chain]->lnLikelihood();
        curLnP[chain] = modelPtr[chain]->lnPriorProbability();
        }
        
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();

    // Metropolis-Hastings algorithm
    auto start = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=numMcmcCycles; n++)
        {
        
        for (int chain=0; chain<numChains; chain++)
            {
            // propose a new value for the chain
            double lnProposalRatio = modelPtr[chain]->update(n);
            
            // calculate the likelihood and prior ratios (natural log scale)
            double newLnL = modelPtr[chain]->lnLikelihood();
            double lnLikelihoodRatio = newLnL - curLnL[chain];
            double newLnP = modelPtr[chain]->lnPriorProbability();
            double lnPriorRatio = newLnP - curLnP[chain];
            
            double myPower = power(chain);
            lnLikelihoodRatio *= myPower;
            lnPriorRatio *= myPower;

            // accept or reject the state
            bool accept = false;
            if ( log(rv->uniformRv()) < lnLikelihoodRatio + lnPriorRatio + lnProposalRatio )
                accept = true;
            
            // print to the screen to give the user some clue where the chain is
            if (chain + 1 == numChains && (n == 1 || n == numMcmcCycles || n % printFrequency == 0))
                {
                auto timePt = std::chrono::high_resolution_clock::now();
                print(n, curLnL, curLnP, accept, timePt, start);
                }
            
            // update the state of the chain
            if (accept == false)
                {
                modelPtr[chain]->reject();
                updateInfo.reject(modelPtr[chain]->getLastUpdate());
                }
            else
                {
                modelPtr[chain]->accept();
                updateInfo.accept(modelPtr[chain]->getLastUpdate());
                curLnL[chain] = newLnL;
                curLnP[chain] = newLnP;
                }
            
            // sample the chain, printing the current state to files
            if (modelPtr[chain] == getColdModel() && (n == 1 || n == numMcmcCycles || n % sampleFrequency == 0))
                sample(n, curLnL[chain], curLnP[chain]);
            }
            
        // attempt to swap
        if (numChains > 1)
            {
            // choose two chains to sawp
            int idx1, idx2;
            chooseModelsToSwap(idx1, idx2);
            double lnR  = power(idx2) * (curLnL[idx1] + curLnP[idx1]) + power(idx1) * (curLnL[idx2] + curLnP[idx2]); // swapped state
                   lnR -= power(idx1) * (curLnL[idx1] + curLnP[idx1]) + power(idx2) * (curLnL[idx2] + curLnP[idx2]); // original condition
            int chainIndex1 = modelPtr[idx1]->getIndex();
            int chainIndex2 = modelPtr[idx2]->getIndex();
            if ( log(rv->uniformRv()) < lnR )
                {
                modelPtr[idx1]->setIndex(chainIndex2);
                modelPtr[idx2]->setIndex(chainIndex1);
                }
            }
        
        }
    auto stop = std::chrono::high_resolution_clock::now();
        
    // print summary information
    std::cout << std::endl;
    std::cout << "   Markov chain Monte Carlo run summary" << std::endl;
    std::cout << "   * Run time = " << formattedTime(start, stop) << std::endl;
    updateInfo.print();
    std::cout << std::endl;
    
    // clean up
    closeOutputFiles();
    delete [] parmValues;
    delete [] curLnL;
    delete [] curLnP;
}

double Mcmc::safeExponentiation(double lnX) {

    if (lnX > 0.0)
        return 1.0;
    else if (lnX < -300.0)
        return 0.0;
    else
        return exp(lnX);
}

void Mcmc::sample(int gen, double lnL, double lnP) {

    UserSettings& settings = UserSettings::userSettings();
    std::string outPath = settings.getOutFile();

    Tree* t = getColdModel()->getTree();
    double tl = t->getTreeLength();
    std::string ts = t->getNewick(6);
    
    if (gen == 1)
        {
        // output header information to the .tsv file containing parameter values
        parmStrm << "Gen\t" << "lnL\t" << "lnP\t" << "TL\t";
        parmStrm << getColdModel()->getParameterHeader();
        parmStrm << std::endl;

        // output leading information for the NEXUS-formatted tree file, with translation block
        std::vector<std::string>& tn = t->getTaxonNames();
        treeStrm << "begin trees;" << std::endl;
        treeStrm << "   translate" << std::endl;
        for (int i=0; i<tn.size(); i++)
            {
            treeStrm << "      " << i+1 << " " << tn[i];
            if (i+1 != tn.size())
                treeStrm << ",";
            else
                treeStrm << ";";
            treeStrm << std::endl;
            }

        // output leading information for the JSON-formatted file containing the alignment
        std::vector<ParameterAlignment*> alns = getColdModel()->getAlignments();
        for (int i=0; i<alns.size(); i++)
            {
            std::string fn = outPath + "." + alns[i]->parmName + ".aln";
            std::ofstream alnStream( fn.c_str(), std::ios::out);
            if (!alnStream)
                Msg::error("Problem opening alignment file " + fn + " for initial output");

            std::string stateSetsStr = getColdModel()->getStateSetsJsonString();
            if (stateSetsStr != "")
                {
                alnStream << "{" << stateSetsStr;
                alnStream << ", \"Samples\": [\n";
                }
            else
                {
                alnStream << "{\"Samples\": [\n";
                }
            alnStream << alns[i]->getJsonString() << "," << std::endl;
            
            alnStream.close();
            }
        }
        
    // output to parameter file
    std::cout << std::fixed << std::setprecision(8);
    getColdModel()->fillParameterValues(parmValues, numParmValues);
    parmStrm << gen << '\t' << lnL << '\t'<< lnP << '\t'<< tl << '\t';
    for (int i=0; i<numParmValues; i++)
        parmStrm << parmValues[i] << '\t';
    parmStrm << std::endl;
    
    getColdModel()->getParameterString();
    
    // output to tree file
    treeStrm << "   tree t_" << gen << " = " << ts << ";";
    treeStrm << std::endl;
    if (gen == numMcmcCycles)
        treeStrm << "end;" << std::endl;
    
    // output to alignment files
    std::vector<ParameterAlignment*> alns = getColdModel()->getAlignments();
    for (int i=0; i<alns.size(); i++)
        {
        std::string fn = outPath + "." + alns[i]->parmName + ".aln";
        std::ofstream alnStream( fn.c_str(), std::ios::app);
        if (!alnStream)
            Msg::error("Problem opening alignment file " + fn + " for appended output");

        alns[i]->jsonStream( alnStream );
        if (gen == numMcmcCycles)
            alnStream << "]\n}" << std::endl;
        else
            alnStream << "," << std::endl;
        alnStream.close();
        }
}
