#include <iomanip>
#include <iostream>
#include "McmcMarginal.hpp"
#include "McmcPhase.hpp"
#include "McmcTimer.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "SteppingStones.hpp"
#include "Stone.hpp"
#include "Update.hpp"
#include "UpdateManager.hpp"
#include "UserSettings.hpp"



McmcMarginal::McmcMarginal(RandomVariable* r, Model* m) : Mcmc(r, m) {

    UserSettings& settings = UserSettings::userSettings();
    mcmcPhases = new McmcPhase(settings.getFirstBurnLength(), settings.getTuneLength(), settings.getBurnLength(), settings.getSampleLength(), settings.getSampleFrequency());
    samples = new SteppingStones(127, 0.3, 1.0);
    
    // for printing nicely to the screen
    numCycles = calculateNumCycles(); // call this before calculatePrintDigits()
    std::pair<int,int> nDigits = calculatePrintDigits();
    numDigits = nDigits.first;
    numDigitsPhase = nDigits.second;
}

McmcMarginal::~McmcMarginal(void) {

    delete mcmcPhases;
    delete samples;
}

int McmcMarginal::calculateNumCycles(void) {

    // numCycles isn't used to control the chain length as it is in the
    // McmcPosterior class. Here, numCycles is used to simply inform the
    // user how many MCMC cycles have passed. Note, however, that we have
    // to calculate numCycles from the MCMC phase information and the number
    // of stones in the path sampling algorithm.
    int nc = 0;
    std::vector<std::string>& phases = mcmcPhases->getPhases();
    for (std::string phase : phases)
        nc += mcmcPhases->getLengthForPhase(1, phase);
    nc *= samples->size();
    nc += mcmcPhases->getFirstBurnLength();
    return nc;
}

std::pair<int,int> McmcMarginal::calculatePrintDigits(void) {

    std::pair<int,int> numDigits;

    size_t maxPhaseLength = 0;
    std::vector<std::string>& phases = mcmcPhases->getPhases();
    for (std::string phase : phases)
        {
        size_t len = mcmcPhases->getLengthForPhase(1, phase);
        if (phase == "burn")
            len += mcmcPhases->getFirstBurnLength();
        if (len > maxPhaseLength)
            maxPhaseLength = len;
        }
    
    numDigits.first = (numCycles == 0) ? 1 : static_cast<int>(log10(numCycles)) + 1;
    numDigits.second = (maxPhaseLength == 0) ? 1 : static_cast<int>(log10(maxPhaseLength)) + 1;;
    return numDigits;
}

void McmcMarginal::printStatus(size_t stoneIdx, size_t cycle, size_t cnt, std::string& phase, double lnL1, double lnL2, McmcTimer* timer) {

    std::cout << "   ";
    std::cout << std::setw(3) << stoneIdx+1 << "/" << samples->size() << " ";
    std::cout << std::setw(numDigits) << cnt << " ";
    std::cout << std::setw(numDigitsPhase) << cycle << " ";
    std::cout << std::fixed << std::setprecision(7) << (*samples)[stoneIdx]->getPower() << " ";
    std::cout << std::setw(6) << phase << " ";
    std::cout << "-- " << std::fixed << std::setprecision(2) << lnL1 << " -> " << lnL2 << " -- ";
    timer->printStatus(cnt, numCycles);
    std::cout << std::endl; 
}

void McmcMarginal::run(void) {

    std::cout << "   Path Sampling Markov chain Monte Carlo (MCMC) Analysis" << std::endl;
    std::cout << "   * Initial Burn Length: " << mcmcPhases->getFirstBurnLength() << std::endl;
    std::cout << "   * Burn Length: " << mcmcPhases->getBurnLength() << std::endl;
    std::cout << "   * Tune Length: " << mcmcPhases->getTuneLength() << std::endl;
    std::cout << "   * Sample Length: " << mcmcPhases->getSampleLength() << std::endl;
    std::cout << "   * Base Output File: " << baseOutputFileName << std::endl;

    double currentLnL = myModel->lnLikelihood();
    double currentLnP = myModel->lnPrior();
    
    McmcTimer timer;
    timer.start();

    std::vector<std::string>& phases = mcmcPhases->getPhases();
    for (size_t powIdx=0, cnt=1; powIdx<samples->size(); powIdx++)
        {
        Stone* stone = (*samples)[powIdx];
        double power = stone->getPower();
        for (std::string phase : phases)
            {
            for (size_t n=1; n<=mcmcPhases->getLengthForPhase(powIdx, phase); n++, cnt++)
                {
                // propose new state
                Update* update = updateMngr->randomlyChooseUpdate();
                double lnProposalRatio = update->update(power);
                updateMngr->updateDependants(update);
                
                // calculate acceptance probability
                double newLnL = myModel->lnLikelihood();
                double newLnP = myModel->lnPrior();
                double lnLikelihoodRatio = (newLnL - currentLnL) * power;
                double lnPriorRatio = newLnP - currentLnP;
                double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

                // accept/reject proposed state
                bool accept = false;
                if (log(rng->uniformRv()) < lnR)
                    accept = true;
                
                // print to the screen
                if (n % printFrequency == 0)
                    printStatus(powIdx, n, cnt, phase, currentLnL, newLnL, &timer);

                // adjust states
                if (accept == true)
                    {
                    currentLnL = newLnL;
                    currentLnP = newLnP;
                    updateMngr->accept(update);
                    }
                else 
                    {
                    updateMngr->reject(update);
                    }
                
                // sample the chain
                if ( (n % mcmcPhases->getSampleFrequency() == 0  || n == mcmcPhases->getSampleLength()) && phase == "Sample" )
                    samples->addSample(powIdx, currentLnL);
                }
            }
        }

    timer.end();
    updateMngr->summary();
    double marginalLnL = samples->marginalLikelihood();
    std::cout << "   * Marginal likelihood = " << marginalLnL << std::endl;
}

