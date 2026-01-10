#include <iomanip>
#include <iostream>
#include "McmcOutput.hpp"
#include "McmcPosterior.hpp"
#include "McmcTimer.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "Update.hpp"
#include "UpdateManager.hpp"



McmcPosterior::McmcPosterior(RandomVariable* r, Model* m) : Mcmc(r, m), tuneFrequency(1000) {

}

McmcPosterior::~McmcPosterior(void) {

}


void McmcPosterior::printStatus(int cycle, int nCycles, double lnL1, double lnL2, McmcTimer* timer) {

    std::cout << "   * " << std::setw(numDigits) << cycle << " -- ";
    std::cout << std::fixed << std::setprecision(2) << lnL1 << " -> " << lnL2 << " -- ";
    timer->printStatus(nCycles, numCycles + burnLength);
    std::cout << std::endl; 
}

void McmcPosterior::run(void) {

    std::cout << "   Markov chain Monte Carlo (MCMC) Analysis" << std::endl;
    std::cout << "   * Chain Length: " << numCycles << std::endl;
    std::cout << "   * Base Output File: " << baseOutputFileName << std::endl;

    McmcOutput output(myModel, baseOutputFileName.c_str());
    output.open();
    
    double currentLnL = myModel->lnLikelihood();
    double currentLnP = myModel->lnPrior();
    
    McmcTimer timer;
    timer.start();
    
    for (int n=-burnLength, len=1; n<=numCycles; n++, len++)
        {
        // propose new state
        Update* update = updateMngr->randomlyChooseUpdate();
        double lnProposalRatio = update->update();
        updateMngr->updateDependants(update);
        
        // calculate acceptance probability
        double newLnL = myModel->lnLikelihood();
        double newLnP = myModel->lnPrior();
        double lnLikelihoodRatio = newLnL - currentLnL;
        double lnPriorRatio = newLnP - currentLnP;
        double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;
        
        // accept/reject proposed state
        bool accept = false;
        if (log(rng->uniformRv()) < lnR)
            accept = true;
            
        // print to the screen
        if (n % printFrequency == 0)
            printStatus(n, len, currentLnL, newLnL, &timer);
            
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
            
        // sample chain
        if (n >= 0 && n % sampleFrequency == 0)
            output.sample(n, currentLnL, currentLnP);
            
        // tune updates
        if (n < 0 && n % tuneFrequency == 0)
            updateMngr->tune();
        if (n == 0)
            updateMngr->zeroOut();
        }

    timer.end();
    updateMngr->summary();
    output.close();
}

