#include <iomanip>
#include <iostream>
#include "Mcmc.hpp"
#include "McmcOutput.hpp"
#include "McmcTimer.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "Update.hpp"
#include "UpdateManager.hpp"
#include "UserSettings.hpp"



Mcmc::Mcmc(RandomVariable* r, Model* m) : rng(r), myModel(m) {

    UserSettings& settings = UserSettings::userSettings();
    numCycles          = settings.getNumMcmcCycles();
    numDigits          = (numCycles == 0) ? 1 : static_cast<int>(log10(numCycles)) + 1;
    printFrequency     = settings.getPrintFrequency();
    sampleFrequency    = settings.getSampleFrequency();
    baseOutputFileName = settings.getOutFile();
    updateMngr         = new UpdateManager(myModel, rng);
}

Mcmc::~Mcmc(void) {

    delete updateMngr;
}

