#include <iostream>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"
#include "Threads.hpp"

void printHeader(void);



int main(int argc, char* argv[]) {

    // print the header
    printHeader();

    // Read the user settings from the command-line arguments. See the
    // header of UserSettings.cpp for a full list of command options.
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();

    // create the thread pool
    ThreadPool threadPool;

    // instantiate the random number generator
    RandomVariable* rv = NULL;
    if (settings.getSeed() == 0)
        rv = new RandomVariable;
    else
        rv = new RandomVariable(settings.getSeed());
        
    // set up the phylogenetic model
    auto chains = settings.getNumChains();
    Model** model = new Model*[chains];
    for (int i=0; i< chains; i++)
        {
        auto m = new Model(rv, threadPool);
        m->setIndex(i);
        model[i] = m;
        }
    
    // run the Markov chain Monte Carlo algorithm
    Mcmc chain(model, rv);
    chain.run();
    
    delete rv;
    for (int i=0; i< chains; i++)
        delete model[i];
    delete [] model;
    
    return 0;
}

void printHeader(void) {

    std::cout << "   TongueTwister 1.0" << std::endl;
    std::cout << "   * Running on " << std::thread::hardware_concurrency() << " threads" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << std::endl;
}
