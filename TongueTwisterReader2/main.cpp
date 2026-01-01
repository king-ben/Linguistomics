#include <string>
#include <iostream>
#include "Analysis.hpp"
#include "AnalysisComparison.hpp"
#include "Container.hpp"
#include "RandomVariable.hpp"
#include "Threads.hpp"
#include "UserSettings.hpp"

void compareAnalyses(std::vector<Analysis*>& analyses);
void printHeader(void);



int main(int argc, char* argv[]) {

    // instantiate objects needed for reading files
    RandomVariable rng;
    ThreadPool threadPool;
    
    // print a salutation
    printHeader();
    
    // read the user settings
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
    
    // read all of the directories with output from TongueTwister analyses
    std::vector<std::string> directories = settings.getInputDirectoryName();
    std::vector<Analysis*> analyses;
    for (size_t i=0; i<directories.size(); i++)
        {
        Analysis* analysis = new Analysis(&rng, &threadPool, directories[i], settings.getBurnFraction());
        analyses.push_back(analysis);

        analysis->print();
        analysis->printSorted();
        analysis->writeNytril(settings.getNytrilOutputFileName());
        analysis->writeR(settings.getROutFile(), analysis->getName(), analysis->modelName(), (int)i);
        }
        
    // compare analyses (if there are more than one)
    AnalysisComparison::compareAnalyses(settings.getROutFile(), analyses);
        
    // clean up
    for (size_t i=0; i<analyses.size(); i++)
        delete analyses[i];
    
    return EXIT_SUCCESS;
}

void printHeader(void) {

    std::cout << "   TongueTwister Reader 2.0" << std::endl;
    std::cout << "   * Running on " << std::thread::hardware_concurrency() << " threads" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << std::endl;
}

