#include <iostream>
#include "FileManager.hpp"
#include "JsonData.hpp"
#include "McmcMarginal.hpp"
#include "McmcPosterior.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "Threads.hpp"
#include "UserSettings.hpp"

void initialize(int argc, char* argv[]);
void mcmc(RandomVariable* rng, Model* model);
void printFarewell(RandomVariable* rng);
void printHeader(void);



int main(int argc, char* argv[]) {

    // print a salutation
    printHeader();

    // instantiate the tread pool
    ThreadPool threadPool;
    
    // initialize the user settings and the JSON data object
    initialize(argc, argv);
    
    // instantiate the random number generator (32-bit Mersenne Twister)
    RandomVariable rng(UserSettings::userSettings().getSeed());

    // set up the phylogenetic model
    Model model(&rng, UserSettings::userSettings().getSubstitutionModel(), &threadPool);
    
    // run the MCMC analysis
    mcmc(&rng, &model);
    
    // print a farewell message
    printFarewell(&rng);
    
    return EXIT_SUCCESS;
}

void initialize(int argc, char* argv[]) {
    
    // read the user settings
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
    
    // check the existence of important files/paths
    if (FileManager::fileExists(settings.getDataFile()) == false)
        Msg::error("Data file, \"" + settings.getDataFile() + ",\" does not exist");
    if (FileManager::directoryExists(FileManager::getParentPath(settings.getOutFile())) == false)
        Msg::error("Output directory, \"" + FileManager::getParentPath(settings.getOutFile()) + ",\" does not exist");

    // read the JSON file
    JsonData& jsonData = JsonData::jsonInstance();
    jsonData.readJsonFile(settings.getDataFile());
}

void mcmc(RandomVariable* rng, Model* model) {

    if (UserSettings::userSettings().getCalculateMarginalLikelihood() == false)
        {
        // approximate the posterior probability distribution
        McmcPosterior mcmc(rng, model);
        mcmc.run();
        }
    else 
        {
        // approximate the marginal likelihood using path sampling
        McmcMarginal mcmc(rng, model);
        mcmc.run();
        }
}

void printFarewell(RandomVariable* rng) {

    std::vector<std::string> farewells = {
        "Convergence achieved. Prior beliefs appropriately modified.",
        "All alignments were marginalized in the making of this posterior.",
        "May the Bayes factor be ever in your favor.",
        "May your posterior be informative and your burnin be brief.",
        "Rooting for you always. Goodbye!",
        "Go Bears!",
        "Convergence achieved. Shutting down."
        // feel free to add more here
    };
    
    std::cout << "   " << farewells[static_cast<int>(rng->uniformRv() * farewells.size())] << std::endl;
    std::cout << std::endl;
}

void printHeader(void) {

    std::cout << "   TongueTwister 2.0" << std::endl;
    std::cout << "   * Running on " << std::thread::hardware_concurrency() << " threads" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << std::endl;
}
