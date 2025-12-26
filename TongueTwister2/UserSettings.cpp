#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include "argparse.hpp"
#include "json.hpp"
#include "Msg.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(void) {

    // default values
    numChains                   = 1;
    temperature                 = 0.1;
    dataFile                    = "";
    outFile                     = "";
    checkPointFrequency         = 10000;
    seed                        = 0;
    numMcmcCycles               = 1000;
    printFrequency              = 100;
    branchLengthLambda          = 50.0;
    substitutionModel           = jc69;
    matrixExpBackend            = autoBackend;
    calculateMarginalLikelihood = false;
    numRateCategories           = 1;
    numIndelCategories          = 1;
    useOnlyCompleteWords        = false;
    useClockConstraint          = false;
    firstBurnLength             = 0;
    preburninLength             = 0;
    tuneLength                  = 0;
    burnLength                  = 0;
    sampleLength                = 20000;
    sampleFrequency             = 100;
}

std::string UserSettings::getVariable(nlohmann::json &settings, const char* name) {

    std::string res = settings[name];
    std::transform(res.begin(), res.end(), res.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return res;
}

void UserSettings::readCommandLineArguments(int argc, char* argv[]) {

    executablePath = argv[0];
    
    // set up the argument parser
    argparse::ArgumentParser program("TongueTwister", "2.0");
    program.add_description("Bayesian phylogenetic inference with statistical alignment");
    
    // required argument: data file
    program.add_argument("-d", "--data")
        .required()
        .help("File with initial alignments of words (and model settings)");
    
    // output file
    program.add_argument("-o", "--output")
        .default_value(std::string("out"))
        .help("File name for output");
    
    // substitution model
    program.add_argument("-m", "--model")
        .default_value(std::string("jc69"))
        .choices("jc69", "f81", "gtr", "custom")
        .help("Substitution model (jc69/f81/gtr/custom)");
    
    // matrix exponential backend
    program.add_argument("-b", "--backend")
        .default_value(std::string("auto"))
        .choices("auto", "cpu-threaded", "cpu-batched", "gpu-batched")
        .help("Matrix exponential compute backend (auto/cpu-threaded/cpu-batched/gpu-batched)");
    
    // branch length prior parameter
    program.add_argument("-l", "--lambda")
        .default_value(50.0)
        .scan<'g', double>()
        .help("Parameter of exponential branch length prior");
    
    // clock constraint
    program.add_argument("-lc", "--clock")
        .default_value(std::string("no"))
        .choices("yes", "no")
        .help("Clock constrained (yes) or unconstrained (no) branch lengths");
    
    // complete words only
    program.add_argument("-c", "--complete")
        .default_value(std::string("no"))
        .choices("yes", "no")
        .help("Use only completely sampled words (yes/no)");
    
    // gamma rate categories for substitution
    program.add_argument("-g", "--gamma-sub")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of gamma rate categories for substitution (1 = no rate variation)");
    
    // gamma rate categories for indels
    program.add_argument("-i", "--gamma-indel")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of gamma indel categories (1 = no indel rate variation)");
    
    // calculate marginal likelihood
    program.add_argument("-z", "--marginal")
        .default_value(std::string("no"))
        .choices("yes", "no")
        .help("Calculate marginal likelihood (yes/no)");
    
    // random seed
    program.add_argument("-e", "--seed")
        .default_value(0u)
        .scan<'u', uint32_t>()
        .help("Seed for pseudo random number generator (0 = use current time)");
    
    // number of chains
    program.add_argument("-nc", "--num-chains")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of Markov chains for (MC)^3");
    
    // temperature
    program.add_argument("-ht", "--heat")
        .default_value(0.1)
        .scan<'g', double>()
        .help("Temperature for (MC)^3");
    
    // number of MCMC cycles
    program.add_argument("-n", "--num-cycles")
        .default_value(1000)
        .scan<'i', int>()
        .help("Number of MCMC cycles");
    
    // print frequency
    program.add_argument("-p", "--print-freq")
        .default_value(100)
        .scan<'i', int>()
        .help("Print-to-screen frequency");
    
    // sample frequency
    program.add_argument("-s", "--sample-freq")
        .default_value(100)
        .scan<'i', int>()
        .help("Chain sample frequency");
    
    // checkpoint frequency
    program.add_argument("-cp", "--checkpoint-freq")
        .default_value(10000)
        .scan<'i', int>()
        .help("Check point file frequency");
    
    // first burn length
    program.add_argument("-fb", "--first-burn")
        .default_value(0)
        .scan<'i', int>()
        .help("Length of initial burnin for path sampling method");
    
    // preburnin length
    program.add_argument("-pbl", "--preburn-length")
        .default_value(0)
        .scan<'i', int>()
        .help("Preburnin length");
    
    // tune length
    program.add_argument("-tl", "--tune-length")
        .default_value(0)
        .scan<'i', int>()
        .help("Tune length");
    
    // burn length
    program.add_argument("-bl", "--burn-length")
        .default_value(0)
        .scan<'i', int>()
        .help("Length of burn period for each stone");
    
    // sample length
    program.add_argument("-sl", "--sample-length")
        .default_value(20000)
        .scan<'i', int>()
        .help("Number of MCMC cycles during which log likelihoods are sampled");
    
    // parse arguments
    try
        {
        program.parse_args(argc, argv);
        }
    catch (const std::exception& err)
        {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
        }
    
    // extract values from parsed arguments
    dataFile            = program.get<std::string>("-d");
    outFile             = program.get<std::string>("-o");
    branchLengthLambda  = program.get<double>("-l");
    seed                = program.get<uint32_t>("-e");
    numChains           = program.get<int>("-nc");
    temperature         = program.get<double>("-ht");
    numMcmcCycles       = program.get<int>("-n");
    printFrequency      = program.get<int>("-p");
    sampleFrequency     = program.get<int>("-s");
    checkPointFrequency = program.get<int>("-cp");
    numRateCategories   = program.get<int>("-g");
    numIndelCategories  = program.get<int>("-i");
    firstBurnLength     = program.get<int>("-fb");
    preburninLength     = program.get<int>("-pbl");
    tuneLength          = program.get<int>("-tl");
    burnLength          = program.get<int>("-bl");
    sampleLength        = program.get<int>("-sl");
    
    // handle string-to-enum/bool conversions
    std::string modelStr = program.get<std::string>("-m");
    if (modelStr == "jc69")
        substitutionModel = jc69;
    else if (modelStr == "f81")
        substitutionModel = f81;
    else if (modelStr == "gtr")
        substitutionModel = gtr;
    else if (modelStr == "custom")
        substitutionModel = custom;
    
    std::string backendStr = program.get<std::string>("-b");
    if (backendStr == "auto")
        matrixExpBackend = autoBackend;
    else if (backendStr == "cpu-threaded")
        matrixExpBackend = cpuThreaded;
    else if (backendStr == "cpu-batched")
        matrixExpBackend = cpuBatched;
    else if (backendStr == "gpu-batched")
        matrixExpBackend = gpuBatched;
    
    std::string clockStr = program.get<std::string>("-lc");
    useClockConstraint = (clockStr == "yes");
    
    std::string completeStr = program.get<std::string>("-c");
    useOnlyCompleteWords = (completeStr == "yes");
    
    std::string marginalStr = program.get<std::string>("-z");
    calculateMarginalLikelihood = (marginalStr == "yes");
    
    // read settings from JSON file (these override command-line arguments)
    readJsonSettings();
}

void UserSettings::readJsonSettings(void) {

    std::string fn = getDataFile();
    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::error("Error parsing JSON file \"" + fn + "\" at byte " + std::to_string(ex.byte));
        }
    
    // read the model settings
    auto it = j.find("OnlyCompleteWords");
    if (it != j.end())
        {
        std::string res = getVariable(j, "OnlyCompleteWords");
        if (res == "no")
            useOnlyCompleteWords = false;
        else if (res == "yes")
            useOnlyCompleteWords = true;
        else
            Msg::error("Unknown option" + res);
        }

    it = j.find("Model");
    if (it != j.end())
        {
        std::string res = getVariable(j, "Model");
        if (res == "jc69")
            substitutionModel = jc69;
        else if (res == "f81")
            substitutionModel = f81;
        else if (res == "gtr")
            substitutionModel = gtr;
        else if (res == "custom")
            substitutionModel = custom;
        else
            Msg::error("Unknown substitution model " + res);
        }

    it = j.find("Clock");
    if (it != j.end())
        {
        std::string res = getVariable(j, "Clock");
        if (res == "no")
            useClockConstraint = false;
        else if (res == "yes")
            useClockConstraint = true;
        else
            Msg::error("Unknown option" + res);
        if (useClockConstraint == true)
            Msg::error("Molecular clock not yet supported");
        }

    it = j.find("BrlenPriorVal");
    if (it != j.end())
        branchLengthLambda = j["BrlenPriorVal"];
   
    // read the MCMC settings
    it = j.find("McmcSettings");
    if (it != j.end())
        {
        nlohmann::json jsonSettings = j["McmcSettings"];

        auto it2 = jsonSettings.find("FileOutput");
        if (it2 != jsonSettings.end())
            {
            Msg::warning("Setting output file path via the configuration file is no longer supported");
            }

        it2 = jsonSettings.find("CalcMarginal");
        if (it2 != jsonSettings.end())
            {
            std::string res = getVariable(jsonSettings, "CalcMarginal");
            if (res == "no")
                calculateMarginalLikelihood = false;
            else if (res == "yes")
                calculateMarginalLikelihood = true;
            else
                Msg::error("Unknown option" + res);
            }
            
        it2 = jsonSettings.find("Seed");
        if (it2 != jsonSettings.end())
            seed = jsonSettings["Seed"];

        it2 = jsonSettings.find("NumCycles");
        if (it2 != jsonSettings.end())
            numMcmcCycles = jsonSettings["NumCycles"];

        it2 = jsonSettings.find("PrintFreq");
        if (it2 != jsonSettings.end())
            printFrequency = jsonSettings["PrintFreq"];

        it2 = jsonSettings.find("SampleFreq");
        if (it2 != jsonSettings.end())
            sampleFrequency = jsonSettings["SampleFreq"];

        it2 = jsonSettings.find("CheckPtFreq");
        if (it2 != jsonSettings.end())
            {
            checkPointFrequency = jsonSettings["CheckPtFreq"];
            Msg::error("Check pointing not yet supported");
            }

        it2 = jsonSettings.find("BrlenPriorVal");
        if (it2 != jsonSettings.end())
            branchLengthLambda = jsonSettings["BrlenPriorVal"];
            
        it2 = jsonSettings.find("FirstBurnLength");
        if (it2 != jsonSettings.end())
            firstBurnLength = jsonSettings["FirstBurnLength"];

        it2 = jsonSettings.find("TuneLength");
        if (it2 != jsonSettings.end())
            tuneLength = jsonSettings["TuneLength"];

        it2 = jsonSettings.find("PreburnLength");
        if (it2 != jsonSettings.end())
            preburninLength = jsonSettings["PreburnLength"];
            
        it2 = jsonSettings.find("BurnLength");
        if (it2 != jsonSettings.end())
            burnLength = jsonSettings["BurnLength"];
            
        it2 = jsonSettings.find("SampleLength");
        if (it2 != jsonSettings.end())
            sampleLength = jsonSettings["SampleLength"];
            
        it2 = jsonSettings.find("NumChains");
        if (it2 != jsonSettings.end())
            {
            numChains = jsonSettings["NumChains"];
            Msg::error("Metropolis coupled MCMC not yet supported");
            }
        it2 = jsonSettings.find("Temperature");
        if (it2 != jsonSettings.end())
            {
            temperature = jsonSettings["Temperature"];
            Msg::error("Metropolis coupled MCMC not yet supported");
            }
        }
}

void UserSettings::print(void) {

    std::cout << "   Settings" << std::endl;
    //std::cout << "   * Executable path                         = " << executablePath << std::endl;
    std::cout << "   * File with initial word alignments       = \"" << dataFile << "\"" << std::endl;
    std::cout << "   * File name for output                    = \"" << outFile << "\"" << std::endl;
    if (useOnlyCompleteWords == true)
        std::cout << "   * Only analyzing completely sampled words = yes" << std::endl;
    else
        std::cout << "   * Only analyzing completely sampled words = no" << std::endl;
    if (substitutionModel == jc69)
        std::cout << "   * Substitution model                      = JC69" << std::endl;
    else if (substitutionModel == f81)
        std::cout << "   * Substitution model                      = F81" << std::endl;
    else if (substitutionModel == gtr)
        std::cout << "   * Substitution model                      = GTR" << std::endl;
    else
        std::cout << "   * Substitution model                      = Custom" << std::endl;
    if (matrixExpBackend == autoBackend)
        std::cout << "   * Matrix exponential backend              = auto" << std::endl;
    else if (matrixExpBackend == cpuThreaded)
        std::cout << "   * Matrix exponential backend              = cpu-threaded" << std::endl;
    else if (matrixExpBackend == cpuBatched)
        std::cout << "   * Matrix exponential backend              = cpu-batched" << std::endl;
    else if (matrixExpBackend == gpuBatched)
        std::cout << "   * Matrix exponential backend              = gpu-batched" << std::endl;
    if (useClockConstraint == false)
        std::cout << "   * Branch length constraints               = unconstrained" << std::endl;
    else
        std::cout << "   * Branch length constraints               = constant rate (clock)" << std::endl;
    //std::cout << "   * Number of gamma rate categories         = " << numRateCategories << std::endl;
    //std::cout << "   * Number of gamma indel categories        = " << numIndelCategories << std::endl;
    std::cout << "   * Branch length prior parameter           = " << branchLengthLambda << std::endl;
    if (seed == 0)
        std::cout << "   * Random number seed                      = Current time" << std::endl;
    else
        std::cout << "   * Random number seed                      = " << seed << std::endl;
    if (calculateMarginalLikelihood == false)
        {
        std::cout << "   * Calculate marginal likelihood           = no" << std::endl;
        //std::cout << "   * Number of chains (MC)^3                 = " << numChains << std::endl;
        //std::cout << "   * Temperature                             = " << temperature << std::endl;
        std::cout << "   * Number of MCMC cycles                   = " << numMcmcCycles << std::endl;
        std::cout << "   * Print-to-screen frequency               = " << printFrequency << std::endl;
        std::cout << "   * Chain sample frequency                  = " << sampleFrequency << std::endl;
        //std::cout << "   * Check point frequency                   = " << checkPointFrequency << std::endl;
        }
   else
        {
        std::cout << "   * Calculate marginal likelihood           = yes" << std::endl;
        std::cout << "   * Pre-burnin length                       = " << preburninLength << std::endl;
        std::cout << "   * Tuning length                           = " << tuneLength << std::endl;
        }
    std::cout << std::endl;
}
