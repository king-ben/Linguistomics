#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include "json.hpp"
#include "Msg.hpp"
#include "UserSettings.hpp"

/* Summary of commands used by this program

   Cmd Arg    Json Key             Command Description
   ------------------------------------------------------------------------------------------------------------------------
   -d         N/A                  Name of the file with the initial alignments of words (and model settings)
   -o         FileOutput           File name for output
   -m         Model                Substitution model (jc69, f81, gtr or custom, making certain the model name is lower case)
   -l         BrlenPriorVal        Parameter of exponential branch length prior
   -lc        Clock                Clock constrained (yes or no, all lower case)
   -c         OnlyCompleteWords    Use only completely sampled words (yes or no, lower case)
   -g         SubRatesNumCats      Number of gamma rate categories for substitution rate variation (1 is no rate variation)
   -i         IndelRatesNumCats    Number of gamma indel categories (1 is no indel rate variation)
   -z         CalcMarginal         Calculate marginal likelihood (yes or no)
   -e         Seed                 Seed for pseudorandom number generator (program uses time if this is not specified)
   -nc        NumChains            Number of Markov chains for (MC)^3
   -h         Temperature          Temperature for (MC)^3
   -n         NumCycles            Number of MCMC cycles (used for posterior analysis)
   -p         PrintFreq            Print-to-screen frequency
   -s         SampleFreq           Chain sample frequency
   -fb        FirstBurnLength      Length of initial burnin for path sampling method
   -pbl       PreburnLength        Preburnin length
   -tl        TuneLength           Tune length
   -bl        BurnLength           Length of burn period for each stone
   -sl        SampleLength         The number of MCMC cycles during which log likelihoods are sampled

   Json commands to be added by Shawn: SubRatesNumCats, IndelRatesNumCats, FirstBurn, PreburnLength, TuneLength, BurnLength, SampleLength */



UserSettings::UserSettings(void) {

    // dafault values
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

    // move the command-line arguments to a good-ol' C++ STL vector of strings
    std::vector<std::string> commands;
    executablePath = argv[0];
    for (int i=0; i<argc; i++)
        commands.push_back(argv[i]);
    
    // check that we have an odd number of commands
    if (commands.size() % 2 == 0 || commands.size() == 1)
        {
        usage();
        Msg::error("Problem with arguments for \"" + executablePath + "\"");
        }
    
    // read the command-line arguments, starting with the second one
    std::string arg = "";
    for (int i=1; i<commands.size(); i++)
        {
        std::string origCmd = commands[i];
        std::string cmd;
        cmd.resize(origCmd.size());

        // Convert the source string to lower case
        // storing the result in destination string
        std::transform(origCmd.begin(), origCmd.end(), cmd.begin(), [](char c) { return (char)std::tolower(c); });
        
        if (arg == "")
            {
            arg = cmd;
            }
        else
            {
            if (arg == "-d")
                dataFile = cmd;
            else if (arg == "-o")
                outFile = cmd;
            else if (arg == "-m")
                {
                if (cmd == "jc69")
                    substitutionModel = jc69;
                else if (cmd == "f81")
                    substitutionModel = f81;
                else if (cmd == "gtr")
                    substitutionModel = gtr;
                else if (cmd == "custom")
                    substitutionModel = custom;
                else
                    Msg::error("Unknown substitution model " + cmd);
                }
            else if (arg == "-e")
                seed = atoi(cmd.c_str());
            else if (arg == "-nc")
                numChains = atoi(cmd.c_str());
            else if (arg == "-h")
                temperature = atof(cmd.c_str());
            else if (arg == "-n")
                numMcmcCycles = atoi(cmd.c_str());
            else if (arg == "-p")
                printFrequency = atoi(cmd.c_str());
            else if (arg == "-s")
                sampleFrequency = atoi(cmd.c_str());
            else if (arg == "-z")
                checkPointFrequency = atoi(cmd.c_str());
            else if (arg == "-l")
                branchLengthLambda = atof(cmd.c_str());
            else if (arg == "-fb")
                firstBurnLength = atoi(cmd.c_str());
            else if (arg == "-pbl")
                preburninLength = atoi(cmd.c_str());
            else if (arg == "-tl")
                tuneLength = atoi(cmd.c_str());
            else if (arg == "-bl")
                burnLength = atoi(cmd.c_str());
            else if (arg == "-sl")
                sampleLength = atoi(cmd.c_str());
            else if (arg == "-z")
                {
                if (cmd == "yes")
                    calculateMarginalLikelihood = true;
                else if (cmd == "no")
                    calculateMarginalLikelihood = false;
                else
                    Msg::error("Unknown option for calculating marginal likelihood");
                }

            else if (arg == "-lc")
                {
                if (cmd == "yes")
                    useClockConstraint = true;
                else if (cmd == "no")
                    useClockConstraint = false;
                else
                    Msg::error("Unknown option for branch length constraints");
                }

            else if (arg == "-c")
                {
                if (cmd == "yes")
                    useOnlyCompleteWords = true;
                else if (cmd == "no")
                    useOnlyCompleteWords = false;
                else
                    Msg::error("Unknown option for use of complete words");
                }
            else if (arg == "-g")
                numRateCategories = atoi(cmd.c_str());
            else if (arg == "-i")
                numIndelCategories = atoi(cmd.c_str());
            else
                {
                usage();
                Msg::error("Unknown command \"" + arg + "\"");
                }
            
            arg = "";
            }
        }
        
    // check if the data file contains the commands. If so, we will use those, overwriting what we
    // just parsed (excepting the path/name to the configuration file)
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
    auto it = j.find("McmcSettings");
    if (it != j.end())
        {
        nlohmann::json jsonSettings = j["McmcSettings"];

        auto it2 = jsonSettings.find("FileOutput");
        if (it2 != jsonSettings.end())
            {
            //outFile = jsonSettings["FileOutput"];
            Msg::warning("Setting output file path via the configuration file is no longer supported");
            }

        it2 = jsonSettings.find("OnlyCompleteWords");
        if (it2 != jsonSettings.end())
            {
            std::string res = getVariable(jsonSettings, "OnlyCompleteWords");
            if (res == "no")
                useOnlyCompleteWords = false;
            else if (res == "yes")
                useOnlyCompleteWords = true;
            else
                Msg::error("Unknown option" + res);
            }

        it2 = jsonSettings.find("Model");
        if (it2 != jsonSettings.end())
            {
            std::string res = getVariable(jsonSettings, "Model");
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
            
        it2 = jsonSettings.find("Clock");
        if (it2 != jsonSettings.end())
            {
            std::string res = getVariable(jsonSettings, "Clock");
            if (res == "no")
                useClockConstraint = false;
            else if (res == "yes")
                useClockConstraint = true;
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
            checkPointFrequency = jsonSettings["CheckPtFreq"];

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
            numChains = jsonSettings["NumChains"];

        it2 = jsonSettings.find("Temperature");
        if (it2 != jsonSettings.end())
            temperature = jsonSettings["Temperature"];
        }
}

void UserSettings::print(void) {

    std::cout << "   Settings" << std::endl;
    std::cout << "   * Executable path                         = " << executablePath << std::endl;
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
    if (useClockConstraint == false)
        std::cout << "   * Branch length constraints               = unconstrained" << std::endl;
    else
        std::cout << "   * Branch length constraints               = constant rate (clock)" << std::endl;
    std::cout << "   * Number of gamma rate categories         = " << numRateCategories << std::endl;
    std::cout << "   * Number of gamma indel categories        = " << numIndelCategories << std::endl;
    std::cout << "   * Branch length prior parameter           = " << branchLengthLambda << std::endl;
    if (seed == 0)
        std::cout << "   * Random number seed                      = Current time" << std::endl;
    else
        std::cout << "   * Random number seed                      = " << seed << std::endl;
    if (calculateMarginalLikelihood == false)
        {
        std::cout << "   * Calculate marginal likelihood           = no" << std::endl;
        std::cout << "   * Number of chains (MC)^3                 = " << numChains << std::endl;
        std::cout << "   * Temperature                             = " << temperature << std::endl;
        std::cout << "   * Number of MCMC cycles                   = " << numMcmcCycles << std::endl;
        std::cout << "   * Print-to-screen frequency               = " << printFrequency << std::endl;
        std::cout << "   * Chain sample frequency                  = " << sampleFrequency << std::endl;
        std::cout << "   * Check point frequency                   = " << checkPointFrequency << std::endl;
        }
   else
        {
        std::cout << "   * Calculate marginal likelihood           = yes" << std::endl;
        std::cout << "   * Pre-burnin length                       = " << preburninLength << std::endl;
        std::cout << "   * Tuning length                           = " << tuneLength << std::endl;
        }
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Usage" << std::endl;
    std::cout << "   * -d                      -- File with initial alignments of words" << std::endl;
    std::cout << "   * -o / FileOutput         -- File name for output" << std::endl;
    std::cout << "   * -m / Model              -- Substitution model (jc69/f81/gtr/custom)" << std::endl;
    std::cout << "   * -l / BrlenPriorVal      -- Parameter of exponential branch length prior" << std::endl;
    std::cout << "   * -lc / Clock             -- Clock constrained (yes) or unconstrained (no) branch lengths" << std::endl;
    std::cout << "   * -c / OnlyCompleteWords  -- Use only completely sampled words (no/yes)" << std::endl;
    std::cout << "   * -g                      -- Number of gamma rate categories (=1 is no rate variation)" << std::endl;
    std::cout << "   * -i                      -- Number of gamma indel categories (=1 is no indel rate variation)" << std::endl;
    std::cout << "   * -z / CalcMarginal       -- Calculate marginal likelihood (no/yes)" << std::endl;
    std::cout << "   * -e / Seed               -- Seed for pseudo random number generator" << std::endl;
    std::cout << "   * -nc / NumChains         -- Number of Markov chains for (MC)^3" << std::endl;
    std::cout << "   * -h / Temperature        -- Temperature for (MC)^3" << std::endl;
    std::cout << "   * -n / NumCycles          -- Number of MCMC cycles" << std::endl;
    std::cout << "   * -p / PrintFreq          -- Print-to-screen frequency" << std::endl;
    std::cout << "   * -s / SampleFreq         -- Chain sample frequency" << std::endl;
    std::cout << "   * -z / CheckPtFreq        -- Check point file frequency" << std::endl;
    std::cout << "   * -nt / NumTunings        -- Number of tunings" << std::endl;
    std::cout << "   * -tl / TuneLength        -- Tune length" << std::endl;
    std::cout << "   * -bl / PreburnLength     -- Preburnin length" << std::endl;
    std::cout << std::endl;
}

