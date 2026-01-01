#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(void) {

    // dafault values
    outFile       = "";
    nytrilOutFile = "";
    rOutFile      = "";
    burnFraction  = 0.0;
    fullOutput    = true;
}

void UserSettings::readCommandLineArguments(int argc, char* argv[]) {

    // move the command-line arguments to a good-ol' C++ STL vector of strings
    std::vector<std::string> commands;
    for (int i=0; i<argc; i++)
        commands.push_back(argv[i]);
    
    // check that we have an odd number of commands
    if (commands.size() % 2 == 0 || commands.size() == 1)
        {
        usage();
        Msg::error("Problem with arguments");
        }
    
    // read the command-line arguments, starting with the second one
    std::string arg = "";
    for (int i=1; i<commands.size(); i++)
        {
        std::string cmd = commands[i];
        if (arg == "")
            {
            arg = cmd;
            }
        else
            {
            if (arg == "-d")
                inFile.push_back(cmd);
            else if (arg == "-o")
                outFile = cmd;
            else if (arg == "-n")
                nytrilOutFile = cmd;
            else if (arg == "-p")
                rOutFile = cmd;
            else if (arg == "-b")
                burnFraction = atof(cmd.c_str());
            else if (arg == "-f")
                {
                std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return std::tolower(c);});
                if (cmd == "true")
                    fullOutput = true;
                else if (cmd == "false")
                    fullOutput = false;
                else
                    Msg::error("Unknown option for full output command");
                }

            else
                {
                usage();
                Msg::error("Unknown command \"" + arg + "\"");
                }
            
            arg = "";
            }
        }
}

void UserSettings::print(void) {

    std::cout << "   Settings" << std::endl;
    if (fullOutput == true)
        std::cout << "   * Full output to file     = True" << std::endl;
    else
        std::cout << "   * Full output to file     = False" << std::endl;
    for (size_t i=0; i<inFile.size(); i++)
        std::cout << "   * Input path name         = \"" << inFile[i] << "\"" << std::endl;
    std::cout << "   * Output file name        = \"" << outFile << "\"" << std::endl;
    std::cout << "   * Nytril Output file name = \"" << nytrilOutFile << "\"" << std::endl;
    std::cout << "   * R Output file name      = \"" << rOutFile << "\"" << std::endl;
    std::cout << "   * Burn in fraction        = " << burnFraction << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Usage" << std::endl;
    std::cout << "   * -i  -- File with data to be summarized" << std::endl;
    std::cout << "   * -o  -- File name for output" << std::endl;
    std::cout << "   * -b  -- Burn-in number" << std::endl;
    std::cout << std::endl;
}

