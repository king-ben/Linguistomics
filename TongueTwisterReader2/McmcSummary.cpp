#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "Alignment.hpp"
#include "FileManager.hpp"
#include "json.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "ParameterStatistics.hpp"
#include "Partition.hpp"
#include "Statistics.hpp"
#include "Threads.hpp"



McmcSummary::McmcSummary(ThreadPool* tp, std::string dn) : pool(tp), fileManager(dn), statePartitions(nullptr) {

    readConfigurationFile();
    readParameterFile();
    readAlignmentFiles();
}

McmcSummary::~McmcSummary(void) {

    for (size_t i=0; i<stats.size(); i++)
        delete stats[i];
    
    for (size_t i=0; i<alignments.size(); i++)
        {
        for (size_t j=0; j<alignments[i].size(); j++)
            delete alignments[i][j];
        }

    if (statePartitions != nullptr)
        delete statePartitions;
}

ParameterStatistics* McmcSummary::operator[](size_t idx) {

    if (idx >= stats.size())
        return nullptr;
    return stats[idx];
}

ParameterStatistics* McmcSummary::operator[](size_t idx) const {

    if (idx >= stats.size())
        return nullptr;
    return stats[idx];
}

bool McmcSummary::hasFrequencies(void) {

    for (size_t i=0; i<size(); i++)
        {
        ParameterStatistics* stat = stats[i];
        std::string parmName = stat->getName();
        if (parmName[0] == 'F')
            return true;
        }
    return false;
}

bool McmcSummary::hasExchangeabilities(void) {

    for (size_t i=0; i<size(); i++)
        {
        ParameterStatistics* stat = stats[i];
        std::string parmName = stat->getName();
        if (parmName[0] == 'R' || parmName[0] == 'N')
            return true;
        }
    return false;
}

bool McmcSummary::hasPartition(void) {

    if (statePartitions != nullptr)
        return true;
    return false;
}

void McmcSummary::printParameterSummary(void) {

    printParameterSummary(0.0);
}

void McmcSummary::printParameterSummary(double burnFraction) {

    // find longest parameter name
    size_t len = 0;
    for (ParameterStatistics* parm : stats)
        {
        if (parm->getName().length() > len)
            len = parm->getName().length();
        }
        
    for (ParameterStatistics* parm : stats)
        {
        std::pair<double,double> mv = Statistics::getMeanAndVariance(parm->getValues(), burnFraction);
        CredibleInterval ci = Statistics::getCredibleInterval(parm->getValues(), burnFraction);

        std::cout << parm->getName();
        for (size_t i=0; i<len-parm->getName().length(); i++)
            std::cout << " ";
        std::cout << " -- ";
        std::cout << std::scientific << std::setprecision(3);
        std::cout << mv.first << " " << mv.second << " (";
        std::cout << ci.lower << "," << ci.upper << ")" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        }
}

void McmcSummary::readAlnFile(std::string fn, size_t n, size_t total) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file " + fn + " at byte " + std::to_string(ex.byte));
        return;
        }

    if (!j.is_array())
        Msg::error("Expected JSON array of sampled alignments in " + fn);

    std::cout << "   * Alignment " << n+1 << "/" << total << ": \"" << fn << "\"" << std::endl;
    std::vector<Alignment*> fileAlignments;
    bool nameRecorded = false;
    
    for (const auto& sample : j)
        {
        std::string name = sample["Name"];
        if (nameRecorded == false)
            {
            alignmentNames[n] = name;
            nameRecorded = true;
            }
        Alignment* aln = new Alignment(sample["Data"]);
        fileAlignments.push_back(aln);
        }
    alignments[n] = fileAlignments;
}

void McmcSummary::readAlignmentFiles(void) {

    std::vector<std::string> alignmentFiles = fileManager.filesWithExtension("json");
    if (alignmentFiles.size() == 0)
        Msg::warning("Did not find any alignment files");

    std::cout << "   Reading Alignment Files" << std::endl;

#   if 0
    alignments.resize(alignmentFiles.size());
    alignmentNames.resize(alignmentFiles.size());
    for (size_t i=0; i<alignmentFiles.size(); i++)
        readAlnFile(alignmentFiles[i], i, alignmentFiles.size());
    std::cout << std::endl;
#   else  

    // set up the reading tasks
    alignments.resize(alignmentFiles.size());
    alignmentNames.resize(alignmentFiles.size());
    std::vector<ReadAlignmentTask*> readTasks(alignmentFiles.size());
    for (size_t i=0; i<alignmentFiles.size(); i++)
        {
        ReadAlignmentTask* task = new ReadAlignmentTask(alignmentFiles[i], i, &alignments, &alignmentNames);
        readTasks[i] = task;
        }
        
    for (ReadAlignmentTask* task : readTasks)
        pool->pushTask(task);
    pool->wait();
    
#   endif
}

void McmcSummary::readConfigurationFile(void) {

    std::vector<std::string> configurationFiles = fileManager.filesWithExtension("config");
    if (configurationFiles.size() != 1)
        Msg::warning("Expecting only one file with the extension \"config\" but found " + std::to_string(configurationFiles.size()) + " such files");
    std::string fn = configurationFiles[0];

    std::cout << "   Reading configuration file \"" << fn << "\"" << std::endl;

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return;
        }

    auto jtaxa = j["Taxa"];
    for (int i = 0; i < jtaxa.size(); i++)
        taxa.push_back(jtaxa[i].get<std::string>());
        
    std::cout << "   * Taxa: ";
    for (int i = 0; i < jtaxa.size(); i++)
        std::cout << taxa[i] << " ";
    std::cout << std::endl;
    
    auto it = j.find("NumberOfStates");
    if (it == j.end())
        Msg::error("Could not find the number of states in teh configuration file");
    else
        numStates = j["NumberOfStates"];
    std::cout << "   * Number of States: " << numStates << std::endl;
        
    it = j.find("PartitionSets");
    if (it == j.end())
        {
        statePartitions = nullptr;
        }
    else
        {
        nlohmann::json jsPart = j["PartitionSets"];
        statePartitions = new Partition(jsPart);
        }
    std::cout << std::endl;
}

void McmcSummary::readParameterFile(void) {

    std::vector<std::string> parameterFiles = fileManager.filesWithExtension("tsv");
    if (parameterFiles.size() != 1)
        Msg::warning("Expecting only one file with the extension \"tsv\" but found " + std::to_string(parameterFiles.size()) + " such files");
    std::string fn = parameterFiles[0];
    
	// open the file
	std::ifstream fstrm(fn.c_str());
	if (!fstrm)
        Msg::error("Cannot open file \"" + fn + "\"");
    std::cout << "   Reading parameter file \"" << fn << "\"" << std::endl;

	std::string lineString = "";
	int line = 0;
	while( getline(fstrm, lineString).good() )
		{
        //std::cout << line << " -- " << lineString << std::endl;
		std::istringstream linestream(lineString);
		int ch;
		std::string word = "";
		int wordNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
            if (word == "")
                continue;
			if (line == 0)
				{
                ParameterStatistics* s = new ParameterStatistics;
                s->setName(word);
                stats.push_back(s);
				}
			else
				{
                double x = atof(word.c_str());
                ParameterStatistics* s = stats[wordNum];
                s->addValue(x);
				}
			wordNum++;
            } while ( (ch=linestream.get()) != EOF );
						
		line++;
		}
	
	// close the file
	fstrm.close();
}
