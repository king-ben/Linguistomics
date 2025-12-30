#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include "Alignment.hpp"
#include "FileManager.hpp"
#include "json.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "ParameterStatistics.hpp"
#include "Partition.hpp"
#include "Statistics.hpp"
#include "Threads.hpp"
#include "Tree.hpp"



McmcSummary::McmcSummary(ThreadPool* tp, std::string dn) : pool(tp), fileManager(dn), statePartitions(nullptr) {

    readConfigurationFile();
    readParameterFile();
    readAlignmentFiles();
    readTreeFile();
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
        
    for (Tree* t : trees)
        delete t;
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

std::vector<std::string> McmcSummary::breakString(std::string str) {

    std::vector<std::string> broken;
    std::string word = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == ',' || c == ';')
            {
            if (word != "")
                broken.push_back(word);
            broken.push_back( std::string(1,c) );
            word = "";
            }
        else
            {
            word +=  c;
            }
        }
    if (word != "")
        broken.push_back(word);
    
    return broken;
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

bool McmcSummary::hasSemicolon(std::string str) {

    std::size_t found = str.find(';');
    if (found != std::string::npos)
        return true;
    return false;
}

bool McmcSummary::hasTrees(void) {

    if (trees.size() > 0)
        return true;
    return false;
}

std::map<int,std::string> McmcSummary::interpretTranslateString(std::vector<std::string> translateTokens) {

    std::map<int,std::string> translateMap;
    bool readingKey = true;
    std::string key = "";
    std::string val = "";
    for (int i=0; i<translateTokens.size(); i++)
        {
        std::string token = translateTokens[i];
        if (token == "," || token == ";")
            {
            readingKey = true;
            int intKey = atoi(key.c_str());
            translateMap.insert( std::make_pair(intKey,val) );
            }
        else
            {
            if (readingKey == true)
                {
                key = token;
                readingKey = false;
                }
            else
                val = token;
            }
        }
        
    return translateMap;
}

std::string McmcSummary::interpretTreeString(std::string str) {

    std::string ns = "";
    bool reading = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '(')
            reading = true;
            
        if (reading == true)
            {
            ns += c;
            }
            
        if (c == ';')
            reading = false;
        }
    
    return ns;
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
        {
        Msg::warning("Did not find any alignment files");
        return;
        }

    int progressWidth = 50;
    std::cout << "   Reading Alignment Files" << std::endl;
    std::cout << "   [";
    for (size_t i=0; i<progressWidth; i++)
        {
        if (i % 10 == 0 && i != 0)
            std::cout << "|";
        else
            std::cout << "-";
        }
    std::cout << "]" << std::endl;
    std::cout << "   [";
    std::cout.flush();

    // shared atomic counter for progress tracking
    std::atomic<int> completedCount(0);
    int totalFiles = static_cast<int>(alignmentFiles.size());

    // set up the reading tasks
    alignments.resize(alignmentFiles.size());
    alignmentNames.resize(alignmentFiles.size());
    std::vector<ReadAlignmentTask*> readTasks(alignmentFiles.size());
    
    for (size_t i=0; i<alignmentFiles.size(); i++)
        {
        ReadAlignmentTask* task = new ReadAlignmentTask(
            alignmentFiles[i], i, &alignments, &alignmentNames, &completedCount);
        readTasks[i] = task;
        }
        
    // push all tasks to the thread pool
    for (ReadAlignmentTask* task : readTasks)
        pool->pushTask(task);
    
    // progress bar in main thread
    int lastProgress = 0;
    
    while (completedCount < totalFiles)
        {
        int completed = completedCount.load();
        int progress = (completed * progressWidth) / totalFiles;
        
        // print new progress characters
        for (int i = lastProgress; i < progress; i++)
            std::cout << "*";
        std::cout.flush();
        
        lastProgress = progress;
        
        // small sleep to avoid busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    
    // ensure we reach 100%
    for (int i=lastProgress; i<progressWidth; i++)
        std::cout << "=";
    std::cout << "] " << totalFiles << " files" << std::endl;
    
    // wait for all tasks to fully complete
    pool->wait();
    
    // clean up tasks
    for (ReadAlignmentTask* task : readTasks)
        delete task;
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

void McmcSummary::readTreeFile(void) {

    std::vector<std::string> treeFiles = fileManager.filesWithExtension("tre");
    if (treeFiles.size() != 1)
        Msg::warning("Expecting only one file with the extension \"tre\" but found " + std::to_string(treeFiles.size()) + " such files");
    std::string fn = treeFiles[0];
    
	// open the file
	std::ifstream fstrm(fn.c_str());
	if (!fstrm)
        Msg::error("Cannot open file \"" + fn + "\"");

	std::string lineString = "";
	int line = 0;
    bool readingTranslateTable = false, readingTree = false;
    std::vector<std::string> translateTokens;
    std::string treeString = "";
    int treeCount = 0;
	while( getline(fstrm, lineString).good() )
		{
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
            if (word == "translate")
                {
                readingTranslateTable = true;
                }
            else if (word == "tree")
                {
                readingTree = true;
                }
            else
                {
                if (readingTranslateTable == true)
                    {
                    std::vector<std::string> brokenWord = breakString(word);
                    if (hasSemicolon(word) == true)
                        {
                        for (int i=0; i<brokenWord.size(); i++)
                            translateTokens.push_back(brokenWord[i]);
                        readingTranslateTable = false;
                        translateMap = interpretTranslateString(translateTokens);
                        }
                    else
                        {
                        for (int i=0; i<brokenWord.size(); i++)
                            translateTokens.push_back(brokenWord[i]);
                        }
                    }
                else if (readingTree == true)
                    {
                    if (hasSemicolon(word) == true)
                        {
                        treeString += word;
                        readingTree = false;
                        treeCount++;

                        std::string newickStr = interpretTreeString(treeString);
                        Tree* t = new Tree(newickStr, translateMap);
                        trees.push_back(t);

                        treeString = "";
                        }
                    else
                        treeString += word;
                    }
                }

			wordNum++;
            } while ( (ch=linestream.get()) != EOF );
						
		line++;
		}

	// close the file
	fstrm.close();
 
}
