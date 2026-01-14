#include <cstring>
#include <cctype>
#include <set>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include "Alignment.hpp"
#include "JsonData.hpp"
#include "McmcOutput.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterExchangeabilities.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "RateMatrix.hpp"
#include "RateMatrixHelper.hpp"
#include "Tree.hpp"



McmcOutput::McmcOutput(Model* m, const char* base) : model(m) {

    strncpy(baseName, base, 255);
    baseName[255] = '\0';
    filesOpen = false;
    
    scalarFile = nullptr;
    treeFile = nullptr;
    alignmentFiles = nullptr;
    alignmentFirstSample = nullptr;
    alignmentCacheCapacity = 64; // conservative default; stays well under typical ulimit -n
    
    // write the configuration file
    writeConfigurationFile();
    
    // find the tree parameter
    treeParm = model->findParameter<ParameterTree>();
    if (treeParm == nullptr)
        Msg::error("Could not find tree parameter for McmcOutput");
    
    // count alignment parameters and store pointers
    const std::vector<Parameter*>& parms = model->getParameters();
    numAlignments = 0;
    for (size_t i=0; i<parms.size(); i++)
        {
        if (dynamic_cast<ParameterAlignment*>(parms[i]) != nullptr)
            numAlignments++;
        }
    
    alignmentParms = new ParameterAlignment*[numAlignments];
    int idx = 0;
    for (size_t i=0; i<parms.size(); i++)
        {
        ParameterAlignment* pa = dynamic_cast<ParameterAlignment*>(parms[i]);
        if (pa != nullptr)
            alignmentParms[idx++] = pa;
        }
    
    // allocate file pointer arrays and first-sample trackers
    alignmentFiles = new FILE*[numAlignments];
    alignmentFirstSample = new bool[numAlignments];
    for (int i=0; i<numAlignments; i++)
        {
        alignmentFiles[i] = nullptr;
        alignmentFirstSample[i] = true;
        }

    alignmentLRUIter.resize(numAlignments, alignmentLRU.end());
        
    // get a pointer to the rate matrix
    rateMatrix = model->getRateMatrix();
}

McmcOutput::~McmcOutput(void) {

    if (filesOpen)
        close();
    delete [] alignmentParms;
    delete [] alignmentFiles;
    delete [] alignmentFirstSample;
}

void McmcOutput::close(void) {

    if (filesOpen == false)
        return;
    closeStandardFiles();
    filesOpen = false;
}

void McmcOutput::closeAlignmentFile(int alnIdx) {

    if (alnIdx < 0 || alnIdx >= numAlignments)
        return;

    // Remove from LRU list if present
    if (alignmentLRUIter[alnIdx] != alignmentLRU.end())
        {
        alignmentLRU.erase(alignmentLRUIter[alnIdx]);
        alignmentLRUIter[alnIdx] = alignmentLRU.end();
        }

    FILE* f = alignmentFiles[alnIdx];
    if (f != nullptr)
        {
        // finalize the JSON array for this file before closing
        fprintf(f, "\n]\n");
        fclose(f);
        alignmentFiles[alnIdx] = nullptr;
        }
}

void McmcOutput::closeStandardFiles(void) {

    if (scalarFile != nullptr)
        {
        fclose(scalarFile);
        scalarFile = nullptr;
        }
    
    if (treeFile != nullptr)
        {
        fprintf(treeFile, "end;\n");
        fclose(treeFile);
        treeFile = nullptr;
        }
    
    // Close any alignment files currently in the LRU cache.
    // (We finalize each JSON array with a closing bracket.)
    while (!alignmentLRU.empty())
        closeAlignmentFile(alignmentLRU.back());
}

FILE* McmcOutput::ensureAlignmentFileOpen(int alnIdx) {

    if (alnIdx < 0 || alnIdx >= numAlignments)
        Msg::error("Alignment index out of range in ensureAlignmentFileOpen");

    // if already open, just update LRU position and return.
    if (alignmentFiles[alnIdx] != nullptr)
        {
        touchAlignmentFileLRU(alnIdx);
        return alignmentFiles[alnIdx];
        }

    // Evict if cache is full.
    evictAlignmentFileIfNeeded();

    char fileName[512];
    snprintf(fileName, 512, "%s.aln%d.json", baseName, alnIdx);

    FILE* f = nullptr;

    if (alignmentFirstSample[alnIdx] == true)
        {
        // first time: overwrite and start JSON array
        f = fopen(fileName, "w");
        if (f == nullptr)
            Msg::error("Could not open alignment output file");
        fprintf(f, "[\n");
        }
    else
        {
        // Re-open an existing JSON array file and remove the closing bracket
        // so we can keep appending while preserving valid JSON on disk.
        f = fopen(fileName, "r+");
        if (f == nullptr)
            Msg::error("Could not re-open alignment output file for append");
        if (reopenJsonArrayForAppend(f) == false)
            Msg::error("Alignment output file did not look like a JSON array");
        }

    alignmentFiles[alnIdx] = f;
    touchAlignmentFileLRU(alnIdx);
    return f;
}

void McmcOutput::evictAlignmentFileIfNeeded(void) {

    while (alignmentLRU.size() >= alignmentCacheCapacity)
        {
        int victim = alignmentLRU.back();
        closeAlignmentFile(victim);
        }
}

void McmcOutput::open(void) {

    if (filesOpen)
        return;
    openStandardFiles();
    filesOpen = true;
}

void McmcOutput::openStandardFiles(void) {

    char fileName[512];
    
    // open scalar file and write header
    snprintf(fileName, 512, "%s.tsv", baseName);
    scalarFile = fopen(fileName, "w");
    if (scalarFile == nullptr)
        Msg::error("Could not open scalar output file");
    writeScalarHeader();
    
    // open tree file
    snprintf(fileName, 512, "%s.tre", baseName);
    treeFile = fopen(fileName, "w");
    if (treeFile == nullptr)
        Msg::error("Could not open tree output file");
    fprintf(treeFile, "begin trees;\n");
    fprintf(treeFile, "   translate\n");
    const std::vector<std::string>& taxonNames = treeParm->getCanonicalTaxonList();
    for (size_t i=0; i<taxonNames.size(); i++)
        {
        fprintf(treeFile, "      %2zu %s", i+1, taxonNames[i].c_str());
        if (i == taxonNames.size() - 1)
            fprintf(treeFile, ";\n");
        else 
            fprintf(treeFile, ",\n");
        }
    
    // Alignment files: do NOT open hundreds of files up front.
    // We lazily open them on demand using an LRU cache.
    alignmentLRU.clear();
    std::fill(alignmentLRUIter.begin(), alignmentLRUIter.end(), alignmentLRU.end());
    for (int i=0; i<numAlignments; i++)
        {
        alignmentFiles[i] = nullptr;
        alignmentFirstSample[i] = true;
        }
}

bool McmcOutput::reopenJsonArrayForAppend(FILE* f) {

    // We assume the file is a JSON array that ends with a closing ']' (possibly
    // followed by whitespace). We remove that closing bracket so future writes
    // can append additional entries, and later we re-add the bracket on close.

    if (f == nullptr)
        return false;

    if (fseek(f, 0, SEEK_END) != 0)
        return false;

    long endPos = ftell(f);
    if (endPos <= 0)
        return false;

    // Walk backwards over trailing whitespace
    long pos = endPos;
    int ch = 0;
    while (pos > 0)
        {
        if (fseek(f, pos - 1, SEEK_SET) != 0)
            return false;
        ch = fgetc(f);
        if (ch == EOF)
            return false;
        if (std::isspace(static_cast<unsigned char>(ch)))
            pos--;
        else
            break;
        }

    if (pos <= 0 || ch != ']')
        return false;

    int fd = fileno(f);
    if (fd < 0)
        return false;

    // truncate at the position of the closing bracket
    if (ftruncate(fd, pos - 1) != 0)
        return false;

    // seek to end so subsequent fprintf appends correctly
    if (fseek(f, 0, SEEK_END) != 0)
        return false;

    return true;
}

void McmcOutput::sample(int generation, double lnL, double lnP) {

    writeScalarSample(generation, lnL, lnP);
    fflush(scalarFile);

    writeTreeSample(generation);
    fflush(treeFile);

    for (int i=0; i<numAlignments; i++)
        writeAlignmentSample(i);
}

void McmcOutput::touchAlignmentFileLRU(int alnIdx) {

    // Remove if present
    if (alignmentLRUIter[alnIdx] != alignmentLRU.end())
        alignmentLRU.erase(alignmentLRUIter[alnIdx]);

    alignmentLRU.push_front(alnIdx);
    alignmentLRUIter[alnIdx] = alignmentLRU.begin();
}

void McmcOutput::writeAlignmentSample(int alnIdx) {

    FILE* f = ensureAlignmentFileOpen(alnIdx);
    ParameterAlignment* aln = alignmentParms[alnIdx];
    
    if (alignmentFirstSample[alnIdx] == false)
        fprintf(f, ",\n");
    alignmentFirstSample[alnIdx] = false;
    
    Alignment* a = aln->getAlignment(0);
    int numTaxa = static_cast<int>(a->getNumTaxa());
    int numSegments = static_cast<int>(a->getNumSegments());
    int gapCode = static_cast<int>(aln->getNumStates());
    const std::vector<std::string>& taxonNames = aln->getTaxonNames();
    size_t longestName = 0;
    for (size_t i=0; i<taxonNames.size(); i++)
        {
        if (taxonNames[i].length() > longestName)
            longestName = taxonNames[i].length();
        }
    
    fprintf(f, "{\"Name\": \"%s\", \"Data\": [\n", aln->getName().c_str());
    for (int i=0; i<numTaxa; i++)
        {
        fprintf(f, "{\"Taxon\": \"%s\", ", taxonNames[i].c_str());
        for (size_t j=0; j<longestName-taxonNames[i].length(); j++)
            fprintf(f, " ");
        fprintf(f, "\"Segments\": [");
        for (int j=0; j<numSegments; j++)
            {
            if (j > 0)
                fprintf(f, ", ");
            int c = (*a)(i, j);
            if (c == gapCode)
                fprintf(f, "%2d", -1);
            else
                fprintf(f, "%2d", c);
            }
        if (i == numTaxa-1)
            fprintf(f, "]}\n");
        else 
            fprintf(f, "]},\n");
        }
    
    fprintf(f, "]}");
}

void McmcOutput::writeConfigurationFile(void) {

    // write the configuration file
    char fileName[512];
    snprintf(fileName, 512, "%s.config", baseName);
    FILE* configFile = fopen(fileName, "w");
    if (configFile == nullptr)
        Msg::error("Could not open configuration file");
    
    JsonData& j = JsonData::jsonInstance();
    fprintf(configFile, "%s\n", j.dump().c_str());
    
    fclose(configFile);
}

void McmcOutput::writeNewick(FILE* f, Tree* t) {

    Node* root = t->getRoot();
    writeNewickNode(f, root, root);
    fprintf(f, ";");
}

void McmcOutput::writeNewickNode(FILE* f, Node* p, Node* root) {

    if (p->getIsLeaf())
        {
        //fprintf(f, "%s:%.6f", p->getName(), p->getBranchLength());
        fprintf(f, "%d:%.6f", p->getIndex()+1, p->getBranchLength());
        }
    else
        {
        fprintf(f, "(");
        bool first = true;
        for (size_t i=0; i<p->numDescendants(); i++)
            {
            Node* d = p->getDescendant(i);
            if (first == false)
                fprintf(f, ",");
            writeNewickNode(f, d, root);
            first = false;
            }
        fprintf(f, ")");
        if (p != root)
            fprintf(f, ":%.6f", p->getBranchLength());
        }
}

void McmcOutput::writeScalarHeader(void) {

    fprintf(scalarFile, "Gen\t");
    fprintf(scalarFile, "lnL\t");
    fprintf(scalarFile, "lnP\t");
    
    // print the insertion/deletion rates header
    ParameterIndelRates* indelRatesParm = model->findParameter<ParameterIndelRates>();
    if (indelRatesParm != nullptr)
        {
        fprintf(scalarFile, "\tLambda");
        fprintf(scalarFile, "\tMu");
        }
        
    // print the stationary frequencies header
    ParameterFrequencies* freqsParm = model->findParameter<ParameterFrequencies>();
    if (freqsParm != nullptr)
        {
        const std::vector<double>& freqs = freqsParm->getFrequencies();
        for (size_t i=0; i<freqs.size(); i++)
            fprintf(scalarFile, "\tF[%zu]", i);
        }

    // print the exchangeability rates header
    ParameterExchangeabilities* exchParm = model->findParameter<ParameterExchangeabilities>();
    if (exchParm != nullptr)
        {
        RateMatrixHelper* helper = rateMatrix->getHelper();
        if (helper == nullptr)
            {
            // format for GTR model (numStates * (numStates-1) / 2 = numRates columns)
            size_t numStates = exchParm->getNumStates();
            for (size_t i=0; i<numStates; i++)
                {
                for (size_t j=i+1; j<numStates; j++)
                    fprintf(scalarFile, "\tR[%zu-%zu]", i, j);
                }
            }
        else
            {
            // format for custom model (numClasses * (numClasses-1) / 2 = numRates columns)
            std::vector<std::string> columnHeader = helper->getLabels();
            size_t numRates = exchParm->getNumRates();
            if (numRates != columnHeader.size())
                Msg::error("Problem printing column header for the Natural Class model");
            for (size_t i=0; i<columnHeader.size(); i++)
                fprintf(scalarFile, "\tNCR[%s]", columnHeader[i].c_str());
            }
        }

    fprintf(scalarFile, "\n");
}

void McmcOutput::writeScalarSample(int gen, double lnL, double lnP) {

    fprintf(scalarFile, "%d\t", gen);
    fprintf(scalarFile, "%.2lf\t", lnL);
    fprintf(scalarFile, "%.2lf\t", lnP);
    
    // print the insertion/deletion rates
    ParameterIndelRates* indelRatesParm = model->findParameter<ParameterIndelRates>();
    if (indelRatesParm != nullptr)
        {
        fprintf(scalarFile, "\t%.6f", indelRatesParm->getInsertionRate());
        fprintf(scalarFile, "\t%.6f", indelRatesParm->getDeletionRate());
        }
        
    // print the stationary frequencies 
    ParameterFrequencies* freqsParm = model->findParameter<ParameterFrequencies>();
    if (freqsParm != nullptr)
        {
        const std::vector<double>& freqs = freqsParm->getFrequencies();
        for (size_t i=0; i<freqs.size(); i++)
            fprintf(scalarFile, "\t%.6f", freqs[i]);
        }

    // print the exchangeability rates 
    ParameterExchangeabilities* exchParm = model->findParameter<ParameterExchangeabilities>();
    if (exchParm != nullptr)
        {
        const std::vector<double>& rates = exchParm->getRates();
        for (size_t i=0; i<rates.size(); i++)
            fprintf(scalarFile, "\t%.6f", rates[i]);
        }
    
    fprintf(scalarFile, "\n");
}

void McmcOutput::writeTreeSample(int gen) {

    fprintf(treeFile, "   tree gen_%d = ", gen);
    writeNewick(treeFile, treeParm->getTree());
    fprintf(treeFile, "\n");
}


