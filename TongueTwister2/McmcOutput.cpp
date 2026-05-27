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
#include "UserSettings.hpp"




McmcOutput::McmcOutput(Model* m, const char* base) : model(m) {

    strncpy(baseName, base, 255);
    baseName[255] = '\0';
    filesOpen = false;
    
    scalarFile = nullptr;
    treeFile = nullptr;
    familyTreeFiles.clear();
    alignmentFiles = nullptr;
    alignmentFirstSample = nullptr;
    alignmentCacheCapacity = 64;
    
    // read alignment sampling setting
    UserSettings& settings = UserSettings::userSettings();
    sampleAlignments = settings.getSampleAlignments();
    
    // write the configuration file
    writeConfigurationFile();
    
    // find the tree parameter for single-family analyses.
    // In multi-family mode, each SubModel has its own tree.
    if (model->getNumFamilies() == 0)
        {
        treeParm = model->findParameter<ParameterTree>();
        if (treeParm == nullptr)
            Msg::error("Could not find tree parameter for McmcOutput");
        }
    else
        {
        treeParm = nullptr;
        }
    
    // count alignment parameters and store pointers
    const std::vector<Parameter*>& parms = model->getParameters();
    numAlignments = 0;
    for (size_t i=0; i<parms.size(); i++)
        {
        if (dynamic_cast<ParameterAlignment*>(parms[i]) != nullptr)
            numAlignments++;
        }
    
    // only set up alignment output infrastructure if alignment sampling is on
    if (sampleAlignments)
        {
        alignmentParms = new ParameterAlignment*[numAlignments];
        int idx = 0;
        for (size_t i=0; i<parms.size(); i++)
            {
            ParameterAlignment* pa = dynamic_cast<ParameterAlignment*>(parms[i]);
            if (pa != nullptr)
                alignmentParms[idx++] = pa;
            }
        
        alignmentFiles = new FILE*[numAlignments];
        alignmentFirstSample = new bool[numAlignments];
        for (int i=0; i<numAlignments; i++)
            {
            alignmentFiles[i] = nullptr;
            alignmentFirstSample[i] = true;
            }
        alignmentLRUIter.resize(numAlignments, alignmentLRU.end());
        }
    else
        {
        alignmentParms = nullptr;
        numAlignments = 0;  // prevents any loop from iterating
        }
    
    // get a pointer to the rate matrix
    rateMatrix = model->getRateMatrix();
}

McmcOutput::~McmcOutput(void) {

    if (filesOpen)
        close();
    if (sampleAlignments)
        {
        delete [] alignmentParms;
        delete [] alignmentFiles;
        delete [] alignmentFirstSample;
        }
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

    for (FILE*& f : familyTreeFiles)
        {
        if (f != nullptr)
            {
            fprintf(f, "end;\n");
            fclose(f);
            f = nullptr;
            }
        }

    familyTreeFiles.clear();
    
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
    

    // open tree file(s)
    if (model->getNumFamilies() == 0)
        {
        // single-family case: original behaviour
        snprintf(fileName, 512, "%s.tre", baseName);
        treeFile = fopen(fileName, "w");
        if (treeFile == nullptr)
            Msg::error("Could not open tree output file");

        fprintf(treeFile, "begin trees;\n");
        fprintf(treeFile, "   translate\n");

        const std::vector<std::string>& taxonNames =
            treeParm->getCanonicalTaxonList();

        for (size_t i=0; i<taxonNames.size(); i++)
            {
            fprintf(treeFile, "      %2zu %s", i+1, taxonNames[i].c_str());

            if (i == taxonNames.size() - 1)
                fprintf(treeFile, ";\n");
            else
                fprintf(treeFile, ",\n");
            }
        }
    else
        {
        // multi-family case: one .tre file per family
        int numFamilies = model->getNumFamilies();
        familyTreeFiles.resize(numFamilies, nullptr);

        for (int f=0; f<numFamilies; f++)
            {
            snprintf(fileName, 512, "%s.family%d.tre", baseName, f+1);

            familyTreeFiles[f] = fopen(fileName, "w");
            if (familyTreeFiles[f] == nullptr)
                Msg::error("Could not open family tree output file");

            fprintf(familyTreeFiles[f], "begin trees;\n");
            fprintf(familyTreeFiles[f], "   translate\n");

            ParameterTree* familyTree =
                model->getSubModel(f)->getTree();

            const std::vector<std::string>& taxonNames =
                familyTree->getCanonicalTaxonList();

            for (size_t i=0; i<taxonNames.size(); i++)
                {
                fprintf(familyTreeFiles[f], "      %2zu %s",
                        i+1,
                        taxonNames[i].c_str());

                if (i == taxonNames.size() - 1)
                    fprintf(familyTreeFiles[f], ";\n");
                else
                    fprintf(familyTreeFiles[f], ",\n");
                }
            }
        }
    
    // only initialise alignment file tracking if alignment sampling is on
    if (sampleAlignments)
        {
        alignmentLRU.clear();
        std::fill(alignmentLRUIter.begin(), alignmentLRUIter.end(), alignmentLRU.end());
        for (int i=0; i<numAlignments; i++)
            {
            alignmentFiles[i] = nullptr;
            alignmentFirstSample[i] = true;
            }
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

    if (model->getNumFamilies() == 0)
        {
        fflush(treeFile);
        }
    else
        {
        for (FILE* f : familyTreeFiles)
            {
            if (f != nullptr)
                fflush(f);
            }
        }

    if (sampleAlignments)
        {
        for (int i=0; i<numAlignments; i++)
            writeAlignmentSample(i);
        }
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

void McmcOutput::writeNewickWithNames(FILE* f, Tree* t) {

    Node* root = t->getRoot();
    writeNewickNodeWithNames(f, root, root);
    fprintf(f, ";");
}

void McmcOutput::writeNewickNodeWithNames(FILE* f, Node* p, Node* root) {

    if (p->getIsLeaf())
        {
        fprintf(f, "%s:%.6f", p->getName(), p->getBranchLength());
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

            writeNewickNodeWithNames(f, d, root);
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
    fprintf(scalarFile, "lnP");

    const std::vector<Parameter*>& parms = model->getParameters();

    // Count how many of each parameter type we have
    int numIndel = 0;
    int numFreqs = 0;
    int numExch  = 0;

    for (Parameter* p : parms)
        {
        if (dynamic_cast<ParameterIndelRates*>(p) != nullptr)
            numIndel++;
        else if (dynamic_cast<ParameterFrequencies*>(p) != nullptr)
            numFreqs++;
        else if (dynamic_cast<ParameterExchangeabilities*>(p) != nullptr)
            numExch++;
        }

    // Indel-rate headers
    int indelIdx = 0;
    for (Parameter* p : parms)
        {
        ParameterIndelRates* indelRatesParm =
            dynamic_cast<ParameterIndelRates*>(p);

        if (indelRatesParm != nullptr)
            {
            indelIdx++;

            if (numIndel == 1)
                {
                fprintf(scalarFile, "\tLambda");
                fprintf(scalarFile, "\tMu");
                }
            else
                {
                fprintf(scalarFile, "\tLambda[%d]", indelIdx);
                fprintf(scalarFile, "\tMu[%d]", indelIdx);
                }
            }
        }

    // Frequency headers
    int freqIdx = 0;
    for (Parameter* p : parms)
        {
        ParameterFrequencies* freqsParm =
            dynamic_cast<ParameterFrequencies*>(p);

        if (freqsParm != nullptr)
            {
            freqIdx++;

            const std::vector<double>& freqs = freqsParm->getFrequencies();

            for (size_t i=0; i<freqs.size(); i++)
                {
                if (numFreqs == 1)
                    fprintf(scalarFile, "\tF[%zu]", i);
                else
                    fprintf(scalarFile, "\tF%d[%zu]", freqIdx, i);
                }
            }
        }

    // Exchangeability-rate headers
    int exchIdx = 0;
    for (Parameter* p : parms)
        {
        ParameterExchangeabilities* exchParm =
            dynamic_cast<ParameterExchangeabilities*>(p);

        if (exchParm != nullptr)
            {
            exchIdx++;

            size_t numRates = exchParm->getNumRates();

            RateMatrixHelper* helper = nullptr;
            if (rateMatrix != nullptr)
                helper = rateMatrix->getHelper();

            if (helper != nullptr)
                {
                // custom/natural class model
                std::vector<std::string> columnHeader = helper->getLabels();

                if (numRates != columnHeader.size())
                    Msg::error("Problem printing column header for the Natural Class model");

                for (size_t i=0; i<columnHeader.size(); i++)
                    {
                    if (numExch == 1)
                        fprintf(scalarFile, "\tNCR[%s]", columnHeader[i].c_str());
                    else
                        fprintf(scalarFile, "\tNCR%d[%s]", exchIdx, columnHeader[i].c_str());
                    }
                }
            else
                {
                // GTR-style model, or fallback
                size_t numStates = exchParm->getNumStates();
                size_t expectedGtrRates = numStates * (numStates - 1) / 2;

                if (numRates == expectedGtrRates)
                    {
                    for (size_t i=0; i<numStates; i++)
                        {
                        for (size_t j=i+1; j<numStates; j++)
                            {
                            if (numExch == 1)
                                fprintf(scalarFile, "\tR[%zu-%zu]", i, j);
                            else
                                fprintf(scalarFile, "\tR%d[%zu-%zu]", exchIdx, i, j);
                            }
                        }
                    }
                else
                    {
                    // generic fallback
                    for (size_t i=0; i<numRates; i++)
                        {
                        if (numExch == 1)
                            fprintf(scalarFile, "\tR[%zu]", i);
                        else
                            fprintf(scalarFile, "\tR%d[%zu]", exchIdx, i);
                        }
                    }
                }
            }
        }

    fprintf(scalarFile, "\n");
}

void McmcOutput::writeScalarSample(int gen, double lnL, double lnP) {

    fprintf(scalarFile, "%d\t", gen);
    fprintf(scalarFile, "%.2lf\t", lnL);
    fprintf(scalarFile, "%.2lf", lnP);

    const std::vector<Parameter*>& parms = model->getParameters();

    // Print all indel-rate parameters
    for (Parameter* p : parms)
        {
        ParameterIndelRates* indelRatesParm =
            dynamic_cast<ParameterIndelRates*>(p);

        if (indelRatesParm != nullptr)
            {
            fprintf(scalarFile, "\t%.6f", indelRatesParm->getInsertionRate());
            fprintf(scalarFile, "\t%.6f", indelRatesParm->getDeletionRate());
            }
        }

    // Print all stationary-frequency parameters
    for (Parameter* p : parms)
        {
        ParameterFrequencies* freqsParm =
            dynamic_cast<ParameterFrequencies*>(p);

        if (freqsParm != nullptr)
            {
            const std::vector<double>& freqs = freqsParm->getFrequencies();

            for (size_t i=0; i<freqs.size(); i++)
                fprintf(scalarFile, "\t%.6f", freqs[i]);
            }
        }

    // Print all exchangeability-rate parameters
    for (Parameter* p : parms)
        {
        ParameterExchangeabilities* exchParm =
            dynamic_cast<ParameterExchangeabilities*>(p);

        if (exchParm != nullptr)
            {
            const std::vector<double>& rates = exchParm->getRates();

            for (size_t i=0; i<rates.size(); i++)
                fprintf(scalarFile, "\t%.6f", rates[i]);
            }
        }

    fprintf(scalarFile, "\n");
}

void McmcOutput::writeTreeSample(int gen) {

    if (model->getNumFamilies() == 0)
        {
        // single-family case
        fprintf(treeFile, "   tree gen_%d = ", gen);
        writeNewick(treeFile, treeParm->getTree());
        fprintf(treeFile, "\n");
        }
    else
        {
        // multi-family case: write one sample to each family tree log
        int numFamilies = model->getNumFamilies();

        if ((int)familyTreeFiles.size() != numFamilies)
            Msg::error("Mismatch between number of families and tree output files");

        for (int f=0; f<numFamilies; f++)
            {
            FILE* out = familyTreeFiles[f];

            if (out == nullptr)
                Msg::error("Family tree output file is not open");

            ParameterTree* familyTree =
                model->getSubModel(f)->getTree();

            fprintf(out, "   tree gen_%d = ", gen);
            writeNewick(out, familyTree->getTree());
            fprintf(out, "\n");
            }
        }
}



