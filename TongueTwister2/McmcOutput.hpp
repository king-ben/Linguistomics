#ifndef McmcOutput_hpp
#define McmcOutput_hpp

#include <cstdio>
#include <list>
#include <vector>
class Model;
class Node;
class ParameterAlignment;
class ParameterTree;
class RateMatrix;
class Tree;
typedef std::vector<std::list<int>::iterator> AlignmentLruIterator;



class McmcOutput {

    public:
                                McmcOutput(void) = delete;
                                McmcOutput(Model* m, const char* baseName);
                               ~McmcOutput(void);
        void                    open(void);
        void                    close(void);
        void                    sample(int generation, double lnL, double lnP);
        
    private:
        void                    closeStandardFiles(void);
        void                    openStandardFiles(void);
        
                                // .config
        void                    writeConfigurationFile(void);
        
                                // scalar variables
        void                    writeScalarHeader(void);
        void                    writeScalarSample(int gen, double lnL, double lnP);
        
                                // tree
        void                    writeNewick(FILE* f, Tree* t);
        void                    writeNewickNode(FILE* f, Node* p, Node* root);
        void                    writeNewickWithNames(FILE* f, Tree* t);
        void                    writeNewickNodeWithNames(FILE* f, Node* p, Node* root);

        void                    writeTreeSample(int gen);

                                // alignments
        void                    closeAlignmentFile(int alnIdx);
        FILE*                   ensureAlignmentFileOpen(int alnIdx);
        void                    evictAlignmentFileIfNeeded(void);
        void                    touchAlignmentFileLRU(int alnIdx);
        void                    writeAlignmentSample(int alnIdx);
        
        bool                    reopenJsonArrayForAppend(FILE* f);
        
        Model*                  model;
        char                    baseName[256];
        bool                    filesOpen;
        
        FILE*                   scalarFile;
        FILE*                   treeFile;
        std::vector<FILE*>      familyTreeFiles;

        FILE**                  alignmentFiles;
        std::list<int>          alignmentLRU;      // most-recently-used at front
        AlignmentLruIterator    alignmentLRUIter;
        size_t                  alignmentCacheCapacity;
        bool*                   alignmentFirstSample;
        int                     numAlignments;
        
        ParameterTree*          treeParm;
        ParameterAlignment**    alignmentParms;
        RateMatrix*             rateMatrix;

        bool                    sampleAlignments;
};

#endif
