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



class McmcOutput {

    public:
                                McmcOutput(void) = delete;
                                McmcOutput(Model* m, const char* baseName);
                               ~McmcOutput(void);
        void                    open(void);
        void                    close(void);
        void                    sample(int generation);
        
    private:
        void                    openStandardFiles(void);
        void                    closeStandardFiles(void);
        
        void                    writeConfigurationFile(void);
        
        void                    writeScalarHeader(void);
        void                    writeScalarSample(int gen);
        void                    writeTreeSample(int gen);
        void                    writeAlignmentSample(int alnIdx);
        
        void                    writeNewick(FILE* f, Tree* t);
        void                    writeNewickNode(FILE* f, Node* p, Node* root);

                                // alignment file descriptor management
        FILE*                   ensureAlignmentFileOpen(int alnIdx);
        void                    touchAlignmentFileLRU(int alnIdx);
        void                    evictAlignmentFileIfNeeded(void);
        void                    closeAlignmentFile(int alnIdx);
        bool                    reopenJsonArrayForAppend(FILE* f);
        
        Model*                  model;
        char                    baseName[256];
        bool                    filesOpen;
        
        FILE*                   scalarFile;
        FILE*                   treeFile;

        FILE**                  alignmentFiles;
        std::list<int>          alignmentLRU;      // most-recently-used at front
        std::vector<std::list<int>::iterator> alignmentLRUIter;
        size_t                  alignmentCacheCapacity;
        bool*                   alignmentFirstSample;
        int                     numAlignments;
        
        ParameterTree*          treeParm;
        ParameterAlignment**    alignmentParms;
        RateMatrix*             rateMatrix;
};

#endif
