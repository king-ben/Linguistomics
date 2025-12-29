#ifndef ParameterAlignment_hpp
#define ParameterAlignment_hpp

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "AlnMatrix.hpp"
#include "Parameter.hpp"
#include "RbBitSet.h"
class Alignment;
class AlnMatrix;
class AlignmentProposal;
class IndelMatrix;
class Model;



class ParameterAlignment : public Parameter {

    public:
                                        ParameterAlignment(void) = delete;
                                        ParameterAlignment(const ParameterAlignment& pa) = delete;
                                        ParameterAlignment(RandomVariable* r, Model* m, Alignment* a, std::string n, int idx);
                                       ~ParameterAlignment(void);
        void                            accept(void);
        void                            fillParameterValues(double* x, int& start, int maxNumValues);
        AlnMatrix*                      getAlignment(void) { return alignment[0]; }
        AlnMatrix*                      getAlignment(int idx) { return alignment[idx]; }
        char                            getCharFromCode(int code);
        int                             getGapCode(void) { return gapCode; }
        std::string                     getHeader(void) { return ""; }
        std::vector<std::vector<int> >  getIndelMatrix(void);
        std::vector<std::vector<int> >  getIndelMatrix(int idx);
        void                            getIndelMatrix(IndelMatrix* indelMat);
        std::vector<std::vector<int> >  getIndelMatrix(std::vector<std::vector<int> >& aln);
        void                            getIndelMatrix(std::vector<std::vector<int> >& aln, std::vector<std::vector<int> >& m);
        int                             getIndex(void) { return index; }
        bool                            getIsCompletelySampled(void) { return completelySampled; }
        std::string                     getJsonString(void);
        int                             getNumSites(void) { return alignment[0]->getNumSites(); }
        int                             getNumStates(void) { return numStates; }
        int                             getNumTaxa(void) { return numTaxa; }
        int                             getNumValues(void) { return 1; }
        int                             getPrintWidth(void) { return printWidth; }
        std::vector<std::vector<int> >  getRawSequenceMatrix(void) { return sequences; }
        std::string                     getString(void) { return ""; }
        char*                           getCString(void) { return parmStr; }
        RbBitSet                        getTaxonMask(void);
        std::string                     getTaxonMaskString(void);
        std::vector<std::string>        getTaxonNames(void) { return taxonNames; }
        std::string                     getUpdateName(int idx);
        void                            jsonStream(std::ofstream& strm);
        int                             lengthOfLongestSequence(void);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(int iter);
    
    protected:
        AlignmentProposal*              alignmentProposal;
        AlnMatrix*                      alignment[2];  // numTaxa X numSites
        std::vector<std::vector<int> >  sequences;     // numTaxa X numSites for taxon i (i.e., left-justified), instantiated once on construction
        bool                            completelySampled;
        double                          tuning;
        double                          exponent;
        double                          gapPenalty;
        int                             numStates;
        int                             gapCode;       // the gap code is simply the maximum number of states possible in the Markov model
        int                             numTaxa;
        int                             printWidth;
        int                             index;
        std::vector<std::string>        taxonNames;
        RbBitSet                        taxonMask;
        std::map<int,int>               taxonMapKeyCanonical;
        std::map<int,int>               taxonMapKeyAlignment;
};

#endif
