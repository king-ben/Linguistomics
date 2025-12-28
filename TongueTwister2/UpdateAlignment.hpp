#ifndef UpdateAlignment_hpp
#define UpdateAlignment_hpp

#include "Update.hpp"
class Model;
class ParameterAlignment;
class ParameterFrequencies;
class ParameterTree;
class TransitionProbabilityManager;



class UpdateAlignment : public Update {

    public:
                                        UpdateAlignment(void) = delete;
                                        UpdateAlignment(Model* m, RandomVariable* r, ParameterAlignment* p, double exponent, double gapPenalty, double extProb);
                                       ~UpdateAlignment(void);
        std::string                     getUpdateName(void);
        std::string                     parameterType(void) { return "ParameterAlignment"; }
        void                            setDependants(void);
        double                          update(void);
        double                          update(double power);
        double                          updateFromPrior(void);
        
                                        // Public method to realign entire alignment (extensionProb = 1.0)
                                        // This is intended for use by UpdateTopology to jointly update
                                        // the tree topology and all alignments together.
                                        // Returns the log proposal ratio for the alignment change.
        double                          realignFull(void);
        
                                        // accessor for the alignment parameter (needed by UpdateTopology)
        ParameterAlignment*             getAlignmentParameter(void) { return myParm; }
    
    private:
        void                            buildScoringMatrix(Node* lftChild, Node* rhtChild);
        double                          calculateReverseProposal(Alignment* aln, int pos, int len);
        int                             countPaths(Alignment* aln, int startCol, int endCol);
        void                            getIndelMatrix(Alignment* aln, int startCol, int len);
        void                            initializeTreeStructure(void);
        double                          propose(void);
        double                          propose(double extProb, bool alignmentsMustBeDifferent);  // overload that accepts extension probability
        
        ParameterAlignment*             myParm;
        ParameterFrequencies*           freqsParm;
        ParameterTree*                  treeParm;
        TransitionProbabilityManager*   tiProbs;
        
        Tree*                           tree;
        unsigned                        taxonMask;
        
        int                             numTaxa;
        int                             numNodes;
        int                             numStates;
        int                             gapCode;
        double                          gapPenalty;
        double                          basis;
        double                          extensionProb;
        
        int                             maxLength;
        int                             maxUnalignDimension;
        
        int***                          profile;
        int*                            xProfile;
        int*                            yProfile;
        double**                        dp;
        double**                        scoring;
        
        int**                           indelMatrix;
        int*                            profileNumber;
        int**                           tempProfile;
        
        Node**                          lftChildren;
        Node**                          rhtChildren;
        int*                            sorter;
        
        int*                            possibles;
        int*                            state;
        int*                            pathPos;
        int*                            pathNewPos;
        int*                            pathMask;
        int*                            pathFinalPos;
        std::vector<int>                pathKey;
        std::map<std::vector<int>, int> dpTable;
        
        enum StateLabels { freeToUse, possible, edgeUsed, used };
};

#endif 
