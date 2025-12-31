#ifndef UpdateAlignment_hpp
#define UpdateAlignment_hpp

#include "Update.hpp"
class Model;
class ParameterAlignment;
class ParameterFrequencies;
class ParameterTree;
class TransitionProbabilities;



class UpdateAlignment : public Update {

    public:
                                        UpdateAlignment(void) = delete;
                                        UpdateAlignment(Model* m, RandomVariable* r, ParameterAlignment* p);
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
                                        // Returns the log proposal ratio for the alignment changes.
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

        static constexpr double         basis = 1.5; 
        static constexpr double         gapPenalty = -5.0;
        static constexpr double         extensionProb = 0.5;
        
                                        // pointers accessed frequently (all 8-byte aligned)
        ParameterAlignment*             myParm;
        ParameterFrequencies*           freqsParm;
        ParameterTree*                  treeParm;
        TransitionProbabilities*        tiProbs;
        Tree*                           tree;
        
                                        // array pointers used in inner loops
        int***                          profile;
        int*                            xProfile;
        int*                            yProfile;
        double**                        dp;
        double**                        scoring;
        int**                           indelMatrix;
        int*                            profileNumber;
        int**                           tempProfile;
        Node**                          lftDescendants;
        Node**                          rhtDescendants;
        int*                            sorter;
        int*                            possibles;
        int*                            state;
        int*                            pathPos;
        int*                            pathNewPos;
        int*                            pathMask;
        int*                            pathFinalPos;
        
                                        // containers
        std::vector<int>                pathKey;
        std::map<std::vector<int>, int> dpTable;
        
                                        // int values grouped together (each 4 bytes, pack efficiently)
        int                             numTaxa;
        int                             numNodes;
        int                             numStates;
        int                             gapCode;
        int                             maxLength;
        int                             maxUnalignDimension;
        
                                        // unsigned at end to avoid mid-struct padding
        unsigned                        taxonMask;
                
        enum StateLabels { freeToUse, possible, edgeUsed, used };
};

#endif 
