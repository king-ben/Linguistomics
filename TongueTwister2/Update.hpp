#ifndef Update_hpp
#define Update_hpp

#include <string>
#include <vector>
class Model;
class Parameter;
class RandomVariable;



class Update {

    public:
                                        Update(void) = delete;
                                        Update(Model* m, RandomVariable* r);
        virtual                        ~Update(void) { }
        
                                        // accessors for dependency flags
        Parameter*                      getUpdatedParameter(void) { return updatedParameter; }
        bool                            getRateMatrixNeedsUpdate(void) { return rateMatrixNeedsUpdate; }
        
                                        // transition probability update types
        bool                            getAllTiprobsNeedUpdate(void) { return allTiprobsNeedUpdate; }
        bool                            getSingleBranchChanged(void) { return singleBranchChanged; }
        double                          getChangedBranchLength(void) { return changedBranchLength; }
        
                                        // For updates that modify multiple parameters (e.g., UpdateTopology
                                        // modifies the tree AND all alignments). These additional parameters
                                        // need to be kept/restored along with the primary updatedParameter.
        virtual std::vector<Parameter*> getAdditionalModifiedParameters(void) { return {}; }
        
        virtual std::string             getUpdateName(void) = 0;
        virtual std::string             parameterType(void) = 0;
        virtual void                    setDependants(void) = 0;
        virtual double                  update(void) = 0;
        virtual double                  update(double power) = 0;
        virtual double                  updateFromPrior(void) = 0;
        
    protected:
        void                            clearDependencyFlags(void);
        double                          priorSampleProb(double power);
        
        Model*                          model;
        RandomVariable*                 rng;
        Parameter*                      updatedParameter;
        
                                        // dependency flags
        bool                            rateMatrixNeedsUpdate;
        bool                            allTiprobsNeedUpdate;     // model changed -> recompute all matrices
        bool                            singleBranchChanged;      // one branch length changed -> recompute one matrix
        double                          changedBranchLength;      // the new branch length (if singleBranchChanged)
};

#endif
