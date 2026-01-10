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
        
        virtual double                  getTuningParameter(void) = 0;      
                                        // accessors for dependency flags
        Parameter*                      getUpdatedParameter(void) { return updatedParameter; }
        virtual uint64_t                getUpdateId(void) = 0;
        bool                            getRateMatrixNeedsUpdate(void) { return rateMatrixNeedsUpdate; }
        
                                        // transition probability update types
        bool                            getAllTiprobsNeedUpdate(void) { return allTiprobsNeedUpdate; }
        bool                            getSingleBranchChanged(void) { return singleBranchChanged; }
        double                          getChangedBranchLength(void) { return changedBranchLength; }
        
                                        /* For updates that modify multiple parameters (e.g., UpdateTopology
                                           modifies the tree AND all alignments). These additional parameters
                                           need to be kept/restored along with the primary updatedParameter. */
        virtual std::vector<Parameter*> getAdditionalModifiedParameters(void) { return {}; }
        
        virtual std::string             getUpdateName(void) = 0;
        virtual std::string             parameterType(void) = 0;
        virtual void                    setDependants(void) = 0;
        virtual double                  update(void) = 0;
        virtual double                  update(double power) = 0;
        
    protected:
        void                            clearDependencyFlags(void);
        uint64_t                        hashUpdateName(const std::string& name);
        uint64_t                        load_tail_le(const unsigned char* p, size_t n);
        uint64_t                        load64_le(const unsigned char* p);
        double                          priorSampleProb(double power);
        uint64_t                        rotl64(uint64_t x, int b);
        uint64_t                        siphash24(const unsigned char* data, size_t len, uint64_t k0, uint64_t k1);
        void                            sip_round(uint64_t& v0, uint64_t& v1, uint64_t& v2, uint64_t& v3);
        
        Model*                          model;
        RandomVariable*                 rng;
        Parameter*                      updatedParameter;
        uint64_t                        updateId;
        double                          tuningParameter;
        
        double                          changedBranchLength;        // the new branch length (if singleBranchChanged)
        bool                            rateMatrixNeedsUpdate;
        bool                            allTiprobsNeedUpdate;       // model changed -> recompute all matrices
        bool                            singleBranchChanged;        // one branch length changed -> recompute one matrix
};

#endif
