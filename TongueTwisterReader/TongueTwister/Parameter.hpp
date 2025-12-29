#ifndef Parameter_hpp
#define Parameter_hpp

#include <string>
class Model;
class RandomVariable;



class Parameter {

    public:
                                    Parameter(void) = delete;
                                    Parameter(RandomVariable* r, Model* m, std::string n);
        virtual                    ~Parameter(void) { }
        virtual void                accept(void) = 0;
        virtual void                fillParameterValues(double* x, int& start, int maxNumValues) = 0;
        virtual std::string         getJsonString(void) = 0;
        virtual std::string         getHeader(void) = 0;
        virtual int                 getNumValues(void) = 0;
        int                         getParameterId(void) { return parameterId; }
        double                      getProposalProbability(void) { return proposalProbability; }
        virtual std::string         getString(void) = 0;
        virtual char*               getCString(void) = 0;
        int                         getParameterStringLength(void) { return parmStrLen; }
        bool                        getUpdateChangesRateMatrix(void) { return updateChangesRateMatrix; }
        bool                        getUpdateChangesTransitionProbabilities(void) { return updateChangesTransitionProbabilities; }
        virtual std::string         getUpdateName(int idx) = 0;
        virtual double              lnPriorProbability(void) = 0;
        virtual void                print(void) = 0;
        virtual void                reject(void) = 0;
        void                        setParameterId(int x) { parameterId = x; }
        void                        setProposalProbability(double x) { proposalProbability = x; }
        void                        setUpdateChangesRateMatrix(bool tf) { updateChangesRateMatrix = tf; }
        void                        setUpdateChangesTransitionProbabilities(bool tf) { updateChangesTransitionProbabilities = tf; }
        virtual double              update(int iter) = 0;
        std::pair<Parameter*,int>   lastUpdateType;
        std::string                 parmName;
        
    protected:
        Model*                      modelPtr;
        RandomVariable*             rv;
        int                         parameterId;
        double                      proposalProbability;
        bool                        updateChangesRateMatrix;
        bool                        updateChangesTransitionProbabilities;
        char*                       parmStr;
        int                         parmStrLen;
};

#endif
