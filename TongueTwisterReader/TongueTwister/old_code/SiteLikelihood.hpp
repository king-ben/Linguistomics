#ifndef SiteLikelihood_hpp
#define SiteLikelihood_hpp

#include <iostream>


class SiteLikelihood {

    public:
                    SiteLikelihood(void) = delete;
                    SiteLikelihood(int nn, int ns);
                   ~SiteLikelihood(void);
        int         getNumNodes(void) { return numNodes; }
        int         getNumStates(void) { return numStates; }
        double**    getProbsH(void) { return probsH; }
        double*     getProbsH(int idx) { return probsH[idx]; }
        double**    getProbsI(void) { return probsI; }
        double*     getProbsI(int idx) { return probsI[idx]; }
        void        print(void);
        void        zeroOutH(void);
        void        zeroOutI(void);
        void        zeroOutH(int idx);
        void        zeroOutI(int idx);
    
    private:
        int         numNodes;
        int         numStates,
                    numStates1;
        double**    probsH;
        double**    probsI;
};

#endif
