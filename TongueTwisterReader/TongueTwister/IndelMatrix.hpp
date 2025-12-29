#ifndef IndelMatrix_hpp
#define IndelMatrix_hpp

#include "Container.hpp"



class IndelMatrix : public IntMatrix {

    public:
                    IndelMatrix(void) = delete;
                    IndelMatrix(int nt, int maxNs);
        int         getNumTaxa(void) { return numTaxa; }
        int         getNumSites(void) { return numSites; }
        void        print(void);
        void        setNumSites(int ns);
    
    private:
        int         numTaxa;
        int         numSites;
        int         maxNumSites;
};

#endif
