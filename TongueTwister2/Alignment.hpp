#ifndef AlnMatrix_hpp
#define AlnMatrix_hpp

#include "Container.hpp"



class Alignment : public IntMatrix {

    public:
                    Alignment(void) = delete;
                    Alignment(const Alignment& a);
                    Alignment(size_t nt, size_t maxNs);
        Alignment&  operator=(const Alignment& rhs);
        bool        operator==(const Alignment& m) const;
        size_t      getMaximumNumberOfSites(void) { return maxNumSites; }
        size_t      getNumSites(void) { return numSites; }
        size_t      getNumTaxa(void) { return numTaxa; }
        void        print(void);
        void        setNumSites(size_t ns);
    
    private:
        size_t      numTaxa;
        size_t      numSites;
        size_t      maxNumSites;
};

#endif
