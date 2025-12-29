#ifndef AlnMatrix_hpp
#define AlnMatrix_hpp

#include "Container.hpp"



class AlnMatrix : public IntMatrix {

    public:
                    AlnMatrix(void) = delete;
                    AlnMatrix(int nt, int maxNs);
        AlnMatrix&  operator=(const AlnMatrix& rhs);
        bool        operator==(const AlnMatrix& m) const;
        int         getNumSites(void) { return numSites; }
        int         getNumTaxa(void) { return (int)numRows; }
        void        print(void);
        void        setNumSites(int ns);
    
    private:
        int         numSites;
        int         maxNumSites;
};

#endif
