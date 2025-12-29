#ifndef Exchangeabilities_hpp
#define Exchangeabilities_hpp

#include "Simplex.hpp"
class McmcSummary;



class Exchangeabilities : public Simplex {

    public:
                                            Exchangeabilities(void) = delete;
                                            Exchangeabilities(McmcSummary& samples, double bf);
        size_t                              getFirstIndex(size_t idx) { return fromTo[idx].first; }
        size_t                              getSecondIndex(size_t idx) { return fromTo[idx].second; }
        size_t                              getNumRates(void) { return dimension; }
        void                                print(void);
        
    private:
        std::vector<std::pair<int,int>>     fromTo;
};

#endif
