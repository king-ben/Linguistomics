#ifndef Exchangeabilities_hpp
#define Exchangeabilities_hpp

#include "json.hpp"
#include "Simplex.hpp"
#include "Statistics.hpp"
class McmcSummary;



class Exchangeabilities : public Simplex {

    public:
                                            Exchangeabilities(void) = delete;
                                            Exchangeabilities(McmcSummary& samples, double bf);
        size_t                              getFirstIndex(size_t idx) { return fromTo[idx].first; }
        size_t                              getSecondIndex(size_t idx) { return fromTo[idx].second; }
        size_t                              getNumRates(void) { return dimension; }
        size_t                              getNumRates(void) const { return dimension; }
        double                              getAverageRate(int r1, int r2) const;
        void                                print(void);
        nlohmann::json                      toJson(void);
        
    private:
        std::vector<std::pair<int,int>>     fromTo;
};

#endif
