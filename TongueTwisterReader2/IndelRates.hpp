#ifndef IndelRates_hpp
#define IndelRates_hpp

#include <vector>
#include "json.hpp"
class McmcSummary;



class IndelRates {

    public:
                                            IndelRates(void) = delete;
                                            IndelRates(McmcSummary& samples, double bf);
        size_t                              numSamples(void) { return values.size(); }
        std::pair<float,float>&             getMean(void) { return mean; }
        std::pair<float,float>&             getLowerCI(void) { return lowerCI; }
        std::pair<float,float>&             getUpperCI(void) { return upperCI; }
        void                                print(void);
        nlohmann::json                      toJson(void);
        
    private:
        double                              burnFraction;
        std::vector<std::pair<float,float>> values;
        std::pair<float,float>              mean;
        std::pair<float,float>              variance;
        std::pair<float,float>              lowerCI;
        std::pair<float,float>              upperCI;
};

#endif
