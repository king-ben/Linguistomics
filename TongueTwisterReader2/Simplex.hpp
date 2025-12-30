#ifndef Simplex_hpp
#define Simplex_hpp

#include <cstdlib>
#include <map>
#include <regex>
#include <vector>
typedef std::vector<std::vector<float>>  SampleVector;

struct StatInfo {

    float           mean;
    float           priorMean;
    float           variance;
    float           lowerCI;
    float           upperCI;
    size_t          index1;
    size_t          index2;
    std::string     name;
};



class Simplex {

    public:
                                        Simplex(void) = delete;
                                        Simplex(double bf);
        float                           getKullbackLiebler(size_t idx) { return kl[idx]; }
        virtual size_t                  getFirstIndex(size_t idx) = 0;
        virtual size_t                  getSecondIndex(size_t idx) = 0;
        float                           getMean(size_t idx) { return mean[idx]; }
        float                           getPriorMean(size_t idx) { return priorMean[idx]; }
        float                           getVariance(size_t idx) { return variance[idx]; }
        float                           getLowerCI(size_t idx) { return lowerCI[idx]; }
        float                           getUpperCI(size_t idx) { return upperCI[idx]; }
        const SampleVector&             getValues(void) const { return values; }
        size_t                          numSamples(void) { return values[0].size(); }
        virtual void                    print(void) = 0;
        size_t                          size(void) { return values.size(); }
        size_t                          size(void) const { return values.size(); }
        std::map<float,StatInfo>        sortByKL(void);
        void                            valuesAtIndex(size_t idx, std::vector<float>& vec);
    
    protected:
        std::pair<int,int>              parseBracketedNumbers(const std::string& s);
        double                          burnFraction;
        size_t                          dimension;
        SampleVector                    values;
        std::vector<float>              mean;
        std::vector<float>              priorMean;
        std::vector<float>              variance;
        std::vector<float>              lowerCI;
        std::vector<float>              upperCI;
        std::vector<float>              kl;
};

#endif
