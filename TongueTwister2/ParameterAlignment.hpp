#ifndef ParameterAlignment_hpp
#define ParameterAlignment_hpp

#include "json.hpp"
#include "Parameter.hpp"
class Alignment;
class ParameterIndelRates;



class ParameterAlignment: public Parameter {

    public:
                                        ParameterAlignment(void) = delete;
                                        ParameterAlignment(Model* m, RandomVariable* r, std::string n);
                                       ~ParameterAlignment(void);
        Alignment*                      getAlignment(void) { return alignment[0]; }
        Alignment*                      getAlignment(size_t idx) { return alignment[idx]; }
        size_t                          getNumStates(void) { return numStates; }
        size_t                          getNumTaxa(void);
        const unsigned&                 getTaxonMask(void) const { return taxonMask; }
        const std::vector<std::string>& getTaxonNames(void) const { return taxonNames; }
        bool                            initializeFromJson(nlohmann::json j, size_t ns, std::vector<std::string> ctl);
        void                            keep(void);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            restore(void);
        void                            setCanonicalTaxonList(std::vector<std::string> cl) { canonicalTaxonList = cl; }
        
    private:
        size_t                          calculateMaximumAlignmentLength(nlohmann::json j);
        size_t                          longestNameLength(void);
        
        ParameterIndelRates*            indelRates;                // 8 bytes
        Alignment*                      alignment[2];              // 8 bytes
        size_t                          gapCode;                   // 8 bytes
        size_t                          numStates;                 // 8 bytes
        unsigned                        taxonMask;                 // 4 bytes
        std::vector<std::string>        canonicalTaxonList;        // containers (24 bytes each due to std::vector guts)
        std::vector<std::string>        taxonNames;
};

#endif
