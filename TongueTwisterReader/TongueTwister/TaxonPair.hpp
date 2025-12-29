#ifndef TaxonPair_hpp
#define TaxonPair_hpp

#include <string>


class TaxonPair {

    public:
                        TaxonPair(void) = delete;
                        TaxonPair(const TaxonPair& p);
                        TaxonPair(std::string t1, std::string t2);
        std::string     getTaxon1(void) { return taxon1; }
        std::string     getTaxon2(void) { return taxon2; }
        std::string     getTaxon1(void) const { return taxon1; }
        std::string     getTaxon2(void) const { return taxon2; }
    
    private:
        std::string     taxon1;
        std::string     taxon2;
        
};

struct CompTaxonPair {

    bool operator()(const TaxonPair& t1, const TaxonPair& t2) const {
        
        if (t1.getTaxon1() < t2.getTaxon1())
            return true;
        else if (t1.getTaxon1() == t2.getTaxon1())
            {
            if (t1.getTaxon2() < t2.getTaxon2())
                return true;
            }
        return false;
        }
};

#endif
