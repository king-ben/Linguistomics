#include "Msg.hpp"
#include "TaxonPair.hpp"



TaxonPair::TaxonPair(const TaxonPair& p) {

    taxon1 = p.taxon1;
    taxon2 = p.taxon2;
}

TaxonPair::TaxonPair(std::string t1, std::string t2) {

    if (t1 < t2)
        {
        taxon1 = t1;
        taxon2 = t2;
        }
    else if (t2 < t1)
        {
        taxon1 = t2;
        taxon2 = t1;
        }
    else
        {
        Msg::error("Making taxon pair with identical taxon names");
        }
}
