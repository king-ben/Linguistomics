#ifndef Simplex_hpp
#define Simplex_hpp

#include <vector>
class RandomVariable;



namespace Simplex {

    double  updateALRMVN(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal);
    double  updateALRMVN(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal, size_t blockSize);
    double  updateCenteredDirichlet(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
    double  updateCenteredDirichlet(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal);
    double  updateMassTransfer(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
    double  updateFromPrior(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, std::vector<double>& alpha);
}

#endif
