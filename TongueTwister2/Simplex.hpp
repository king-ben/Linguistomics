#ifndef Simplex_hpp
#define Simplex_hpp

#include <vector>
class RandomVariable;



namespace Simplex {

    double  updateALRMVN(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal);
    double  updateALRMVN(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal, size_t blockSize);
    double  updateCenteredDirichlet(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
    double  updateCenteredDirichlet(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal);
    double  updateMassTransfer(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
    double  updateFromPrior(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, std::vector<double>& alpha);
}

#endif
