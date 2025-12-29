#ifndef Statistics_hpp
#define Statistics_hpp

#include <utility>
#include <vector>

struct BetaFitResult {

    float   alpha;
    float   beta;
    int    iterations;
    bool   converged;
};

struct CredibleInterval {
    double lower;
    double upper;
    double median;

     CredibleInterval(void) {
        lower = 0.0;
        upper = 0.0;
        median = 0.0;
    }
   CredibleInterval(double low, double up, double med) {
        lower = low;
        upper = up;
        median = med;
    }
};



namespace Statistics {

    double                      calcKlToPrior(const std::vector<float>& x, size_t n, size_t dimension, BetaFitResult* fit=nullptr);
    CredibleInterval            getCredibleInterval(const std::vector<float>& vals);
    CredibleInterval            getCredibleInterval(const std::vector<float>& vals, float burnFraction);
    std::pair<float,float>      getMeanAndVariance(const std::vector<float>& x);
    std::pair<float,float>      getMeanAndVariance(const std::vector<float>& x, float burnFraction);

    std::pair<double,double>    betaMomInit(const std::vector<float>& xs);
    double                      clamp01Open(double x, double eps=1e-12);
    double                      digamma(double x);
    BetaFitResult               estimatePosteriorBetaMl(const std::vector<float>& x);
    BetaFitResult               fitBetaMl(const std::vector<float>& xRaw, int maxIter=200, double tol=1e-10);
    std::pair<double,double>    flatDirichletMarginalPriorBeta(size_t K);
    double                      klBeta(double a, double b, double c, double d);
    double                      logBeta(double a, double b);
    double                      trigamma(double x);
}

#endif
