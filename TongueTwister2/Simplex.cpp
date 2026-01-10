#include <set>
#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Simplex.hpp"



static inline void alrForward(const std::vector<double>& x, std::vector<double>& y) {

    const size_t K = x.size();
    y.resize(K - 1);

    double xK = x[K-1];
    for (size_t i=0; i<K-1; i++)
        y[i] = std::log(x[i] / xK);
}

static inline void alrInverse(const std::vector<double>& y, std::vector<double>& x) {

    const size_t d = y.size();
    const size_t K = d + 1;

    // stable inverse: subtract max
    double m = y[0];
    for (size_t i=1; i<d; i++)
        {
        if (y[i] > m)
            m = y[i];
        }

    std::vector<double> e(d, 0.0);
    double sumE = 0.0;
    for (size_t i=0; i<d; i++)
        {
        e[i] = std::exp(y[i] - m);
        sumE += e[i];
        }

    double expNegM = std::exp(-m);
    double denom = expNegM + sumE;

    x.assign(K, 0.0);
    for (size_t i=0; i<d; i++)
        x[i] = e[i] / denom;
    x[K-1] = expNegM / denom;
}

static inline bool isValidSimplex(const std::vector<double>& x, double eps) {

    double s = 0.0;
    for (size_t i=0; i<x.size(); i++)
        {
        if (!(x[i] > 0.0))
            return false;
        s += x[i];
        }
    return std::fabs(s - 1.0) < eps;
}

static inline double sumLogVec(const std::vector<double>& x) {

    double s = 0.0;
    for (size_t i=0; i<x.size(); i++)
        {
        if (x[i] <= 0.0)
            return -1e300;
        s += std::log(x[i]);
        }
    return s;
}

double Simplex::updateALRMVN(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal) {

    const size_t K = oldVec.size();
    if (oldVec.size() != newVec.size())
        Msg::error("updateSimplexALRMVN: unequal vector lengths");
    if (K < 2)
        Msg::error("updateSimplexALRMVN: K < 2");
    if (!(sigma > 0.0))
        Msg::error("updateSimplexALRMVN: sigma <= 0");
    if (!isValidSimplex(oldVec, 1e-10))
        Msg::error("updateSimplexALRMVN: oldVec not valid simplex");

    // forward ALR transform
    std::vector<double> y;
    alrForward(oldVec, y);

    // MVN random walk in ALR space (isotropic sigma^2 I)
    for (size_t i=0; i<y.size(); i++)
        y[i] += sigma * Probability::Normal::rv(rng);

    // inverse transform
    alrInverse(y, newVec);

    // For this kernel, do NOT "normalize-with-floor" (that changes the mapping).
    // Instead, hard-reject if any component violates minVal.
    if (minVal > 0.0)
        {
        for (size_t i=0; i<newVec.size(); i++)
            {
            if (newVec[i] < minVal)
                {
                newVec = oldVec;
                return -1e300;
                }
            }
        }

    if (!isValidSimplex(newVec, 1e-10))
        Msg::error("updateSimplexALRMVN: newVec not valid simplex");

    // Hastings/Jacobian correction:
    // proposal symmetric in ALR space => only Jacobian ratio remains
    return sumLogVec(newVec) - sumLogVec(oldVec);
}

double Simplex::updateFromPrior(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, std::vector<double>& alpha) {

    Probability::Dirichlet::rv(rng, alpha, newVec);
    return Probability::Dirichlet::lnPdf(alpha, oldVec) - Probability::Dirichlet::lnPdf(alpha, newVec);
}

double Simplex::updateMassTransfer(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal) {

    const size_t K = oldVec.size();
    if (oldVec.size() != newVec.size())
        Msg::error("proposeSimplex_MassTransfer: unequal vector lengths");
    if (K < 2)
        Msg::error("proposeSimplex_MassTransfer: K < 2");
    if (!(alpha0 > 0.0))
        Msg::error("proposeSimplex_MassTransfer: alpha0 <= 0");
    if (!isValidSimplex(oldVec, 1e-10))
        Msg::error("proposeSimplex_MassTransfer: x_curr not valid simplex");

    // pick a random pair (i,j), i<j
    size_t i = (size_t)std::floor(rng->uniformRv() * (double)K);
    size_t j = (size_t)std::floor(rng->uniformRv() * (double)(K-1));
    if (j >= i)
        j += 1;
    if (j < i)
        std::swap(i, j);

    double xi = oldVec[i];
    double xj = oldVec[j];
    double s  = xi + xj;

    if (!(s > 0.0))
        Msg::error("proposeSimplex_MassTransfer: s <= 0 (should not happen)");

    double u = xi / s;

    // Beta parameters centered at u
    // To avoid a,b becoming too tiny (can cause numerical ugliness), clamp u away from 0/1 a hair.
    const double eps = 1e-15;
    double uClamped = std::min(1.0 - eps, std::max(eps, u));

    double a = alpha0 * uClamped;
    double b = alpha0 * (1.0 - uClamped);

    double uProp = Probability::Beta::rv(rng, a, b);

    // construct proposal
    newVec[i] = s * uProp;
    newVec[j] = s * (1.0 - uProp);

    // (tiny numerical safety) renormalize in case of roundoff
    Probability::Helper::normalize(newVec, minVal);

    // reverse proposal is Beta centered at u_prop (from the proposed state)
    double uRev = newVec[i] / (newVec[i] + newVec[j]); // should equal u_prop up to rounding
    double uRevClamped = std::min(1.0 - eps, std::max(eps, uRev));

    double aRev = alpha0 * uRevClamped;
    double bRev = alpha0 * (1.0 - uRevClamped);

    // Hastings ratio (reverse - forward)
    // Forward density: u_prop ~ Beta(a,b)
    // Reverse density: u ~ Beta(a_rev, b_rev)
    double log_q_fwd = Probability::Beta::lnPdf(a, b, uProp);
    double log_q_rev = Probability::Beta::lnPdf(aRev, bRev, u);
    double lnHastings = log_q_rev - log_q_fwd;

    if (!isValidSimplex(newVec, 1e-10))
        Msg::error("proposeSimplex_MassTransfer: x_prop invalid");
    return lnHastings;
}

double Simplex::updateCenteredDirichlet(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal) {
    
    std::vector<double> alphaForward(oldVec.size());
    for (size_t i=0, n=alphaForward.size(); i<n; i++)
        alphaForward[i] = oldVec[i] * alpha0 + 1.0;
        
    Probability::Dirichlet::rv(rng, alphaForward, newVec);
    Probability::Helper::normalize(newVec, minVal);
    
    std::vector<double> alphaReverse(newVec.size());
    for (size_t i=0, n=alphaReverse.size(); i<n; i++)
        alphaReverse[i] = newVec[i] * alpha0 + 1.0;

    double lnForwardProb = Probability::Dirichlet::lnPdf(alphaReverse, oldVec) - Probability::Dirichlet::lnPdf(alphaForward, newVec);
    return lnForwardProb;
}

double Simplex::updateCenteredDirichlet(RandomVariable* rng, std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal) {
    
    // choose which elements to update
    size_t n = oldVec.size();
    if (n == k)
        Msg::error("Cannot update that many simplex variables!!!");
        
    std::set<size_t> indices;
    while (indices.size() < k)
        {
        size_t idx = (size_t)(rng->uniformRv() * n);
        indices.insert(idx);
        }
    
    // fill in old values for the selected indices plus the "remainder" element
    std::vector<double> oldValues(k+1);
    double sum = 0.0;
    int i = 0;
    for (size_t idx : indices)
        {
        double x = oldVec[idx];
        oldValues[i++] = x;
        sum += x;
        }
    oldValues[i] = 1.0 - sum;  // remainder element
    
    // construct the forward Dirichlet parameters
    std::vector<double> alphaForward(k + 1);
    for (size_t j=0, m=alphaForward.size(); j<m; j++)
        alphaForward[j] = oldValues[j] * alpha0 + 1.0;
        
    // draw new values from Dirichlet
    std::vector<double> newValues(k+1);
    Probability::Dirichlet::rv(rng, alphaForward, newValues);
    
    // FIX 1: normalize newValues (the k+1 element vector), NOT newVec (the full n-element vector)
    Probability::Helper::normalize(newValues, minVal);
    
    // construct the reverse Dirichlet parameters
    std::vector<double> alphaReverse(k + 1);
    for (size_t j=0, m=alphaReverse.size(); j<m; j++)
        alphaReverse[j] = newValues[j] * alpha0 + 1.0;
        
    // fill in the full newVec from the updated (newValues) vector
    // The non-selected elements are scaled by the ratio of remainders
    double factor = newValues[k] / oldValues[k];
    
    // FIX 2: Initialize newVec from oldVec explicitly, then scale
    for (size_t j=0, m=newVec.size(); j<m; j++)
        newVec[j] = oldVec[j] * factor;
    
    // overwrite the selected indices with their new values
    i = 0;
    for (size_t idx : indices)
        newVec[idx] = newValues[i++];

    // FIX 3: Compute Hastings ratio using oldValues and newValues (both size k+1),
    // NOT oldVec and newVec (which are size n and would cause dimension mismatch)
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
    lnProposalProb += (n - k - 1) * log(factor);
    return lnProposalProb;
}
