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

double Simplex::updateALRMVN(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal) {

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

    // MVN random walk in ALR space (sigma^2 I)
    for (size_t i=0; i<y.size(); i++)
        y[i] += sigma * Probability::Normal::rv(rng);

    // inverse transform
    alrInverse(y, newVec);

    /* Do NOT "normalize-with-floor" (that changes the mapping).
       Instead, hard-reject if any component violates minVal. */
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

    /* Hastings/Jacobian correction:
       proposal symmetric in ALR space -> only Jacobian ratio remains */
    return sumLogVec(newVec) - sumLogVec(oldVec);
}

double Simplex::updateALRMVN(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal, size_t blockSize) {

    const size_t K = oldVec.size();
    if (oldVec.size() != newVec.size())
        Msg::error("updateSimplexALRMVN_blocked: unequal vector lengths");
    if (K < 2)
        Msg::error("updateSimplexALRMVN_blocked: K < 2");
    if (!(sigma > 0.0))
        Msg::error("updateSimplexALRMVN_blocked: sigma <= 0");
    if (blockSize < 2)
        Msg::error("updateSimplexALRMVN_blocked: blockSize < 2");
    if (blockSize > K)
        Msg::error("updateSimplexALRMVN_blocked: blockSize > K");
    if (!isValidSimplex(oldVec, 1e-10))
        Msg::error("updateSimplexALRMVN_blocked: oldVec not valid simplex");

    // default: start with no change
    //newVec = oldVec; // not needed because both vectors start every MCMC cycle identical

    // choose a random block of indices 
    // We want a set of distinct indices of size blockSize (Fisher-Yates partial shuffle)
    std::vector<size_t> idx(K);
    for (size_t i=0; i<K; i++)
        idx[i] = i;

    for (size_t i=0; i<blockSize; i++)
        {
        size_t j = i + (size_t)std::floor( rng->uniformRv() * (double)(K - i) );
        if (j >= K)
            j = K - 1;
        std::swap(idx[i], idx[j]);
        }

    idx.resize(blockSize);

    // extract block mass and block proportions
    double blockMass = 0.0;
    for (size_t t=0; t<blockSize; t++)
        blockMass += oldVec[idx[t]];

    if (!(blockMass > 0.0))
        {
        // should never happen for a valid simplex, but be safe
        return -1e300;
        }

    // proportions p live on an m-simplex (m = blockSize)
    std::vector<double> p_old(blockSize, 0.0);
    for (size_t t=0; t<blockSize; t++)
        p_old[t] = oldVec[idx[t]] / blockMass;

    // if you enforce a global minVal, then inside the block the implied minimum on p is minVal / blockMass.
    double pMin = 0.0;
    if (minVal > 0.0)
        pMin = minVal / blockMass;

    // propose new block proportions with your existing ALR-MVN code
    std::vector<double> p_new(blockSize, 0.0);

    /* This returns the Hastings correction for the move in p-space.
       Because x_block = blockMass * p, the constant blockMass cancels in Jacobians and in proposal densities.
       So the Hastings correction for the overall x-move is IDENTICAL to the one for p. */
    double logH = updateALRMVN(rng, p_old, p_new, sigma, pMin);

    // if updateALRMVN hard-rejected (returns LOG_ZERO), propagate rejection.
    if (logH <= -1e250)
        {
        newVec = oldVec;
        return logH;
        }

    // write back block and validate global minVal if requested
    for (size_t t=0; t<blockSize; t++)
        newVec[idx[t]] = blockMass * p_new[t];

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
        Msg::error("updateSimplexALRMVN_blocked: newVec not valid simplex");

    return logH;
}

double Simplex::updateFromPrior(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, std::vector<double>& alpha) {

    Probability::Dirichlet::rv(rng, alpha, newVec);
    return Probability::Dirichlet::lnPdf(alpha, oldVec) - Probability::Dirichlet::lnPdf(alpha, newVec);
}

double Simplex::updateMassTransfer(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal) {

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

    // renormalize in case of roundoff
    Probability::Helper::normalize(newVec, minVal);

    // reverse proposal is Beta centered at u_prop (from the proposed state)
    double uRev = newVec[i] / (newVec[i] + newVec[j]); // should equal u_prop up to rounding
    double uRevClamped = std::min(1.0 - eps, std::max(eps, uRev));

    double aRev = alpha0 * uRevClamped;
    double bRev = alpha0 * (1.0 - uRevClamped);

    /* Hastings ratio (reverse - forward)
       forward density: u_prop ~ Beta(a,b)
       reverse density: u ~ Beta(a_rev, b_rev) */
    double log_q_fwd = Probability::Beta::lnPdf(a, b, uProp);
    double log_q_rev = Probability::Beta::lnPdf(aRev, bRev, u);
    double lnHastings = log_q_rev - log_q_fwd;

    if (!isValidSimplex(newVec, 1e-10))
        Msg::error("proposeSimplex_MassTransfer: x_prop invalid");
    return lnHastings;
}

double Simplex::updateCenteredDirichlet(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal) {
    
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

double Simplex::updateCenteredDirichlet(RandomVariable* rng, const std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal) {
    
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
    
    // normalize newValues (the k+1 element vector)
    Probability::Helper::normalize(newValues, minVal);
    
    // construct the reverse Dirichlet parameters
    std::vector<double> alphaReverse(k + 1);
    for (size_t j=0, m=alphaReverse.size(); j<m; j++)
        alphaReverse[j] = newValues[j] * alpha0 + 1.0;
        
    // fill in the full newVec from the updated (newValues) vector
    // The non-selected elements are scaled by the ratio of remainders
    double factor = newValues[k] / oldValues[k];
    
    // initialize newVec from oldVec explicitly, then scale
    for (size_t j=0, m=newVec.size(); j<m; j++)
        newVec[j] = oldVec[j] * factor;
    
    // overwrite the selected indices with their new values
    i = 0;
    for (size_t idx : indices)
        newVec[idx] = newValues[i++];

    // compute Hastings ratio using oldValues and newValues (both size k+1),
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
    lnProposalProb += (n - k - 1) * log(factor);
    return lnProposalProb;
}
