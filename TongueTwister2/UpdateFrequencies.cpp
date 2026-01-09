#include <map>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateFrequencies.hpp"

double UpdateFrequencies::minVal = 0.00001;



UpdateFrequencies::UpdateFrequencies(Model* m, RandomVariable* r, ParameterFrequencies* p) : Update(m, r), myParm(p) {

    numStates = myParm->getNumStates();

    updateNames[0] = "Stationary Frequencies: Dirichlet k=1";
    updateNames[1] = "Stationary Frequencies: Mass Transfer";
    updateNames[2] = "Stationary Frequencies: ALR MVN";
    updateNames[3] = "Stationary Frequencies: Prior";
    tuningValues[0] = 100.0;
    tuningValues[1] = 300.0;
    tuningValues[2] = 0.001;
    tuningValues[3] = 0.0;
    for (size_t i=0; i<4; i++)
        updateHashes[i] = hashUpdateName(updateNames[i]);
}

uint64_t UpdateFrequencies::getUpdateId(void) {

    return updateHashes[lastUpdate];
}

void UpdateFrequencies::setDependants(void) {

    clearDependencyFlags();
    
    updatedParameter      = myParm;
    allTiprobsNeedUpdate  = true;   // Model changed → all matrices need recomputation
    singleBranchChanged   = false;
    
    // Rate matrix needs update for GTR, but not for F81
    if (model->getRateMatrix() == nullptr)
        rateMatrixNeedsUpdate = false;
    else
        rateMatrixNeedsUpdate = true;
}

double UpdateFrequencies::update(void) {

    double u = rng->uniformRv();
    if (u < 0.5)
        return updateMassTransfer(tuningValues[1]);
    else
        return updateDirichlet(1);
}

double UpdateFrequencies::updateDirichlet(int k) {
    
    lastUpdate = 0;
    
    std::vector<double>& freqs = myParm->getFrequencies();
    double bounce = 1.0;
    double lnP = 0.0;
    if (k == 1)
        {
        // resize vectors for move
        oldValues.resize(2, 0.0);
        newValues.resize(2, 0.0);
        alphaForward.resize(2, 0.0);
        alphaReverse.resize(2, 0.0);
        
        // parameterize dirichlet for forward move
        double alpha0 = tuningValues[0];
        int indexToUpdate = static_cast<int>(rng->uniformRv()*numStates);
        oldValues[0] = freqs[indexToUpdate];
        oldValues[1] = 1.0 - oldValues[0];
        alphaForward[0] = oldValues[0] * alpha0 + bounce;
        alphaForward[1] = oldValues[1] * alpha0 + bounce;
        if (alphaForward[0] < 0.0)
            Msg::error("Negative alpha[0] in k=1 proposal for frequencies");
        if (alphaForward[1] < 0.0)
            Msg::error("Negative alpha[1] in k=1 proposal for frequencies");
        
        // propose new values
        bool err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
        if (err == true)
            Msg::error("Problem proposing dirichlet values in k=1 for frequencies");
        double factor = newValues[1] / oldValues[1];
        for (size_t i=0; i<numStates; i++)
            freqs[i] *= factor;
        freqs[indexToUpdate] = newValues[0];
        Probability::Helper::normalize(newValues, minVal);

        // parameterize dirichlet for reverse move
        alphaReverse[0] = newValues[0] * alpha0 + bounce;
        alphaReverse[1] = newValues[1] * alpha0 + bounce;

        // calculate proposal probability
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < static_cast<int>(numStates))
        {
        double alpha0 = 1000.0;
        
        std::vector<int> indicesToUpdate(numStates);
        for (size_t i=0; i<numStates; i++)
            indicesToUpdate[i] = static_cast<int>(i);
        for (int i=0; i<k; i++)
            {
            int idx = static_cast<int>(rng->uniformRv() * (numStates-i));
            indicesToUpdate[idx] = indicesToUpdate[numStates-i-1];
            indicesToUpdate.pop_back();
            }
        if (indicesToUpdate.size() != (size_t)k)
            Msg::error("Problem selecting indices to update");
        
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
            
        oldValues.resize(k+1, 0.0);
        newValues.resize(k+1, 0.0);
        alphaForward.resize(k+1, 0.0);
        alphaReverse.resize(k+1, 0.0);
        
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += freqs[it->first];
            else
                oldValues[k] += freqs[i];
            }
        
        for (int i=0; i<k+1; i++)
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
        
        // draw a new value for the reduced vector
        bool err = false;
        do
            {
            err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
            if (err == true)
                Msg::warning("Trying to recover from Dirichlet random variable error in state frequencies update (k=intermediate)");
            } while (err == true);
        Probability::Helper::normalize(newValues, minVal);
        
        // fill in the Dirichlet parameters for the reverse probability calculations
        for (int i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;
        
        // fill in the full vector
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                freqs[i] = newValues[it->second];
            else
                freqs[i] *= factor;
            }
        
        // Hastings ratio
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - k - 1) * log(factor); // Jacobian
        }
    else
        {
        // update all of the rates
        double alpha0 = 100000.0;
        oldValues.resize(numStates, 0.0);
        newValues.resize(numStates, 0.0);
        alphaForward.resize(numStates, 0.0);
        alphaReverse.resize(numStates, 0.0);

        for (size_t i=0; i<numStates; i++)
            {
            oldValues[i] = freqs[i];
            alphaForward[i] = oldValues[i] * alpha0 + bounce;
            }
        
        bool err = false;
        do
            {
            err = Probability::Dirichlet::rv(rng, alphaForward, newValues);
            if (err == true)
                Msg::warning("Trying to recover from Dirichlet random variable error in state frequencies update (k=N)");
            } while (err == true);
        Probability::Helper::normalize(newValues, minVal);
        
        for (size_t i=0; i<numStates; i++)
            alphaReverse[i] = newValues[i] * alpha0 + bounce;

        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        
        freqs = newValues;
        }

    setDependants();
    
    return lnP;
}

double UpdateFrequencies::updateMassTransfer(double alpha0) {

    lastUpdate = 1;

    std::vector<double>& newFreqs = myParm->getFrequencies(); // 0
    std::vector<double>& oldFreqs = myParm->getFrequencies(1);

    const size_t K = newFreqs.size();
    if (newFreqs.size() != oldFreqs.size())
        Msg::error("Stationary frequencies mass transfer: unequal vector lengths");
    if (K < 2)
        Msg::error("Stationary frequencies mass transfer: k < 2");
    if (!(alpha0 > 0.0))
        Msg::error("Stationary frequencies mass transfer: alpha0 < 0");

    // pick a random pair (i,j), i<j
    size_t i = (size_t)std::floor(rng->uniformRv() * (double)K);
    size_t j = (size_t)std::floor(rng->uniformRv() * (double)(K-1));
    if (j >= i)
        j += 1;
    if (j < i)
        std::swap(i, j);

    double xi = oldFreqs[i];
    double xj = oldFreqs[j];
    double s  = xi + xj;

    if (!(s > 0.0))
        Msg::error("Stationary frequencies mass transfer: s <= 0 (should not happen)");

    double u = xi / s;

    // Beta parameters centered at u
    // To avoid a,b becoming too tiny (can cause numerical ugliness), clamp u away from 0/1 a hair.
    const double eps = 1e-15;
    double uClamped = std::min(1.0 - eps, std::max(eps, u));

    double a = alpha0 * uClamped;
    double b = alpha0 * (1.0 - uClamped);
    double uProp = Probability::Beta::rv(rng, a, b);

    // construct proposal
    newFreqs[i] = s * uProp;
    newFreqs[j] = s * (1.0 - uProp);

    // (tiny numerical safety) renormalize in case of roundoff
    Probability::Helper::normalize(newFreqs, minVal);

    // reverse proposal is Beta centered at u_prop (from the proposed state)
    double uRev = newFreqs[i] / (newFreqs[i] + newFreqs[j]); // should equal u_prop up to rounding
    double uRevClamped = std::min(1.0 - eps, std::max(eps, uRev));

    double aRev = alpha0 * uRevClamped;
    double bRev = alpha0 * (1.0 - uRevClamped);

    // Hastings ratio (reverse - forward)
    // Forward density: u_prop ~ Beta(a,b)
    // Reverse density: u ~ Beta(a_rev, b_rev)
    double log_q_fwd = Probability::Beta::lnPdf(a, b, uProp);
    double log_q_rev = Probability::Beta::lnPdf(aRev, bRev, u);
    double lnHastings = log_q_rev - log_q_fwd;

    if (!Probability::Helper::isValidSimplex(newFreqs, 1e-10))
        Msg::error("Stationary frequencies mass transfer: invalid simplex proposed");

    setDependants();

    return lnHastings;
}

double UpdateFrequencies::update(double power) {

    if (rng->uniformRv() < priorSampleProb(power))
        return updateFromPrior();
    return update();
}

double UpdateFrequencies::updateFromPrior(void) {

    lastUpdate = 3;

    // draw from the prior distribution, which is a flat Dirichlet distribution
    std::vector<double>& freqs = myParm->getFrequencies();
    std::vector<double>& alpha = myParm->getAlpha();
    Probability::Dirichlet::rv(rng, alpha, freqs);

    setDependants();

    return 0.0;
}
