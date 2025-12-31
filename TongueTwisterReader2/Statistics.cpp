#include "Msg.hpp"
#include "Statistics.hpp"



CredibleInterval Statistics::getCredibleInterval(const std::vector<float>& vals) {

    return getCredibleInterval(vals, 0.0);
}

CredibleInterval Statistics::getCredibleInterval(const std::vector<float>& vals, float burnFraction) {

    if (burnFraction < 0.0 || burnFraction > 0.95)
        Msg::error("Burn fraction must be between 0 and 0.95");
        
    CredibleInterval ci(0.0, 0.0, 0.0);
    if (vals.size() == 0)
        {
        }
    else
        {
        size_t startIdx = static_cast<size_t>(vals.size() * burnFraction);
        size_t numSamples = vals.size() - startIdx;
        if (numSamples == 0)
            return ci;
        std::vector<float> tempValues(numSamples);
        for (size_t i=startIdx, k=0, n=vals.size(); i<n; i++)
            tempValues[k++] = vals[i];
        std::sort(tempValues.begin(), tempValues.end());
        size_t size = tempValues.size();
        
        size_t lp = size / 40;
        size_t up = size - lp - 1;

        ci.lower = tempValues[lp];
        ci.upper = tempValues[up];
        if (size % 2 == 0)
            ci.median = (tempValues[size / 2 - 1] + tempValues[size / 2]) / 2.0;
        else
            ci.median = tempValues[size / 2];
        }
        
    return ci;
}

std::pair<float,float> Statistics::getMeanAndVariance(const std::vector<float>& vals) {
    
    return getMeanAndVariance(vals, 0.0);
}

std::pair<float,float> Statistics::getMeanAndVariance(const std::vector<float>& vals, float burnFraction) {

    if (burnFraction < 0.0 || burnFraction > 0.95)
        Msg::error("Burn fraction must be between 0 and 0.95");
    if (vals.size() * (1.0 - burnFraction) < 3)
        Msg::error("Can only calclate statistics for at least three samples");

    size_t startIdx = static_cast<size_t>(vals.size() * burnFraction);

    // Welford's online algorithm for mean and variance
    int n = 0;
    double mean = 0.0;
    double M2 = 0.0;
    
    for (size_t i=startIdx; i<vals.size(); i++)
        {
        n++;
        double x =vals[i];
        double delta = x - mean;
        mean += delta / n;
        double delta2 = x - mean;
        M2 += delta * delta2;
        }

    double variance = 0.0;
    if (n > 1)
        variance = M2 / (n - 1);
        
    std::pair<float,float> meanVariance;
    meanVariance.first = mean;
    meanVariance.second = variance;
    return meanVariance;
}

std::pair<double,double> Statistics::betaMomInit(const std::vector<float>& xs) {

    // method-of-moments init

    if (xs.size() < 2)
        return std::make_pair(1.0, 1.0);

    double mean = 0.0;
    for (size_t i=0; i<xs.size(); i++)
        mean += xs[i];
    mean /= (double)xs.size();

    double var = 0.0;
    for (size_t i=0; i<xs.size(); i++)
        {
        double d = xs[i] - mean;
        var += d * d;
        }
    var /= (double)(xs.size() - 1);

    // basically constant data -> huge concentration
    if (var < 1e-14)
        {
        double conc = 1e6;
        double a = std::max(1e-6, mean * conc);
        double b = std::max(1e-6, (1.0 - mean) * conc);
        return std::make_pair(a, b);
        }

    double common = mean * (1.0 - mean) / var - 1.0;
    if (!std::isfinite(common) || common <= 0.0)
        return std::make_pair(1.0, 1.0);

    double a = std::max(1e-6, mean * common);
    double b = std::max(1e-6, (1.0 - mean) * common);

    return std::make_pair(a, b);
}

double Statistics::calcKlToPrior(const std::vector<float>& x, size_t n, size_t dimension, BetaFitResult* fit) {

    if (dimension < 2)
        throw std::runtime_error("Simplex::calcKlToPrior: dimension must be >= 2");
    if (n >= dimension)
        throw std::runtime_error("Simplex::calcKlToPrior: n must be dimension");

    BetaFitResult r = fitBetaMl(x);
    if (fit != nullptr)
        *fit = r;

    // flat Dirichlet(1,...,1) marginal prior for one component: Beta(1, K-1)
    double c = static_cast<double>(n);
    double d = static_cast<double>(dimension - n);

    return klBeta(r.alpha, r.beta, c, d);
}

double Statistics::clamp01Open(double x, double eps) {

    if (x <= eps)
        return eps;
    if (x >= 1.0 - eps)
        return 1.0 - eps;
    return x;
}

double Statistics::digamma(double x) {

    if (x <= 0.0)
        throw std::runtime_error("digamma: x must be > 0");

    double result = 0.0;

    // recurrence to push x up
    while (x < 8.0)
        {
        result -= 1.0 / x;
        x += 1.0;
        }

    // asymptotic expansion
    double inv  = 1.0 / x;
    double inv2 = inv * inv;

    result += std::log(x) - 0.5 * inv
            - inv2 / 12.0
            + (inv2 * inv2) / 120.0
            - (inv2 * inv2 * inv2) / 252.0
            + (inv2 * inv2 * inv2 * inv2) / 240.0;

    return result;
}

BetaFitResult Statistics::estimatePosteriorBetaMl(const std::vector<float>& x) {

    return fitBetaMl(x);
}

BetaFitResult Statistics::fitBetaMl(const std::vector<float>& xRaw, int maxIter, double tol) {

    if (xRaw.size() == 0)
        throw std::runtime_error("fitBetaMl: empty data");

    std::vector<float> x;
    x.reserve(xRaw.size());
    for (size_t i=0; i<xRaw.size(); i++)
        x.push_back(clamp01Open(xRaw[i]));

    size_t n = x.size();

    double sumLogX   = 0.0;
    double sumLog1mX = 0.0;
    for (size_t i=0; i<n; i++)
        {
        sumLogX   += std::log(x[i]);
        sumLog1mX += std::log(1.0 - x[i]);
        }

    std::pair<double,double> init = betaMomInit(x);
    double a = init.first;
    double b = init.second;

    auto logLik = [&](double aa, double bb) -> double {

        return (double)n * (std::lgamma(aa + bb) - std::lgamma(aa) - std::lgamma(bb))
             + (aa - 1.0) * sumLogX
             + (bb - 1.0) * sumLog1mX;
    };

    double ll = logLik(a, b);
    bool converged = false;
    int it = 0;

    for (it=0; it<maxIter; it++)
        {
        double psiAB = digamma(a + b);
        double psiA  = digamma(a);
        double psiB  = digamma(b);

        // gradient
        double g1 = (double)n * (psiAB - psiA) + sumLogX;
        double g2 = (double)n * (psiAB - psiB) + sumLog1mX;

        // Hessian
        double triAB = trigamma(a + b);
        double h11 = (double)n * (triAB - trigamma(a));
        double h22 = (double)n * (triAB - trigamma(b));
        double h12 = (double)n * triAB;

        double det = h11 * h22 - h12 * h12;
        if (!std::isfinite(det) || fabs(det) < 1e-30)
            break;

        // solve H * delta = g
        double da = ( h22 * g1 - h12 * g2) / det;
        double db = (-h12 * g1 + h11 * g2) / det;

        if (std::max(fabs(da), fabs(db)) < tol * (1.0 + std::max(a, b)))
            {
            converged = true;
            break;
            }

        // backtracking line search (keep positive and improve ll)
        double step = 1.0;
        bool accepted = false;
        for (int ls=0; ls<30; ls++)
            {
            double na = a + step * da;
            double nb = b + step * db;

            if (na <= 1e-12 || nb <= 1e-12 || !std::isfinite(na) || !std::isfinite(nb))
                {
                step *= 0.5;
                continue;
                }

            double nll = logLik(na, nb);
            if (std::isfinite(nll) && nll >= ll)
                {
                a = na;
                b = nb;
                ll = nll;
                accepted = true;
                break;
                }

            step *= 0.5;
            }

        if (accepted == false)
            break;
        }

    BetaFitResult r;
    r.alpha = a;
    r.beta  = b;
    r.iterations = it;
    r.converged  = converged;
    return r;
}

std::pair<double,double> Statistics::flatDirichletMarginalPriorBeta(size_t K) {

    if (K < 2)
        throw std::runtime_error("K must be >= 2");

    return std::make_pair(1.0, (double)(K - 1));
}

double Statistics::klBeta(double a, double b, double c, double d) {

    double psiAB = digamma(a + b);

    return (logBeta(c, d) - logBeta(a, b))
         + (a - c) * (digamma(a) - psiAB)
         + (b - d) * (digamma(b) - psiAB);
}

double Statistics::logBeta(double a, double b) {

    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

double Statistics::trigamma(double x) {

    if (x <= 0.0)
        throw std::runtime_error("trigamma: x must be > 0");

    double result = 0.0;

    // recurrence: psi1(x) = psi1(x+1) + 1/x^2
    while (x < 8.0)
        {
        result += 1.0 / (x * x);
        x += 1.0;
        }

    // asymptotic expansion
    double inv  = 1.0 / x;
    double inv2 = inv * inv;
    double inv3 = inv2 * inv;
    double inv5 = inv3 * inv2;
    double inv7 = inv5 * inv2;
    double inv9 = inv7 * inv2;

    result += inv
            + 0.5 * inv2
            + (1.0 / 6.0) * inv3
            - (1.0 / 30.0) * inv5
            + (1.0 / 42.0) * inv7
            - (1.0 / 30.0) * inv9;

    return result;
}
