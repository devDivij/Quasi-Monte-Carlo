#include "sobol.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
#include <numeric>
#include <iomanip>

// Model parameters

double S0 = 100.0;
double K = 100.0;
double r_f = 0.05;
double T = 1.0;
double sigma = 0.20;
int time_steps = 12; // Dimension of the problem (number of time steps in the Asian option)

// Acklam's approximation to Phi^{-1}(U)

double inverseNormalCDF(double p)
{
    if (p <= 0.0 || p >= 1.0)
        return 0.0;

    static const double a[] = {
        -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {
        -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01};
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {
        7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00};

    const double p_low = 0.02425, p_high = 1.0 - p_low;
    double z;
    if (p < p_low)
    {
        double q = std::sqrt(-2.0 * std::log(p));
        z = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    else if (p > p_high)
    {
        double q = std::sqrt(-2.0 * std::log(1.0 - p));
        z = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    else
    {
        double q = p - 0.5, r = q * q;
        z = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    return z;
}

// Asian arithmetic-average call payoff, discounted.
// U is a d-dimensional point in (0,1)^d.

double asianCallPayoff(const std::vector<double> &U)
{
    const double dt = T / time_steps;
    const double drift = (r_f - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);

    double S = S0, sumS = 0.0;
    for (int i = 0; i < time_steps; ++i)
    {
        S *= std::exp(drift + vol * inverseNormalCDF(U[i]));
        sumS += S;
    }
    return std::exp(-r_f * T) * std::max(sumS / time_steps - K, 0.0);
}

// Standard Monte Carlo

std::pair<double, double> standardMC(int N, std::mt19937 &rng)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double sum = 0.0, sumSq = 0.0;
    std::vector<double> U(time_steps);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < time_steps; ++j)
            U[j] = dist(rng);
        double p = asianCallPayoff(U);
        sum += p;
        sumSq += p * p;
    }
    double mean = sum / N;
    double estimator_var = sumSq / N - mean * mean; // sigma^2, independent of N

    // Glasserman metric: N * Var[estimator] = N * (sigma^2/N) = sigma^2
    return {mean, estimator_var};
}

// Randomized QMC with digital shift [random shift of point with a d-dimensional (here d = time_steps) U vector modulo 1]
//
// Each replicate:
// 1. Take the fixed N-point Sobol net
// 2. Add an independent uniform random shift (mod 1) to every point
// 3. Compute the sample mean over the shifted net
//
// Variance is estimated across R replicate means.
// Metric used: N * Var[replicate means]
//
// Key: each replicate uses the FULL N-point net (not a sub-batch).
// The R replicates exist only to estimate variance — not to improve accuracy.

std::pair<double, double> randomizedQMC(int N, int R, std::mt19937 &rng, const std::vector<Sobolparams> &params)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<std::vector<double>> base = Sobolpts(1, N, time_steps, params);
    std::vector<double> replicateMeans(R);

    for (int rep = 0; rep < R; ++rep)
    {
        std::vector<double> shift(time_steps);
        for (int j = 0; j < time_steps; ++j)
            shift[j] = dist(rng);

        double sum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            std::vector<double> U(time_steps);
            for (int j = 0; j < time_steps; ++j)
                U[j] = std::fmod(base[i][j] + shift[j], 1.0);
            sum += asianCallPayoff(U);
        }
        replicateMeans[rep] = sum / N;
    }

    double mean = 0.0;
    for (double m : replicateMeans)
        mean += m;
    mean /= R;

    double var = 0.0;
    for (double m : replicateMeans)
        var += (m - mean) * (m - mean);
    var /= (R - 1);

    return {mean, N * var};
}

int main(int argc, char *argv[])
{
    if (argc >= 7)
    {
        S0 = std::stod(argv[1]);
        K = std::stod(argv[2]);
        r_f = std::stod(argv[3]);
        T = std::stod(argv[4]);
        sigma = std::stod(argv[5]);
        time_steps = std::stoi(argv[6]);
    }

    std::mt19937 rng(42);
    const int R = 30;

    std::vector<Sobolparams> sobolParams = extract_Sobolparams("new-joe-kuo-6.21201", time_steps);

    // Sobol sequences have good properties for powers of 2, but we can also test some intermediate values to see the trend.
    std::vector<int> Nvals = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288};

    std::cout << "N\tmc_price\tmc_nvar\trqmc_price\trqmc_nvar\n";

    for (int N : Nvals)
    {
        auto [mcPrice, mcNVar] = standardMC(N, rng);
        auto [rqmcPrice, rqmcNVar] = randomizedQMC(N, R, rng, sobolParams);

        std::cout << std::fixed << std::setprecision(6) << N << "\t" << mcPrice << "\t" << mcNVar << "\t" << rqmcPrice << "\t" << rqmcNVar << std::endl;
    }

    return 0;
}