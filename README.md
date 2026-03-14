# Quasi-Monte Carlo for Asian Option Pricing

A C++ implementation of **Randomized Quasi-Monte Carlo (RQMC)** using Sobol low-discrepancy sequences, benchmarked against standard Monte Carlo for pricing Asian arithmetic-average call options.

---

## Repository Structure

```
Quasi-Monte-Carlo/
├── sobol.cpp / sobol.hpp    # Sobol sequence generator using Joe-Kuo direction 
├── option_pricing.cpp       # Main logic behind pricing Asian Options using both methods
├── sobol_dump.cpp           # Writes Sobol points to a .tsv file for visualization
├── new-joe-kuo-6.21201      # Direction numbers for sobol sequence implementation
├── QMC_Asian_Options.ipynb  # Full analysis notebook: compiles C++, runs benchmark, plots
├── dimensional_comparison/  # Saved plots for d=4 and d=12 time steps
└── README.md
```

---

## Quickstart

**Requirements:** `g++` with C++20, Python 3+, and `pandas matplotlib scipy`.

```bash
git clone https://github.com/devDivij/Quasi-Monte-Carlo.git
cd Quasi-Monte-Carlo
```
```bash
pip install pandas matplotlib scipy
```
The recommended way to run everything is through the Jupyter notebook, which compiles the C++ binaries, runs the benchmark, and generates all plots in sequence. Open `QMC_Asian_Options.ipynb` in VS Code (with the Jupyter extension) or run:
```bash
pip install notebook ipykernel jupyter
jupyter notebook QMC_Asian_Options.ipynb
```

All simulation parameters — spot price, strike, volatility, time steps — are defined as Python variables at the top of the notebook. Changing them and re-running passes the new values as command-line arguments to the binaries, so there is no need to touch any C++ file to explore different scenarios.

To compile and run manually without the notebook: (Simulation constants may be changed inside `option_pricing.cpp` as well)

```bash
g++ -O2 -std=c++20 option_pricing.cpp sobol.cpp -o option_pricing
./option_pricing 100.0 100.0 0.05 1.0 0.20 30
```
OR simulation constants may be changed inside `option_pricing.cpp` as well, with running just 
```bash
g++ -O2 -std=c++20 option_pricing.cpp sobol.cpp -o option_pricing
./option_pricing
```
---
## Results

The core claim of QMC theory is that a Sobol sequence covers the integration domain more uniformly than pseudo-random samples, and this translates directly into a faster convergence rate — roughly $O((\log N)^d / N)$ vs standard MC's $O(1/\sqrt{N})$. The plots below make this concrete.

**$n\sigma^2$ metric** has been used to compare the two simulation methods. For standard MC, the product of sample size and estimator variance equals the process variance $\sigma^2$ — a constant completely independent of $N$. This gives a flat horizontal baseline. For a good QMC method, this quantity should *decrease* as $N$ grows, because the estimator variance falls faster than $1/N$. The gap between the two lines is the variance reduction factor.

![Convergence](dimensional_comparison/convergence(time_steps=4).png)

The variance reduction factor (VRF) tells you how many times fewer samples RQMC needs to match MC's accuracy. A VRF of 30 at $N = 1024$ means RQMC with 1,024 points is as accurate as MC with roughly 30,000 points.

![Variance Reduction](dimensional_comparison/variance_reduction(time_steps=4).png)

---
## Method

**The option.** An Asian arithmetic-average call option is one of the canonical test cases for high-dimensional numerical integration in finance, because its payoff depends on the average of the asset price across all $d$ monitoring dates, not just the terminal value. Pricing it requires integrating over a $d$-dimensional Gaussian — exactly the kind of problem where QMC's uniform coverage advantage shows up most clearly.

The asset price follows geometric Brownian motion, discretised with the log-Euler scheme:

$$S_{i+1} = S_i \exp\!\left[\left(r - \tfrac{1}{2}\sigma^2\right)\Delta t + \sigma\sqrt{\Delta t}\;\Phi^{-1}(U_i)\right]$$

where $U_i \in (0,1)$ is the $i$-th uniform sample and $\Phi^{-1}$ is the inverse normal CDF, computed via Acklam's rational approximation with error below $1.15 \times 10^{-9}$. The discounted payoff is:

$$V = e^{-rT} \max\!\left(\frac{1}{d}\sum_{i=1}^d S_i - K,\; 0\right)$$

**The comparison.** Standard MC draws each $U_i$ from a Mersenne Twister pseudo-random generator. RQMC replaces this with a Sobol low-discrepancy sequence, randomized by a **random digital shift** (Cranley-Patterson rotation): a single uniform random vector $\Delta \in [0,1)^d$ is added modulo 1 to every point in the net. This preserves the low-discrepancy structure of each shifted copy while making independent replicates statistically independent — which is what allows variance estimation from a deterministic sequence.

**The variance estimator.** $R = 30$ independent random shifts are applied to the same $N$-point Sobol net, producing $R$ independent price estimates. The variance of those replicate means is the RQMC estimator variance. Multiplied by $N$, this gives Glasserman's $n\sigma^2$ metric. Crucially, each replicate uses the *full* $N$-point net — the replicates exist only to measure uncertainty, not to improve the price estimate itself.

**Why $n\sigma^2$ and not just the standard error?** Because different QMC methods have natural sample sizes — Sobol works best at powers of 2, Faure sequences at powers of their base, and so on. Reporting $n\sigma^2$ lets you compare methods at their respective natural sizes while producing a number on a common scale. MC's flat $\sigma^2$ baseline then serves as a universal reference: any QMC method whose $n\sigma^2$ is below the MC baseline is provably more efficient per sample, and any method whose $n\sigma^2$ is *decreasing* in $N$ is converging faster than MC's $O(1/\sqrt{N})$ rate.

---
 
## Sobol Sequence Implementation
 
The sequence generator in `sobol.cpp` is implemented from scratch following the algorithm described in Glasserman Chapter 5, using the Joe-Kuo direction numbers from `new-joe-kuo-6.21201` which support up to 21,201 dimensions.
 
**Direction numbers and the generating matrix.** Each dimension $i$ of a Sobol sequence is defined by a $r \times r$ binary generating matrix $V^{(i)}$, where $r = \lceil \log_2 N \rceil$. The columns of this matrix are the *direction numbers* — carefully chosen binary fractions that encode how the sequence fills each bit position. For dimension 1 the matrix is the identity (giving the Van der Corput sequence in base 2); for all other dimensions it is constructed from the primitive polynomial and initial direction numbers in the Joe-Kuo file via the three-term recurrence:
 
$$m_j = 2^q m_{j-q} \oplus \left(2^{q-1} a_1 m_{j-1} \oplus \cdots \oplus 2^1 a_{q-1} m_{j-q+1}\right) \oplus m_{j-q}$$
 
where $a_1, \ldots, a_{q-1}$ are the coefficients of the primitive polynomial of degree $q$, and $\oplus$ denotes bitwise XOR. This is what `Sobolmat()` computes — it fills the lower-triangular part of $V^{(i)}$ by applying this recurrence to extend the initial $m_1, \ldots, m_q$ values from the direction numbers file.
 
**Gray code ordering.** Rather than iterating $k = 0, 1, 2, \ldots$ directly, the implementation uses **Gray code ordering** — $g(k) = k \oplus (k \gg 1)$ — which has the property that consecutive integers differ in exactly one bit. This means each new Sobol point differs from the previous one by XORing in a single column of the generating matrix, rather than recomputing the full matrix-vector product from scratch. The `rightmostbit()` function identifies which column to XOR at each step by finding the position of the rightmost zero bit of $k-1$, which is precisely where $g(k)$ and $g(k-1)$ differ.
 
**The point computation.** The $k$-th Sobol point in dimension $i$ is produced by accumulating the running Gray-code-ordered XOR vector $y^{(i)}$ and converting the resulting $r$-bit integer to a real number in $(0,1)$ by weighting each bit $y_j$ as $y_j \cdot 2^{-(j+1)}$. This is equivalent to computing $V^{(i)} \cdot g(k) \pmod{2}$ in $\mathbb{F}_2^r$ and then interpreting the result as a binary fraction, but the incremental Gray code update makes it $O(r)$ per point instead of $O(r^2)$.
 
**Arbitrary starting point.** The `Sobolpts(n0, npts, d, S)` interface supports starting from an arbitrary index `n0` rather than always beginning at 0. This is useful for the convergence benchmark, where Glasserman recommends starting each $N$-point net at the $N$-th point of the sequence (so successive runs use non-overlapping, well-spaced sub-sequences). The implementation initialises $y$ by dotting the full generating matrix against the Gray code of `n0 - 1`, then resumes the incremental updates from there.
 
---

## References

Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*. Springer. Chapter 5: Quasi-Monte Carlo.