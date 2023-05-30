```@meta
CurrentModule = MultihistogramAnalysis
```

# MultihistogramAnalysis

Documentation for
[MultihistogramAnalysis](https://github.com/cocoa1231/MultihistogramAnalysis.jl).

!!! warning
    At the moment the library will only work with systems that have a discreet
    energy spectra. Generation of histograms neccesarily means that we bin data,
    and thus throw away data points if the energy spectra is not discreet. Use
    with caution on systems with continuous energy spectra.

The multihistogram method solves for the density of states by solving for the
partition functions at the temperatures simulated at. The main equation that is
solved iteratively is.

```math
Z(\beta_k) = \sum_{i}\sum_E \frac{N_i(E)}{\sum_j n_j Z^{-1}(\beta_j) e^{(\beta_k
- \beta_j)E}}
```

Where ``i,j`` indexes the simulation number, ``E`` are the energies sampled in
the i-th simulation, ``N_i`` is the energy histogram (i.e, number of times an
energy occurs in the simulation), and `n_i` are the number of samples in the
i-th histogram. For this library, ``n_i = 0`` for all ``i``. Since the partition
functions depend exponentially on the energy, it is better to use the logarithms
(``F(\beta_k) = F_k = \log(Z(\beta_k)``)
of the partition functions to avoid float overflow. Thus the equation becomes.

```math
F_k =\log\left( \sum_i \sum_E \frac{N_i(E)}{\sum_j \log(-F_k + (\beta_k - \beta
)E} \right)
```

To avoid overflow further, we normalize all partition functions by the following
factor.

```math
A = \frac{1}{\sqrt{Z_{max}Z_{min}}} \implies \log(A) = -\frac{1}{2} (F_{max} +
F_{min})
```

To interpolate the free energy we use a natural extension of the above equation.

```math
F(\beta) = \log\left( \sum_i \sum_E \frac{N_i(E)}{\sum_j \log(-F_k + (\beta -
\beta)E} \right)
```

And the calculate the thermal average we use the following.

```math
\langle Q \rangle = \frac{1}{Z(\beta)} \sum_{i, (Q, E)} Q \frac{N_i(Q,
E)}{\sum_j A Z^{-1}_j e^{(\beta - \beta_j)E}}
```

## Logsum Methods

The logsum trick allows us to calculate the logarithm of a sum of a series of
logarithms that we want to exponentiate. Specifically, if ``l_i`` are
logarithms, and we want to calculate ``\log(\sum_i e^{l_i})`` without ever
exponentiating the logarithms, this trick lets us do that. If ``l_0`` is the
largest logarithm, then one can calculate the sum as follows.

```math
\log(\sum_{i=0} \exp(l_i)) = \log(e^{l_0} \left[1 + \sum_{i=1} e^{l_i -
l_0}\right]) = l_0 + \log(1 +  \sum_{i=1} e^{l_i - l_0})
```

This has the advantage of avoiding numerical overflow and calculate observables
for very large systems. A function name that ends with `_logsum` uses this
method. Note that this has also been implemented for `ComplexF64` vectors.

## Interpolatioon Functions

After importing the library, you can write `interpolate_` and try to
autocomplete with `TAB` to see all the interpolation methods.
Counterintuatively, the `MultihistogramData` structs go as the last argument
rather than the first. I did not have a lot of time to write this library
cleanly, will be fixed in the future.
