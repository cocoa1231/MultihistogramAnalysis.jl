module MultihistogramAnalysis

using DataFrames
using Tables
using ProgressMeter

export MultihistogramData

# Export a billion things because I CAN'T WRITE CODE
export getrange,
    marginalvalue,
    marginalize,
    free_energy_iteration,
    free_energy_iteration_logsum,
    calculate_free_energies!,
    interpolate_free_energy,
    interpolate_free_energy_logsum,
    interpolate_energy_logsum,
    interpolate_energy,
    interpolate_energy_logsum,
    interpolate_energy_second_moment,
    interpolate_energy_second_moment_logsum,
    interpolate_observable,
    interpolate_observable_logsum,
    interpolate_observable_abs,
    interpolate_observable_abs_logsum,
    interpolate_observable_second_moment,
    interpolate_observable_second_moment_logsum,
    interpolate_observable_second_moment_abs,
    interpolate_observable_second_moment_abs_logsum

include("MultihistogramDataStructures.jl")
include("MultihistogramInterpolation.jl")

@doc raw"""
One iteration of the multihistogram (wham) free energy equations with the logsum
trick. Inputs are

- `u`: One dimensional array-like type. Guess for the free energies.
- `data::MultihistogramData`: Dataset for Multihistogram interpolation.
"""
function free_energy_iteration_logsum(u, data::MultihistogramData)
    J = length(data.parameter_values)
    A = -1/2 * (maximum(u) + minimum(u))
    Evec = data.marginal_energy_histograms
    Tvec = data.parameter_values

    nterms = sum([length(H) for H in Evec])
    Fk = zeros(J)
    denom_arrays = [zeros(J) for _ in 1:Threads.nthreads()]
    terms = zeros(nterms)
    
    for k in 1:J
        offset = 0
        for i in 1:J
            iterthing = data.marginal_energy_histiters[i]
            Threads.@threads for idx in 1:length(iterthing)
                for j in 1:J
                    denom_arrays[Threads.threadid()][j] = -A-u[j]+(1/Tvec[k] - 1/Tvec[j])*first(iterthing[idx])
                end
                terms[offset + idx] = log(last(iterthing[idx])) - logsum!(denom_arrays[Threads.threadid()])
            end
            offset += length(iterthing)
        end
        Fk[k] += logsum!(terms)
    end
    
    return Fk
end

@doc raw"""
One iteration of the multihistogram (wham) free energy equations. Inputs are

- `u`: One dimensional array-like type. Guess for the free energies.
- `data::MultihistogramData`: Dataset for Multihistogram interpolation.
"""
function free_energy_iteration(u, data::MultihistogramData)
    J = length(data.parameter_values)
    A = -1/2 * (maximum(u) + minimum(u))
    Evec = data.marginal_energy_histograms
    Tvec = data.parameter_values
    
    Fk = zeros(J)
    for k in 1:J
        total_k = 0.
        for i in 1:J
            for (E, f) in pairs(Evec[i])
                num = f
                den = sum([exp(-A-u[j]+(1/Tvec[k] - 1/Tvec[j])E) for j in 1:J])
                total_k += num/den
            end
        end
        Fk[k] += log(total_k)
    end
    return Fk
end

@doc raw"""
Iterate the free energy multihistogram equation until the error is below a
threshold. Error is calculated as.

```math
\Delta^2 = \sum_{i} \left(\frac{F^(m)_i - F^(m-1)_i}{F^(m-1)_i}\right)^2
```

Inputs are

- `data::MultihistogramData`: Multihistogram dataset.
- `rtol` (default `1e-10`): Tolerance until the error is minimized. Value of
  ``\Delta``.
- `logsum` (default `false`): Use the logsum method. Recommended for systems
  with energies larger than 300.
"""
function calculate_free_energies!(data::MultihistogramData; rtol = 1e-10, logsum = false)
    P = ProgressThresh(rtol^2, "Minimizing: ")
    F = Fnew = zeros(length(data.parameter_values))

    f = free_energy_iteration
    if logsum
        f = free_energy_iteration_logsum
    end

    while true
        Fnew = f(F, data)
        err = sum( ((Fnew .- F) ./ F ) .^2 )
        ProgressMeter.update!(P, err)
        if err < rtol^2
            break
        end
        F = Fnew
    end
    data.free_energies = Fnew
    return Fnew
end

end # module MultihistogramAnalysis
