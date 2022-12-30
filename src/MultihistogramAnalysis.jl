module MultihistogramAnalysis

using DataFrames
using Tables
using ProgressMeter

export MultihistogramData

export getrange,
    marginalvalue,
    marginalize,
    free_energy_iteration,
    calculate_free_energies!,
    interpolate_free_energy,
    interpolate_energy,
    interpolate_energy_second_moment,
    interpolate_observable,
    interpolate_observable_second_moment,
    interpolate_observable_second_moment_abs

include("MultihistogramDataStructures.jl")

function free_energy_iteration(u, data::MultihistogramData)
    J = length(data.parameter_values)
    A = -1/2 * (maximum(u) + minimum(u))
    Evec = data.marginal_energy_histograms
    Tvec = data.parameter_values
    
    Fk = zeros(J)
    for k in 1:J
        total_k = 0.
        for i in 1:J
            for (E, f) in Evec[i]
                num = f
                den = sum([exp(-A-u[j]+(1/Tvec[k] - 1/Tvec[j])E) for j in 1:J])
                total_k += num/den
            end
        end
        Fk[k] += log(total_k)
    end
    return Fk
end

function calculate_free_energies!(data::MultihistogramData; rtol = 1e-10)
    P = ProgressThresh(rtol^2, "Minimizing: ")
    F = Fnew = zeros(length(data.parameter_values))
    while true
        Fnew = free_energy_iteration(F, data)
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

include("MultihistogramInterpolation.jl")

end
