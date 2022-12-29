mutable struct MultihistogramData
    observables::Vector
    parameter_values::Vector
    simulated_histograms::Vector
    observable_ranges::Vector
    marginal_energy_histograms::Vector
    free_energies::Vector
    tuple_iterators::Vector
end

function getrange(df, col)
    vec = getproperty(df, col)
    return minimum(vec):minimum(filter(t -> t != 0, abs.(diff(vec)))):maximum(vec)
end

function marginalvalue(df, col, val)
    return sum(df[getproperty(df, col) .== val, :].f)
end

function marginalize(df, col)
    vals = collect(getrange(df, col))
    freqs = [marginalvalue(df, col, val) for val in vals]
    return sort(Dict(vals .=> freqs))
end

function MultihistogramData(nparams::Integer, simulated_parameters::Vector, simulated_histograms::Vector)
    @assert allequal(names.(simulated_histograms))
    
    for (idx, histogram) in enumerate(simulated_histograms)
        @assert length(names(histogram)) == nparams+1 "Histogram $idx does not have correct number of parameters"
    end
    
    observables = Symbol.(names(first(simulated_histograms)))[1:end-1]
    
    if !(:U in observables)
        @error "Energy not found in observables. Please ensure energy is labelled with :U"
        throw(KeyError(:U))
    end
    
    observable_ranges = [
        Dict([C => getrange(histogram, C) for C in observables]...)
        for histogram in simulated_histograms]
    
    marginal_energy_histograms = [marginalize(histogram, :U) for histogram in simulated_histograms]
    tuple_iterators = Tables.namedtupleiterator.(simulated_histograms)
    return MultihistogramData(observables, simulated_parameters, simulated_histograms, observable_ranges, marginal_energy_histograms, [], tuple_iterators)
end