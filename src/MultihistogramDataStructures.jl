using Dictionaries
using DataStructures

mutable struct MultihistogramData
    observables::Vector
    parameter_values::Vector
    simulated_histograms::Vector
    observable_ranges::Vector
    marginal_energy_histograms::Vector
    nEbins::Int64
    nbins::Int64
    free_energies::Vector
    tuple_iterators::Vector
end

function getrange(df, col)
    vec = getproperty(df, col)
    if allequal(vec)
        dv = abs(vec[1])
    else
        dv = minimum(filter(t -> t != 0, abs.(diff(vec))))
    end
    return minimum(vec):dv:maximum(vec)
end

function marginalvalue(df, col, val; atol)
    return sum(df[(x -> isapprox(x, val, atol=atol)).(getproperty(df, col)), :].f)
end

function marginalize(df, col; atol = 1e-3)
    vals = collect(getrange(df, col))
    freqs = [marginalvalue(df, col, val; atol = atol) for val in vals]
    histogram = filter(x -> last(x) != 0, Dict(vals .=> freqs)) |> SortedDict
    return Dictionary(histogram)
end

function MultihistogramData(nparams::Integer, simulated_parameters::Vector, simulated_histograms::Vector; atol = 1e-3)
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

    @info "Generating marginal histograms and tuple iterators"    
    marginal_energy_histograms = [marginalize(histogram, :U; atol = atol) for histogram in simulated_histograms]
    tuple_iterators = collect.(Tables.namedtupleiterator.(simulated_histograms))
    nEbins = sum([length(H) for H in marginal_energy_histograms])
    nbins = sum([length(H) for H in tuple_iterators])
    return MultihistogramData(observables, simulated_parameters, simulated_histograms, observable_ranges, marginal_energy_histograms, nEbins, nbins, [], tuple_iterators)
end