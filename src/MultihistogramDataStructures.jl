using Dictionaries
using DataStructures

@doc raw"""
Primary data structure of the library. Holds all the information needed to run
all of the functions provided. The constructor provided for this will precompute
all the required datasets to ensure fastest interpolation. The fields are

- `observables::Vector{Symbol}`: A vector of symbols comprising of all the
  observables of the system. For instance, for an Ising system, this will be the
  internal energy and magnetization. It is mandatory that the column for internal
  energy be labelled with `:U`. Magnetization can be labelled whatever, so for
  instance a valid observables vector is `Symbol[:U, :M]`.
- `parameter_values::Vector{pType}`: List of parameters at which the system was
  simulated. For instnace, this can be the temperatures at which MCMC simulations
  were done.
- `simulated_histograms::Vector{DataFrame}`: Flattened vector of all simulated
  histograms. The DataFrame must have a column for internal energy labelled `:U`
  and the frequency of occurance labelled `:f`. Other columns can be named
  whatever.
- `observable_ranges::Vector{orType}`: Vector of the ranges spanned by the
  observables.
- `marginal_energy_histograms::Vector{ehType}`: Dictionary of marginal energy
  histograms, mapping energy to frequency of occurance.
- `marginal_energy_histiters::Vector{ehiType}`: Array of collected pairs of
  `marginal_energy_histograms`. Collected they take less memory.
- `nEbins::Int64`: Number of energy bins sampled in total.
- `nbins::Int64`: Number of histogram bins in total.
- `free_energies::Vector{pType}`: Current best guess of free energies.
- `tuple_iterators::{ntiType}`: Named tuple iterators from the histograms for
  fast iterations.
"""
mutable struct MultihistogramData{pType, orType, ehType, ehiType, ntiType}
    observables::Vector{Symbol}
    parameter_values::Vector{pType}
    simulated_histograms::Vector{DataFrame}
    observable_ranges::Vector{orType}
    marginal_energy_histograms::Vector{ehType}
    marginal_energy_histiters::Vector{ehiType}
    nEbins::Int64
    nbins::Int64
    free_energies::Vector{pType}
    tuple_iterators::Vector{ntiType}
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

function MultihistogramData(nparams::Int64, simulated_parameters::Vector{pType}, simulated_histograms::Vector{hType}; atol = 1e-3) where {pType, hType}
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
    marginal_energy_histiters = [collect(pairs(H)) for H in marginal_energy_histograms]
    tuple_iterators = collect.(Tables.namedtupleiterator.(simulated_histograms))
    nEbins = sum([length(H) for H in marginal_energy_histograms])
    nbins = sum([length(H) for H in tuple_iterators])
    return MultihistogramData(observables, simulated_parameters, simulated_histograms, observable_ranges,
        marginal_energy_histograms, marginal_energy_histiters, nEbins, nbins, pType[], tuple_iterators)
end
