function logsum(logterms)
    lm = maximum(logterms)
    logterms_removed = filter(!=(lm), logterms)
    return lm + log1p(sum(exp.(logterms_removed .- lm)))
end

function logsum(logterms::Vector{ComplexF64})
    lr = real.(logterms)
    maxidx = findfirst(==(maximum(lr)), lr)
    maxterm = logterms[maxidx]
    logterms_filtered = filter(x -> x != maxterm, logterms)
    return maxterm + log1p(sum(exp.(logterms_filtered .- maxterm)))
end

### Free Energy

function interpolate_free_energy(parameter_value, data::MultihistogramData; isbeta = false)
    total = 0.
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    
    for i in 1:J
        for (E, f) in pairs(Evec[i])
            num = f
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total += num/den
        end
    end
    return log(total)
end

function interpolate_free_energy_logsum(parameter_value::T, data::MultihistogramData; isbeta = false) where T
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    nterms = sum([length(H) for H in Evec])
    
    logterms = zeros(T, nterms)
    offset = 0
    for i in 1:J
        iterthing = collect(pairs(Evec[i]))

        Threads.@threads for idx in 1:length(iterthing)
            E = first(iterthing[idx])
            f = last(iterthing[idx])
            logterms[offset + idx] = log(f) - logsum([
                -A-u[j]+(β - 1/Tvec[j])*E for j in 1:J
            ])
        end

        offset += length(iterthing)
    end

    return logsum(logterms)
end


### Internal Energy

function interpolate_energy(parameter_value, data::MultihistogramData; isbeta = false)
    total = 0.
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    F    = interpolate_free_energy(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    for i in 1:J
        for (E, f) in pairs(Evec[i])

            num = f * (E - Omin)
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total += num/den
        end
    end
    return exp(log(total) - F) + Omin
end

function interpolate_energy_logsum(parameter_value::T, data::MultihistogramData; isbeta = false) where T
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    F    = interpolate_free_energy_logsum(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
        
    logterms = zeros(T, data.nEbins)
    offset = 0
    for i in 1:J
        iterthing = collect(pairs(Evec[i]))
        Threads.@threads for idx in 1:length(iterthing)
            E = first(iterthing[idx])
            f = last(iterthing[idx])
            num = f * (E - Omin)
            den = [-A-u[j]+(β - 1/Tvec[j])*E for j in 1:J]
            logterms[offset + idx] = log(num) - logsum(den)
        end
        offset += length(iterthing)
    end
    return exp(logsum(logterms) - F) + Omin
end

function interpolate_energy_second_moment(parameter_value, data::MultihistogramData;
        isbeta = false, returnlinear = false)
    total_linear = 0.
    total_square = 0.
    
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    F    = interpolate_free_energy(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    for i in 1:J
        for (E, f) in pairs(Evec[i])

            num_li = f * (E - Omin)
            num_sq = num_li * (E - Omin)
            
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total_linear += num_li/den
            total_square += num_sq/den
        end
    end
    
    Oli = exp(log(total_linear) - F) + Omin
    Osq = exp(log(total_square) - F) + 2 * Oli * Omin - Omin^2
    
    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
end

function interpolate_energy_second_moment_logsum(parameter_value::T, data::MultihistogramData;
        isbeta = false, returnlinear = false) where T
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    F    = interpolate_free_energy_logsum(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms_linear = zeros(T, data.nEbins)
    logterms_square = zeros(T, data.nEbins)
    offset = 0
    for i in 1:J
        iterthing = collect(pairs(Evec[i]))
        for idx in 1:length(iterthing)
            E = first(iterthing[idx])
            f = last(iterthing[idx])
            num_li = f * (E - Omin)
            num_sq = num_li * (E - Omin)
            
            den = logsum([-A-u[j]+(β - 1/Tvec[j])*E for j in 1:J])
            logterms_linear[offset + idx] = log(num_li) - den
            logterms_square[offset + idx] = log(num_sq) - den
        end
        offset += length(iterthing)
    end
    
    Oli = exp(logsum(logterms_linear) - F) + Omin
    Osq = exp(logsum(logterms_square) - F) + 2 * Oli * Omin - Omin^2
    
    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
end

### Arbitrary Observable ⟨M⟩

function interpolate_observable(parameter_value, observable, data::MultihistogramData; isbeta = false)
    total = 0.
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[observable]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num = state.f * (getproperty(state, observable) - Omin)
            den = sum([ exp(-A-u[j]+(β - 1/Tvec[j]) * state.U) for j in 1:J])
            total += num/den
        end
    end
    
    return exp(log(total) - F) + Omin
end


function interpolate_observable_logsum(parameter_value::T, observable, data::MultihistogramData; isbeta = false) where T
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy_logsum(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[observable]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms = zeros(T, data.nbins)
    offset = 0
    for i in 1:J
        t = data.tuple_iterators[i]
        Threads.@threads for idx in 1:length(t)
            state = t[idx]
            num = state.f * (getproperty(state, observable) - Omin)
            den = [ -A-u[j]+(β - 1/Tvec[j]) * state.U for j in 1:J]
            logterms[offset + idx] = log(num) - logsum(den)
        end
        offset += length(t)
    end
    
    return exp(logsum(logterms) - F) + Omin
end

function interpolate_observable_abs(parameter_value, observable, data::MultihistogramData; isbeta = false)
    total = 0.
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))

    
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num = state.f * abs(getproperty(state, observable))
            den = sum([ exp(-A-u[j]+(β - 1/Tvec[j]) * state.U) for j in 1:J])
            total += num/den
        end
    end
    
    return exp(log(total) - F)
end

function interpolate_observable_abs_logsum(parameter_value::T, observable, data::MultihistogramData; isbeta = false) where T
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy_logsum(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))

    logterms = zeros(T, data.nbins)
    offset = 0
    for i in 1:J
        t = data.tuple_iterators[i]
        for idx in 1:length(t)
            state = t[idx]
            num = state.f * abs(getproperty(state, observable))
            den = [ -A-u[j]+(β - 1/Tvec[j]) * state.U for j in 1:J]
            logterms[offset + idx] = log(num) - logsum(den)
        end
        offset += length(t)
    end
    
    return exp(logsum(logterms) - F)
end

function interpolate_observable_second_moment(parameter_value, observable, data::MultihistogramData;
        isbeta = false, returnlinear = false)
    total_linear = 0.
    total_square = 0.
    
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[observable]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num_li = state.f * (getproperty(state, observable) - Omin)
            num_sq = num_li * (getproperty(state, observable) - Omin)
            den = sum([ exp(-A-u[j]+(β - 1/Tvec[j]) * state.U) for j in 1:J])
            
            total_linear += num_li/den
            total_square += num_sq/den
        end
    end
    
    Oli = exp(log(total_linear) - F) + Omin
    Osq = exp(log(total_square) - F) + 2 * Oli * Omin - Omin^2
    
    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
    
end

function interpolate_observable_second_moment_logsum(parameter_value::T, observable, data::MultihistogramData;
        isbeta = false, returnlinear = false) where T
    
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy_logsum(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[observable]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms_linear = zeros(T, data.nbins)
    logterms_square = zeros(T, data.nbins)
    offset = 0
    for i in 1:J
        t = data.tuple_iterators[i]
        for idx in 1:length(t)
            state = t[idx]
            num_li = state.f * (getproperty(state, observable) - Omin)
            num_sq = num_li * (getproperty(state, observable) - Omin)
            den = logsum([-A-u[j]+(β - 1/Tvec[j] * state.U) for j in 1:J])
            
            logterms_linear[offset + idx] = log(num_li) - den
            logterms_square[offset + idx] = log(num_sq) - den
        end
        offset += length(t)
    end
    
    Oli = exp(logsum(logterms_linear) - F) + Omin
    Osq = exp(logsum(logterms_square) - F) + 2 * Oli * Omin - Omin^2
    
    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
    
end


function interpolate_observable_second_moment_abs(parameter_value, observable, data::MultihistogramData;
        isbeta = false, returnlinear = false)
    total_linear = 0.
    total_square = 0.

    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end

    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))

    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num_li = state.f * abs(getproperty(state, observable))
            num_sq = num_li * abs(getproperty(state, observable))
            den = sum([ exp(-A-u[j]+(β - 1/Tvec[j]) * state.U) for j in 1:J])
            
            total_linear += num_li/den
            total_square += num_sq/den
        end
    end

    Oli = exp(log(total_linear) - F)
    Osq = exp(log(total_square) - F)

    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
end

function interpolate_observable_second_moment_abs_logsum(parameter_value::T, observable, data::MultihistogramData;
        isbeta = false, returnlinear = false) where T
    total_linear = 0.
    total_square = 0.

    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end

    Tvec = data.parameter_values
    F    = interpolate_free_energy_logsum(parameter_value, data, isbeta = isbeta)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))

    logterms_linear = zeros(T, data.nbins)
    logterms_square = zeros(T, data.nbins)
    offset = 0
    for i in 1:J
        t = data.tuple_iterators[i]
        for idx in 1:length(t)
            state = t[idx]
            num_li = state.f * abs(getproperty(state, observable))
            num_sq = num_li * abs(getproperty(state, observable))
            den = logsum([ -A-u[j]+(β - 1/Tvec[j]) * state.U for j in 1:J])
            
            logterms_linear[offset + idx] = log(num_li) - den
            logterms_square[offset + idx] = log(num_sq) - den
        end
        offset += length(t)
    end

    Oli = exp(logsum(logterms_linear) - F)
    Osq = exp(logsum(logterms_square) - F)

    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
end