function logsum(logterms)
    lm = maximum(logterms)
    logterms_removed = filter(!=(lm), logterms)
    return lm + log1p(sum(exp.(logterms_removed .- lm)))
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
        for (E, f) in Evec[i]
            if f == 0
                continue
            end
            num = f
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total += num/den
        end
    end
    return log(total)
end

function interpolate_free_energy_logsum(parameter_value, data::MultihistogramData; isbeta = false)
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
    
    logterms = Float64[]
    for i in 1:J
        for (E, f) in Evec[i]
            if f == 0
                continue
            end
            num = f
            den = [-A-u[j]+(β - 1/Tvec[j])*E for j in 1:J]
            push!(logterms, log(num) - logsum(den))
        end
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
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    for i in 1:J
        for (E, f) in Evec[i]
            if f == 0
                continue
            end

            num = f * (E - Omin)
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total += num/den
        end
    end
    return exp(log(total) - F) + Omin
end

function interpolate_energy_logsum(parameter_value, data::MultihistogramData; isbeta = false)
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms = Float64[]
    for i in 1:J
        for (E, f) in Evec[i]
            if f == 0
                continue
            end

            num = f * (E - Omin)
            den = [-A-u[j]+(β - 1/Tvec[j])*E for j in 1:J]
            push!(logterms, log(num) - logsum(den))
        end
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
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    for i in 1:J
        for (E, f) in Evec[i]
            if f == 0
                continue
            end

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

function interpolate_energy_second_moment_logsum(parameter_value, data::MultihistogramData;
        isbeta = false, returnlinear = false)    
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    Evec = data.marginal_energy_histograms
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[:U]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms_linear = Float64[]
    logterms_square = Float64[]
    for i in 1:J
        for (E, f) in Evec[i]
            if f == 0
                continue
            end
            num_li = f * (E - Omin)
            num_sq = num_li * (E - Omin)
            
            den = logsum([-A-u[j]+(β - 1/Tvec[j])*E for j in 1:J])
            push!(logterms_linear, log(num_li) - den)
            push!(logterms_square, log(num_sq) - den)
        end
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
    F    = interpolate_free_energy(parameter_value, data)
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

function interpolate_observable_logsum(parameter_value, observable, data::MultihistogramData; isbeta = false)
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[observable]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms = Float64[]
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num = state.f * (getproperty(state, observable) - Omin)
            den = [ -A-u[j]+(β - 1/Tvec[j]) * state.U for j in 1:J]
            push!(logterms, log(num) - logsum(den))
        end
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
    F    = interpolate_free_energy(parameter_value, data)
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

function interpolate_observable_abs_logsum(parameter_value, observable, data::MultihistogramData; isbeta = false)
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))

    logterms = Float64[]
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num = state.f * abs(getproperty(state, observable))
            den = [ -A-u[j]+(β - 1/Tvec[j]) * state.U for j in 1:J]
            push!(logterms, log(num) - logterms(den))
        end
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
    F    = interpolate_free_energy(parameter_value, data)
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

function interpolate_observable_second_moment_logsum(parameter_value, observable, data::MultihistogramData;
        isbeta = false, returnlinear = false)
    
    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end
    
    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))
    Oranges = (x -> x[observable]).(data.observable_ranges)
    Omin = minimum(minimum.(Oranges))
    
    logterms_linear = Float64[]
    logterms_square = Float64[]
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num_li = state.f * (getproperty(state, observable) - Omin)
            num_sq = num_li * (getproperty(state, observable) - Omin)
            den = logsum([-A-u[j]+(β - 1/Tvec[j] * state.U) for j in 1:J])
            
            push!(logterms_linear, log(num_li) - den)
            push!(logterms_square, log(num_sq) - den)    
        end
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
    F    = interpolate_free_energy(parameter_value, data)
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

function interpolate_observable_second_moment_abs_logsum(parameter_value, observable, data::MultihistogramData;
        isbeta = false, returnlinear = false)
    total_linear = 0.
    total_square = 0.

    if !isbeta
        β = 1/parameter_value
    else
        β = parameter_value
    end

    Tvec = data.parameter_values
    F    = interpolate_free_energy(parameter_value, data)
    u    = data.free_energies
    J    = length(Tvec)
    A    = -1/2 * (maximum(u) + minimum(u))

    logterms_linear = Float64[]
    logterms_square = Float64[]
    for i in 1:J
        for state in data.tuple_iterators[i]
            if state.f == 0
                continue
            end
            num_li = state.f * abs(getproperty(state, observable))
            num_sq = num_li * abs(getproperty(state, observable))
            den = logsum([ -A-u[j]+(β - 1/Tvec[j]) * state.U for j in 1:J])
            
            push!(logterms_linear, log(num_li) - den)
            push!(logterms_square, log(num_sq) - den)
        end
    end

    Oli = exp(logsum(logterms_linear) - F)
    Osq = exp(logsum(logterms_square) - F)

    if returnlinear
        return Oli, Osq
    else
        return Osq
    end
end