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
            num = f
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total += num/den
        end
    end
    return log(total)
end

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
            num = f * (E - Omin)
            den = sum([exp(-A-u[j]+(β - 1/Tvec[j])*E) for j in 1:J])
            total += num/den
        end
    end
    return exp(log(total) - F) + Omin
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

function interpolate_observable(parameter_value, observable, data::MultihistogramData; isbeta = false)
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
    Evec = data.marginal_energy_histograms
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
    Evec = data.marginal_energy_histograms
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