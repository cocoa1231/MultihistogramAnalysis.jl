using MultihistogramAnalysis
using Test
using CSV
using DataFrames
using ProgressMeter

@testset "MultihistogramAnalysis.jl" begin
    # Open Datasets
    temperatures = Float64[]
    histograms = DataFrame[]
    dataset_dir = "./sample_dataset/"
    cd(dataset_dir)
    for dir in filter(isdir, readdir())
        push!(temperatures, parse(Float64, dir))
        push!(histograms, CSV.read(dir * "/histogram.csv", DataFrame))
    end

    # Create MHData object
    MHData = MultihistogramData(2, temperatures, histograms)
    calculate_free_energies!(MHData; logsum = true)

    N = 684
    Trange = minimum(temperatures):0.01:maximum(temperatures)
    βvec   = 1 ./ Trange
    Uvec   = zeros(length(Trange))
    Ufirstvec = zeros(length(Trange))
    Usqvec = zeros(length(Trange))
    Mvec   = zeros(length(Trange))
    Mfirstvec = zeros(length(Trange))
    Msqvec = zeros(length(Trange))

    @showprogress "Interpolating data" for (idx, T) in enumerate(Trange)
        U, Usq = interpolate_energy_second_moment_logsum(T, MHData, returnlinear = true)
        M, Msq = interpolate_observable_second_moment_abs_logsum(T, :M, MHData, returnlinear = true)
        
        Ufirstvec[idx] = interpolate_energy_logsum(T, MHData)
        Mfirstvec[idx] = interpolate_observable_abs_logsum(T, :M, MHData)
        Uvec[idx] = U
        Mvec[idx] = M
        Usqvec[idx] = Usq
        Msqvec[idx] = Msq
    end
    
    specheat = @. βvec / N * (Usqvec - Uvec^2)
    magneticsus = @. βvec * N * (Msqvec/N^2 - (Mvec/N)^2)

    specmax_idx = findfirst(==(maximum(specheat)), specheat)
    susmax_idx  = findfirst(==(maximum(magneticsus)), magneticsus)

    @test Trange[specmax_idx] == 1.66
    @test Trange[susmax_idx]  == 1.66
    @test all(Ufirstvec .== Uvec)
    @test all(Mfirstvec .== Mvec)
    
end
