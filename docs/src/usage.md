# Usage

Using this library is fairly simple. The only things to keep in mind are

1. Your input datasets should be *flattened* histograms in the form of a
   DataFrame.
2. Each DataFrame should have a `:U` column for the energy of the system and
   `:f` column for the corresponding frequency of that state. A sample dataset
   will be uploaded once a file hosting service is settled on.
3. **This library is not thread safe**. While multithreading has been
   implemented for the `logsum` methods, it may not neccesarily be faster, and
   there are diminishing gains as the number of threads increase.

An example usage of this library can be seen on the following dataset. Say we
simulate an Ising spin system for a number of different temperatures and output
the histograms to their CSVs. Our directory structure can look something like

```base
sample_dataset
├── 1.0
│   ├── histogram.csv
│   └── metadata.ini
├── 1.1
│   ├── histogram.csv
│   └── metadata.ini
├── 1.2
│   ├── histogram.csv
│   └── metadata.ini
├── 1.3
│   ├── histogram.csv
│   └── metadata.ini
├── 1.4
│   ├── histogram.csv
│   └── metadata.ini
```

and so on. The directory names are the temperatures and the histogram is always
stored in `histogram.csv`. Additionally some metadata about the simulation can
also be provided, however it is not needed in using this library. We can import
our dataset and create a vector of temperatures and histograms.

```julia
julia> using DataFrames, CSV, MultihistogramAnalysis

julia> datadir = "./sample_dataset/"
"./sample_dataset/"

julia> temperatures = Float64[]; histograms = DataFrame[]
DataFrame[]

julia> for simdir in filter(x -> isdir(joinpath(datadir, x)), readdir(datadir))
           try
               push!(temperatures, parse(Float64, simdir))
               push!(histograms, CSV.read(joinpath(datadir, simdir, "histogram.csv"), DataFrame))
           catch e
               if e isa ArgumentError
                   continue
               else
                   rethrow(e)
               end
           end
       end
```

Lets take a look at one of our histograms

```julia
julia> first(histograms[1], 5)
5×3 DataFrame
 Row │ M        U        f     
     │ Float64  Float64  Int64 
─────┼─────────────────────────
   1 │   578.0  -1024.0      0
   2 │   580.0  -1024.0      0
   3 │   582.0  -1024.0      0
   4 │   584.0  -1024.0      0
   5 │   586.0  -1024.0      0
```

Here `:M` is the magnetization, `:U` is the internal energy, `:f` is a
frequency. We can pass these arrays to the constructor for `MultihistogramData`
and it will generate the required fields for us.

```julia
julia> MHData = MultihistogramData(2, temperatures, histograms);
[ Info: Generating marginal histograms and tuple iterators
```

Next we calculate the free energies. This lattice is small, so it does not
matter if we use logsum or not, but we will to be safe.

```julia
julia> calculate_free_energies!(MHData; logsum = true)
Minimizing:  Time: 0:00:08 (980 iterations)
17-element Vector{Float64}:
  263.0645763744068
  174.58765084969377
  102.55740192250276
   43.384527986582974
   -5.480625687123756
  -45.85692706557473
  -79.05717590780614
 -106.2338328894399
 -128.57926632130915
 -147.1687162811869
 -162.82148053593113
 -176.13878462046813
 -187.57157581937713
 -197.4650510733422
 -206.0875257072719
 -213.65028480975798
 -220.32169662231635
```

We can interpolate the free energy, internal energy, and absolute magnetization
now.

```julia
julia> interpolate_free_energy_logsum(3.1415, MHData)
-245.75444371350906

julia> interpolate_energy_logsum(3.1415, MHData)
-343.5499513557231

julia> interpolate_observable_abs_logsum(3.1415, :M, MHData)
50.15985429475783
```

Note that in the last method, we have to specify which observable in case there
are multiple. Here the observable was magnetization and it's column in our
histogram dataset was `:M` so we give `:M` to this function. We can also
calculate second moments of energy and magnetization, and since while
calculating these the function calculates the first moment anyway, we can ask it
to return that to avoid duplicating our work.

```julia
julia> interpolate_energy_second_moment_logsum(3.1415, MHData; returnlinear = true)
(-343.5499513557231, 119237.72141127102)

julia> interpolate_observable_second_moment_abs_logsum(3.1415, :M, MHData; returnlinear = true)
(50.15985429475783, 3930.3703852912727)
```

And that is mostly all there is to using this library. Any suggestions or
changes are welcome ❤️.
