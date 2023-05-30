var documenterSearchIndex = {"docs":
[{"location":"apiref.html#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"apiref.html","page":"API Reference","title":"API Reference","text":"","category":"page"},{"location":"apiref.html","page":"API Reference","title":"API Reference","text":"Modules = [MultihistogramAnalysis]","category":"page"},{"location":"apiref.html#MultihistogramAnalysis.MultihistogramData","page":"API Reference","title":"MultihistogramAnalysis.MultihistogramData","text":"Primary data structure of the library. Holds all the information needed to run all of the functions provided. The constructor provided for this will precompute all the required datasets to ensure fastest interpolation. The fields are\n\nobservables::Vector{Symbol}: A vector of symbols comprising of all the observables of the system. For instance, for an Ising system, this will be the internal energy and magnetization. It is mandatory that the column for internal energy be labelled with :U. Magnetization can be labelled whatever, so for instance a valid observables vector is Symbol[:U, :M].\nparameter_values::Vector{pType}: List of parameters at which the system was simulated. For instnace, this can be the temperatures at which MCMC simulations were done.\nsimulated_histograms::Vector{DataFrame}: Flattened vector of all simulated histograms. The DataFrame must have a column for internal energy labelled :U and the frequency of occurance labelled :f. Other columns can be named whatever.\nobservable_ranges::Vector{orType}: Vector of the ranges spanned by the observables.\nmarginal_energy_histograms::Vector{ehType}: Dictionary of marginal energy histograms, mapping energy to frequency of occurance.\nmarginal_energy_histiters::Vector{ehiType}: Array of collected pairs of marginal_energy_histograms. Collected they take less memory.\nnEbins::Int64: Number of energy bins sampled in total.\nnbins::Int64: Number of histogram bins in total.\nfree_energies::Vector{pType}: Current best guess of free energies.\ntuple_iterators::{ntiType}: Named tuple iterators from the histograms for fast iterations.\n\n\n\n\n\n","category":"type"},{"location":"apiref.html#MultihistogramAnalysis.calculate_free_energies!-Tuple{MultihistogramData}","page":"API Reference","title":"MultihistogramAnalysis.calculate_free_energies!","text":"Iterate the free energy multihistogram equation until the error is below a threshold. Error is calculated as.\n\nDelta^2 = sum_i left(fracF^(m)_i - F^(m-1)_iF^(m-1)_iright)^2\n\nInputs are\n\ndata::MultihistogramData: Multihistogram dataset.\nrtol (default 1e-10): Tolerance until the error is minimized. Value of Delta.\nlogsum (default false): Use the logsum method. Recommended for systems with energies larger than 300.\n\n\n\n\n\n","category":"method"},{"location":"apiref.html#MultihistogramAnalysis.free_energy_iteration-Tuple{Any, MultihistogramData}","page":"API Reference","title":"MultihistogramAnalysis.free_energy_iteration","text":"One iteration of the multihistogram (wham) free energy equations. Inputs are\n\nu: One dimensional array-like type. Guess for the free energies.\ndata::MultihistogramData: Dataset for Multihistogram interpolation.\n\n\n\n\n\n","category":"method"},{"location":"apiref.html#MultihistogramAnalysis.free_energy_iteration_logsum-Tuple{Any, MultihistogramData}","page":"API Reference","title":"MultihistogramAnalysis.free_energy_iteration_logsum","text":"One iteration of the multihistogram (wham) free energy equations with the logsum trick. Inputs are\n\nu: One dimensional array-like type. Guess for the free energies.\ndata::MultihistogramData: Dataset for Multihistogram interpolation.\n\n\n\n\n\n","category":"method"},{"location":"usage.html#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Using this library is fairly simple. The only things to keep in mind are","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Your input datasets should be flattened histograms in the form of a DataFrame.\nEach DataFrame should have a :U column for the energy of the system and :f column for the corresponding frequency of that state. A sample dataset will be uploaded once a file hosting service is settled on.\nThis library is not thread safe. While multithreading has been implemented for the logsum methods, it may not neccesarily be faster, and there are diminishing gains as the number of threads increase.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"An example usage of this library can be seen on the following dataset. Say we simulate an Ising spin system for a number of different temperatures and output the histograms to their CSVs. Our directory structure can look something like","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"sample_dataset\n├── 1.0\n│   ├── histogram.csv\n│   └── metadata.ini\n├── 1.1\n│   ├── histogram.csv\n│   └── metadata.ini\n├── 1.2\n│   ├── histogram.csv\n│   └── metadata.ini\n├── 1.3\n│   ├── histogram.csv\n│   └── metadata.ini\n├── 1.4\n│   ├── histogram.csv\n│   └── metadata.ini","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"and so on. The directory names are the temperatures and the histogram is always stored in histogram.csv. Additionally some metadata about the simulation can also be provided, however it is not needed in using this library. We can import our dataset and create a vector of temperatures and histograms.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"julia> using DataFrames, CSV, MultihistogramAnalysis\n\njulia> datadir = \"./sample_dataset/\"\n\"./sample_dataset/\"\n\njulia> temperatures = Float64[]; histograms = DataFrame[]\nDataFrame[]\n\njulia> for simdir in filter(x -> isdir(joinpath(datadir, x)), readdir(datadir))\n           try\n               push!(temperatures, parse(Float64, simdir))\n               push!(histograms, CSV.read(joinpath(datadir, simdir, \"histogram.csv\"), DataFrame))\n           catch e\n               if e isa ArgumentError\n                   continue\n               else\n                   rethrow(e)\n               end\n           end\n       end","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Lets take a look at one of our histograms","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"julia> first(histograms[1], 5)\n5×3 DataFrame\n Row │ M        U        f     \n     │ Float64  Float64  Int64 \n─────┼─────────────────────────\n   1 │   578.0  -1024.0      0\n   2 │   580.0  -1024.0      0\n   3 │   582.0  -1024.0      0\n   4 │   584.0  -1024.0      0\n   5 │   586.0  -1024.0      0","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Here :M is the magnetization, :U is the internal energy, :f is a frequency. We can pass these arrays to the constructor for MultihistogramData and it will generate the required fields for us.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"julia> MHData = MultihistogramData(2, temperatures, histograms);\n[ Info: Generating marginal histograms and tuple iterators","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Next we calculate the free energies. This lattice is small, so it does not matter if we use logsum or not, but we will to be safe.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"julia> calculate_free_energies!(MHData; logsum = true)\nMinimizing:  Time: 0:00:08 (980 iterations)\n17-element Vector{Float64}:\n  263.0645763744068\n  174.58765084969377\n  102.55740192250276\n   43.384527986582974\n   -5.480625687123756\n  -45.85692706557473\n  -79.05717590780614\n -106.2338328894399\n -128.57926632130915\n -147.1687162811869\n -162.82148053593113\n -176.13878462046813\n -187.57157581937713\n -197.4650510733422\n -206.0875257072719\n -213.65028480975798\n -220.32169662231635","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"We can interpolate the free energy, internal energy, and absolute magnetization now.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"julia> interpolate_free_energy_logsum(3.1415, MHData)\n-245.75444371350906\n\njulia> interpolate_energy_logsum(3.1415, MHData)\n-343.5499513557231\n\njulia> interpolate_observable_abs_logsum(3.1415, :M, MHData)\n50.15985429475783","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Note that in the last method, we have to specify which observable in case there are multiple. Here the observable was magnetization and it's column in our histogram dataset was :M so we give :M to this function. We can also calculate second moments of energy and magnetization, and since while calculating these the function calculates the first moment anyway, we can ask it to return that to avoid duplicating our work.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"julia> interpolate_energy_second_moment_logsum(3.1415, MHData; returnlinear = true)\n(-343.5499513557231, 119237.72141127102)\n\njulia> interpolate_observable_second_moment_abs_logsum(3.1415, :M, MHData; returnlinear = true)\n(50.15985429475783, 3930.3703852912727)","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"And that is mostly all there is to using this library. Any suggestions or changes are welcome ❤️.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"CurrentModule = MultihistogramAnalysis","category":"page"},{"location":"index.html#MultihistogramAnalysis","page":"Home","title":"MultihistogramAnalysis","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Documentation for MultihistogramAnalysis.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"warning: Warning\nAt the moment the library will only work with systems that have a discreet energy spectra. Generation of histograms neccesarily means that we bin data, and thus throw away data points if the energy spectra is not discreet. Use with caution on systems with continuous energy spectra.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The multihistogram method solves for the density of states by solving for the partition functions at the temperatures simulated at. The main equation that is solved iteratively is.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Z(beta_k) = sum_isum_E fracN_i(E)sum_j n_j Z^-1(beta_j) e^(beta_k\n- beta_j)E","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Where ij indexes the simulation number, E are the energies sampled in the i-th simulation, N_i is the energy histogram (i.e, number of times an energy occurs in the simulation), and n_i are the number of samples in the i-th histogram. For this library, n_i = 0 for all i. Since the partition functions depend exponentially on the energy, it is better to use the logarithms (F(beta_k) = F_k = log(Z(beta_k)) of the partition functions to avoid float overflow. Thus the equation becomes.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"F_k =logleft( sum_i sum_E fracN_i(E)sum_j log(-F_k + (beta_k - beta\n)E right)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"To avoid overflow further, we normalize all partition functions by the following factor.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"A = frac1sqrtZ_maxZ_min implies log(A) = -frac12 (F_max +\nF_min)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"To interpolate the free energy we use a natural extension of the above equation.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"F(beta) = logleft( sum_i sum_E fracN_i(E)sum_j log(-F_k + (beta -\nbeta)E right)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"And the calculate the thermal average we use the following.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"langle Q rangle = frac1Z(beta) sum_i (Q E) Q fracN_i(Q\nE)sum_j A Z^-1_j e^(beta - beta_j)E","category":"page"},{"location":"index.html#Logsum-Methods","page":"Home","title":"Logsum Methods","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The logsum trick allows us to calculate the logarithm of a sum of a series of logarithms that we want to exponentiate. Specifically, if l_i are logarithms, and we want to calculate log(sum_i e^l_i) without ever exponentiating the logarithms, this trick lets us do that. If l_0 is the largest logarithm, then one can calculate the sum as follows.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"log(sum_i=0 exp(l_i)) = log(e^l_0 left1 + sum_i=1 e^l_i -\nl_0right) = l_0 + log(1 +  sum_i=1 e^l_i - l_0)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This has the advantage of avoiding numerical overflow and calculate observables for very large systems. A function name that ends with _logsum uses this method. Note that this has also been implemented for ComplexF64 vectors.","category":"page"},{"location":"index.html#Interpolatioon-Functions","page":"Home","title":"Interpolatioon Functions","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"After importing the library, you can write interpolate_ and try to autocomplete with TAB to see all the interpolation methods. Counterintuatively, the MultihistogramData structs go as the last argument rather than the first. I did not have a lot of time to write this library cleanly, will be fixed in the future.","category":"page"}]
}
