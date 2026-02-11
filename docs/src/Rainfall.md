# Modelling the rainfall occurrence in France

Let us use the `SB` distribution to model rain occurrence across France.
The data used is extracted from the [ECAD](https://www.ecad.eu/dailydata/index.php) database. 

It consists in daily rainfall occurrence at D=37 locations in France.


```@example rainfall
using SpatialBernoulli, Optimization,OptimizationOptimJL, Optim
using CSV,DataFrames
using Random
import ForwardDiff
Random.seed!(1234)

```

## Getting the data
```@example rainfall
station_50Q = CSV.read("data/transformedECAD_stations.csv",DataFrame);
Yobs=Matrix(CSV.read("data/transformedECAD_Yobs.csv",header=false,DataFrame));
my_distance =Matrix(CSV.read("data/transformedECAD_locsdistances.csv",header=false,DataFrame));
```

First, we split the data in months. Explicit seasonality in the sense of periodic parameterisation is not treated here, see [here](https://github.com/caroline-cognot/HMM_Spatial_codepaper) for using this type of model in a periodic setting, or [here](https://github.com/dmetivie/SmoothPeriodicStatsModels.jl) for other parameterisations.

```@example rainfall

my_locations = hcat(station_50Q.LON_idx, station_50Q.LAT_idx)
nlocs = length(my_locations[:, 1])

select_month = function (m::Int64, dates, Y::AbstractMatrix)
    indicesm = findall(month.(dates) .== m)
    return Y[:, indicesm]
end

using Dates
date_start = Date(1973)


date_end = Date(2024) - Day(1)

every_year = date_start:Day(1):date_end
dates = every_year


Ymonths = [select_month(m, every_year, Yobs) for m in 1:12];
```

## Inference for each month
We can now estimate model parameters for each month, using the `fit_mle_vfast` method. The results are saved to the `data` folder. This is easily multi-thread-able.

```@example rainfall
using JLD2
using LineSearches

for imonth in 1:12
    y = Ymonths[imonth]
    n = length(y[1, :])


    init_range = 600
    init_order = 0.5
    init_lambda = fill(0.4, nlocs)
    init_d = SB(init_range, 1.0, init_order, init_lambda, my_distance)
   
    #pairwise solution
    tdist = maximum(my_distance) / 3
    wp = 1.0 .* (my_distance .< tdist)

    @time sol_fixednu = fit_mle_vfast(init_d, y, wp;solver = Optim.LBFGS(
    linesearch = LineSearches.BackTracking()
)
,  order=1/2,maxiters = 100, return_sol=false)

    # Save to JLD2 file
    save("data/vfastfitted_month_QMC100" * string(imonth) * ".jld2", Dict("d" => sol_fixednu))
    
end
```

## Visualise the results  : parameters 

```@example rainfall

vec_models_vfast = Vector{SB}(undef, 12)
for imonth in 1:12
    vec_models_vfast[imonth] = load("data/vfastfitted_month_QMC100" * string(imonth) * ".jld2")["d"]
end
ntot = length(Yobs[1,:])
using CairoMakie

function plot_parameters(models::Vector{SB})
    nmonths = length(models)

    ρ = [m.range for m in models]
    λ = hcat([m.λ for m in models]...)  # (n_sites × n_months)

    fig = Figure(size = (900, 500))

    # --- Panel 1: range parameter ρ ---
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Month",
        ylabel = L"\rho",
        title  = L"Estimated range parameter $\rho$"
    )

    scatter!(ax1, 1:nmonths, ρ; markersize = 10)
    lines!(ax1, 1:nmonths, ρ)

    # --- Panel 2: site probabilities λ_s ---
    ax2 = Axis(
        fig[1, 2],
        xlabel = "Month",
        ylabel = L"\lambda_s",
        title  = L"Estimated rain probabilities $\lambda_s$"
    )

    for s in 1:size(λ, 1)
        lines!(ax2, 1:nmonths, λ[s, :], linewidth = 1, alpha = 0.7)
    end

    fig
end



end
```

```@example rainfall
p1,p2=Plot( vec_models_vfast)
Plots.plot(p1,p2)
```

## 

# there is probably a better way to do this ?
Nb = 500
nlocs= length(my_locations[:,1])
begin
    Ys = zeros(Bool, nlocs, ntot, Nb)
    @time "Simulations  Y" for i in 1:Nb
        for t in 1:ntot
        m = month(every_year[t])
        Ys[:, t, i] = rand(vec_models[m]);
        end
    end
end

begin
    Ysvf = zeros(Bool, nlocs, ntot, Nb)
    @time "Simulations  Y" for i in 1:Nb
        for t in 1:ntot
        m = month(every_year[t])
        Ysvf[:, t, i] = rand(vec_models_vfast[m]);
        end
    end
end


Ys
Ysvf
Yobs

include("../11SpatialBernoulli/plot_validation.jl")
p=compare_ROR_density(Yobs,Ys)
savefig(p, "./11SpatialBernoulli/ROR_500sim.png")

p=compare_ROR_histogram(Yobs,Ys)
Plots.plot!(p, title= "Rain Occurrence Ratio for 37 stations, 500 simulations",size=(1000,500))
savefig(p, "./11SpatialBernoulli/ROR_500sim_histo.png")


p=compare_ROR_histogram(Yobs,Ysvf)
Plots.plot!(p, title= "Rain Occurrence Ratio for 37 stations, 500 simulations",size=(1000,500))
savefig(p, "./11SpatialBernoulli/vfastROR_500sim_histo.png")



savefig(pp, "./11SpatialBernoulli/vfastROR_500sim_andparams.png")
