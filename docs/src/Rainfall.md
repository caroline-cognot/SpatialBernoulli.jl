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
    months=1:12
month_labels = ["Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec"]
    ρ = [m.range for m in models]
    λ = hcat([m.λ for m in models]...)  # (n_sites × n_months)

    fig = Figure(size = (900, 500))
    ax1 = Axis(fig[1, 1],
        xlabel = "Month",
        ylabel = L"\rho",
        title  = L"Estimated range parameter $\rho$",
            xticks = (months, month_labels)

    )
    scatter!(ax1, months, ρ)
    lines!(ax1, months, ρ)

    ax2 = Axis(fig[1, 2],
        xlabel = "Month",
        ylabel = L"\lambda_s",
        title  = L"Estimated rain probabilities $\lambda_s$",
            xticks = (months, month_labels)

    )
    for s in 1:size(λ, 1)
        lines!(ax2, months, λ[s, :], linewidth = 1, alpha = 0.6)
    end

    fig
end




plot_parameters(vec_models_vfast)

```


## Visualise the results  : comparing simulations to observations

The goal of stochastic weather generators is to be able to simulate quickly many plausible sequences of the meteorological variables, sharing the statistical properties of the observations. 

An objective of this secific model is to accurately reproduce large-scale spatial events, such as widespread dry or wet days.
We define the `ROR` indicator as
    `ROR(n) = \frac{1}{D}\sum_{s=1}^D Y_s^{(n)}`
A low (high) value of `ROR(n)` denotes a dry (wet) day for many locations at the same time. 

The distribution in the observations is compared to the distribution evaluated from many simulations, as well as the autocorrelation of this indicator, for each season.

Unsuprisingly, there is a clear lack of temporal dependance when using a simple SpatialBernoulli, suggesting the need for a structured model (such as HMMs! ). However, the distribution of the indicator itself is not that far off.

```@example rainfall
Nb = 10
D = 37
using StatsBase
using LinearAlgebra, NaNMath

begin
    Ysvf = zeros(Bool, nlocs, ntot, Nb)
    @time "Simulations  Y" for i in 1:Nb
        for t in 1:ntot
        m = month(every_year[t])
        Ysvf[:, t, i] = rand(vec_models_vfast[m]);
        end
    end
end

RRmax = 0
RORo = [mean(r .> RRmax) for r in eachcol(Yobs)]
RORs = [[mean(r .> RRmax) for r in eachcol(rr)] for rr in eachslice(Ysvf, dims=3)]

maxlag = 10

include("assets/utilities.jl")
begin
    # Makie ROR distribution and autocorrelation (2×4 grid)
JJA = [6, 7, 8]
MAM = [3, 4, 5]
SON = [9, 10, 11]
DJF = [12, 1, 2]
SEASONS = [DJF, MAM, JJA, SON]
seasonname = ["DJF", "MAM", "JJA", "SON"]
idx_seasons = [findall(month.(every_year) .∈ tuple(season)) for season in SEASONS]
    fig_ROR = Figure(fontsize=19)
    wwww = 200
    hhhh = 150
    # Row 1: Distribution plots
    for m in eachindex(idx_seasons)
        row = 1
        col = m

        # Distribution subplot
        ax_dist = Axis(fig_ROR[row, col],
            xlabel="ROR",
            ylabel=col == 1 ? "Distribution" : "",
            title=seasonname[m],
            width=wwww,
            height=hhhh)
        xax = 0:(1/D):1.0
        xaxbin = vcat(xax, [1.01])
        errorlinehist!(ax_dist, [RORs[i][idx_seasons[m]] for i in 1:Nb];
            label="",
            color=:gray,
            secondarycolor=:gray, normalization=:probability,
            bins=xaxbin,
            errortype=:percentile,
            percentiles=[0, 100],
            alpha=0.5,
            secondaryalpha=0.2,
            centertype=:median)
     
        ylims!(ax_dist, 0, 0.06)
        col > 1 && hideydecorations!(ax_dist, grid=false)

        # Autocorrelation subplot
        row = 2
        ax_acf = Axis(fig_ROR[row, col],
            xlabel="Lag",
            ylabel=col == 1 ? "ACF" : "",
            width=wwww,
            height=hhhh,
        )

        rorsim = [RORs[i][idx_seasons[m]] for i in 1:Nb]
        acf_sim = [autocor(rorsim[i], 0:maxlag) for i in 1:length(rorsim)]
        errorline!(ax_acf, 0:maxlag, stack(acf_sim, dims=1)',
            color=:gray,
            secondarycolor=:gray,
            errortype=:percentile,
            percentiles=[0, 100],
            secondaryalpha=0.2,
            centertype=:median)
    
        # Observations
        acf_obs = autocor(RORo, 0:maxlag)
        scatter!(ax_acf, 0:maxlag, acf_obs, color=:blue, markersize=7)
        col > 1 && hideydecorations!(ax_acf, grid=false)
    end
 
    Legend(
    fig_ROR[:, 5],
    [
        [LineElement(color=:gray), PolyElement(color=:gray, alpha=0.2)],

        MarkerElement(color=:blue, marker=:circle, markersize=8)
    ],
    [
        "SB",
        "Observations"
    ]
)

    resize_to_layout!(fig_ROR)
    fig_ROR
end
```



