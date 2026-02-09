# Example

```@example ex1
using SpatialBernoulli 
using Random
using OptimizationOptimJL
rng = MersenneTwister(1234)
solver = Optim.LBFGS()
```
## Parameters

```@example ex1



# define locations in the unit square
my_locations = vcat(([x y] for x in 0:0.2:1 for y in 0:0.2:1)...)
nlocs = length(my_locations[:, 1])


my_distance = [sqrt(sum(abs2, my_locations[i, :] - my_locations[j, :])) for i in axes(my_locations, 1), j in axes(my_locations, 1)]

# randomly generate a SpatialBernoulli
my_λ = rand(nlocs)# [0.5+ 0.1*i for i in 1:nlocs] 
my_range = 0.3
my_sill = 1.0
my_order = 1 / 2

```

## Create the distribution

```@example ex1
d = SB(my_range, my_sill, my_order, my_λ, my_distance)
```

use it like any Distributions.jl distribution



## Random draws

```@example ex1

n = 2000
y=rand(rng,d,n)

using CairoMakie
yplot=y[:,2]
begin
fig = Figure()
ax = Axis(fig[1, 1])

scatter!(
    ax,
    my_locations[yplot .== 0, 1],
    my_locations[yplot .== 0, 2];
    color = :black,
    label = "y = 0"
)

scatter!(
    ax,
    my_locations[yplot .== 1, 1],
    my_locations[yplot .== 1, 2];
    color = :white,
    strokecolor = :black,
    label = "y = 1"
)

axislegend(ax)
fig
end
```

## Fitting 
```@example ex1
init_range = 0.5
init_order = 1.5
init_lambda = fill(0.4, nlocs)
init_d = SB(init_range, 1.0, init_order, init_lambda, my_distance)

tdist = maximum(my_distance) / 1
wp = 1.0 .* (my_distance .< tdist)

@timed sol3 = fit_mle_vfast(init_d, y, wp;solver, order=my_order, return_sol=true, maxiters = 2000)

```