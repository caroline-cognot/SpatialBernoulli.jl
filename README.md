# SpatialBernoulli.jl

A package to define Spatial Bernoulli as 

$$ {Y} \sim \mathcal{SB}(C_{Y},{\lambda})$$

$$ {X_{Y}}_{1},...,{X_{Y}}_{D} \sim \mathcal{N}(0,C_{Y}) $$

$$ \forall s,~  Y_s =  
\begin{cases} 
    1 & \text{if } {X_{Y}}_{s} \leq  \Phi^{-1}(\lambda_s)  \\ 
    0 & \text{else} 
\end{cases} $$
with $\Phi$ the cumulative distribution function (CDF) of the standard normal distribution.$$

```julia
# Function to generate binary field given λ and range parameter
function generate_binary_field(my_λ, ρ)
    my_sill, my_order = 1.0, 0.5
    d = SpatialBernoulli(ρ, my_sill, my_order, my_λ, my_distance)
    u = rand(MvNormal(d.ΣU))
    thresholds = quantile.(Normal(), d.λ)
    ys = u .< thresholds
    return reshape(ys, dy, dx)
end


"""
    make_squares(λs, ρs, fields)
Binary fields for different ρ, λ values and fields
"""
function make_square(λs, ρs, fields, xg, yg)
    fig = Figure()
    axs = [Axis(fig[j, i], width=120, height=120, limits=(0, 2, 0, 2)) for (j, λ) in enumerate(λs), (i, ρ) in enumerate(ρs)]
    for (j, λ) in enumerate(λs)
        for (i, ρ) in enumerate(ρs)
            heatmap!(axs[j, i], xg, yg, fields[j, i], colormap=:binary, colorrange=(0, 1))
            if j < length(λs)
                hidexdecorations!(axs[j, i])
            end
            if i > 1
                hideydecorations!(axs[j, i])
            end
        end
    end
    for (j, λ) in enumerate(λs)
        Label(fig[j, 0], L"\lambda = %$(λ[1])", tellheight=false)
    end
    for (i, ρ) in enumerate(ρs)
        Label(fig[0, i], L"\rho = %$(ρ)", tellwidth=false)
    end
    Label(fig, bbox=BBox(-115, 200, 50, 300), L"$Y=1: \blacksquare$\n$Y=0: □$")
    colgap!(fig.layout, (25))
    resize_to_layout!(fig)
    # end
    return fig
end

# Grid
xg = 0:0.03:2
yg = 0:0.03:2
dx, dy = length(xg), length(yg)
my_locations = vcat(([xx yy] for xx in xg for yy in yg)...)
my_distance = [sqrt(sum(abs2, my_locations[i, :] - my_locations[j, :]))
               for i in axes(my_locations, 1), j in axes(my_locations, 1)]
nlocs = size(my_locations, 1)

λs = [0.2, 0.5]
ρs = [0.01, 0.1, 0.5]
fields = [generate_binary_field(fill(λ, nlocs), ρ) for λ in λs, ρ in ρs]
fig_row = make_square(λs, ρs, fields, xg, yg)

```

[![CI](https://github.com/caroline-cognot/SpatialBernoulli.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/caroline-cognot/SpatialBernoulli.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://caroline-cognot.github.io/SpatialBernoulli.jl/dev)


## Installation

The package runs on julia 1.11 and above.
In a Julia session switch to `pkg>` mode to add the package:

```julia
julia>] # switch to pkg> mode
pkg> add https://github.com/caroline-cognot/SpatialBernoulli.jl
```
To run an example , see documentation.


## Automatically generated  
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caroline-cognot
.github.io/SpatialBernoulli.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caroline-cognot
.github.io/SpatialBernoulli.jl/dev/)
[![Build Status](https://github.com/caroline-cognot
/SpatialBernoulli.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/caroline-cognot
/SpatialBernoulli.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/caroline-cognot
/SpatialBernoulli.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/caroline-cognot
/SpatialBernoulli.jl)
