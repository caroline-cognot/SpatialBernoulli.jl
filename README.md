[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caroline-cognot.github.io/SpatialBernoulli.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caroline-cognot.github.io/SpatialBernoulli.jl/dev/)
[![Build Status](https://github.com/caroline-cognot/SpatialBernoulli.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/caroline-cognot/SpatialBernoulli.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/caroline-cognot/SpatialBernoulli.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/caroline-cognot/SpatialBernoulli.jl)


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

Below are examples of generation from the Spatial Bernoulli with constant $\lambda$ in space and the exponential covariance function $C_Y(h)=\exp(-h/\rho)$.

![Model illustration](docs/src/assets/6fd2db61-1.png
)



## Installation

The package runs on julia 1.11 and above.
In a Julia session switch to `pkg>` mode to add the package:

```julia
julia>] # switch to pkg> mode
pkg> add https://github.com/caroline-cognot/SpatialBernoulli.jl
```
To run an example, see documentation.

For now, model fitting relies on Optimization.


## How did I build this package ?

Followed instructions from
- https://julialang.org/contribute/developing_package/ for the creation of the package using PkgTemplates
- https://documenter.juliadocs.org/stable/man/guide/ for setting up the documentation
- https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/ for publishing the documentation

