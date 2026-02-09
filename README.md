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

<p align="center">
![logo](https://www.raspberrypi.org/app/uploads/2018/03/RPi-Logo-Reg-SCREEN-199x250.png "Raspberry pi")
</p>

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
