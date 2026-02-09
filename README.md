# SpatialBernoulli

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

# MyExtendedExtremes.jl

A dummy package for learning how to make packages
It uses the very real package Extended Extremes, and create a custom distribution.

[![CI](https://github.com/caroline-cognot/MyExtendedExtremes.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/caroline-cognot/MyExtendedExtremes.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://caroline-cognot.github.io/MyExtendedExtremes.jl/dev)


## Installation

The package runs on julia 1.11 and above.
In a Julia session switch to `pkg>` mode to add the package:

```julia
julia>] # switch to pkg> mode
pkg> add https://github.com/caroline-cognot/MyExtendedExtremes.jl
```
To run an example 

```julia

using Distributions
using MyExtendedExtremes 

uniform_dist = Uniform(0.0, 1.0)
tail_dist = Exponential(1.0)
mixed_dist = MixedUniformTail(0.3, uniform_dist, tail_dist, 0.1, 0.5)

```
