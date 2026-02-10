# SpatialBernoulli package

*[SpatialBernoulli.jl](https://github.com/caroline-cognot/SpatialBernoulli.jl)* is a package to define a spatially correlated Bernoulli variable `Y` using a latent Gaussian construction.

## Model

The Spatial Bernoulli model is defined as follows:

- `Y ~ SB(C_Y, λ)`
- `(X_{Y,1}, …, X_{Y,D}) ~ N(0, C_Y)`
- For all spatial locations `s`:

  - `Y_s = 1` if `X_{Y,s} ≤ Φ⁻¹(λ_s)`
  - `Y_s = 0` otherwise

Here:

- `C_Y` is the covariance matrix of the latent Gaussian field
- `λ = (λ_s)_{s=1,…,D}` is the vector of marginal probabilities
- `Φ` denotes the cumulative distribution function (CDF) of the standard normal distribution

This model is used in the paper *to be inserted later* to model precipitation occurrence across a large region.

## Features

The package provides several methods, including:

- Model definition
- Probability density functions (`pdf`) and their logarithm (`logpdf`), computed using bivariate Gaussian integrals with  
  *[MvNormalCDF.jl](https://github.com/PharmCat/MvNormalCDF.jl)*
- Maximum likelihood estimation (not recommended for high-dimensional samples)
- Maximum pairwise likelihood estimation with bivariate integrals computed using  
  *[MvNormalCDF.jl](https://github.com/PharmCat/MvNormalCDF.jl)*
- Fast maximum pairwise likelihood estimation with bivariate integrals computed using the approximation of  
  [Tsay (2023)](https://doi.org/10.1080/03610918.2021.1884718)

## Example

The documentation also provides a simulated example.

Below are examples of generation from the Spatial Bernoulli model with:

- a constant marginal probability `λ` in space, and
- an exponential covariance function  
  `C_Y(h) = exp(-h / ρ)`

![Model illustration](assets/6fd2db61-1.png)






