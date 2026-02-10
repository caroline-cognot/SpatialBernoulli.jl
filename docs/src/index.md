# SpatialBernoulli package


*[SpatialBernoulli.jl](https://github.com/caroline-cognot/SpatialBernoulli.jl)* is a package to define a spatially correlated Bernoulli variable $Y$ using the following construction :


$$ {Y} \sim \mathcal{SB}(C_{Y},{\lambda})$$

$$ {X_{Y}}_{1},...,{X_{Y}}_{D} \sim \mathcal{N}(0,C_{Y}) $$

$$ \forall s,~  Y_s =  
\begin{cases} 
    1 & \text{if } {X_{Y}}_{s} \leq  \Phi^{-1}(\lambda_s)  \\ 
    0 & \text{else} 
\end{cases} $$
with $\Phi$ the cumulative distribution function (CDF) of the standard normal distribution.

It is used in the paper *to be inserted later* to model precipitation occurrence across a large region.

The package provides several methods, such as 
* Model definition 
* Probability density functions (pdf) and their logarithm (logpdf) computed using with bivariate integrals computed using *[MvNormalCDF.jl](https://github.com/PharmCat/MvNormalCDF.jl)*
* Maximum likelihood estimation (not recommended for high-dimensional samples) 
* Maximum pairwise likelihood estimation with bivariate integrals computed using  *[MvNormalCDF.jl](https://github.com/PharmCat/MvNormalCDF.jl)*
* Fast maximum pairwise likelihood estimation with bivariate integrals computed using the approximation of [Tsay (2023)](https://doi.org/10.1080/03610918.2021.1884718)


The documentation also provides a simulated example. 

Below are examples of generation from the Spatial Bernoulli with constant $\lambda$ in space and the exponential covariance function $C_Y(h)=\exp(-h/\rho)$.

![Model illustration](assets/6fd2db61-1.png)







