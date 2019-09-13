# Package ArbitrageFreeSABR

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/fabienlefloch/ArbitrageFreeSABR.jl.svg?branch=master)](https://travis-ci.org/fabienlefloch/ArbitrageFreeSABR.jl) | [![codecov.io](http://codecov.io/github/fabienlefloch/ArbitrageFreeSABR.jl/coverage.svg?branch=master)](http://codecov.io/github/fabienlefloch/ArbitrageFreeSABR.jl?branch=master) |


Julia package to accompany the paper [Finite Difference Techniques for Arbitrage-Free SABR](https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr). The code is not meant for production purpose and does not catter for corner cases. It is only meant to illustrate the main techniques described in the paper.

##Installation
In a Julia REPL, enter `pkg` mode (by pressing `]`) and run:

```julia
(v1.0) pkg> add ArbitrageFreeSABR
```

[Julia](https://julianlang.org) version 1.0 or higher is required.

##Usage

Start by creating a `SABRMaturity` structure, which represents the parameters of the SABR model for a given option maturity.

Then make a `TransformedDensity` structure via `makeTransformedDensityLawsonSwayne`. This function will solve the arbitrage-free PDE on a grid using the specified number of space steps and time steps, and spanning `nd` standard deviations below and above. `nd=4` is a reasonable choice.

And price vanilla call and puts of a given strike by evaluating the `TransformedDensity`.


##Examples
###Price from the paper Hagan example
We use the same parameters as the example of negative density with the standard SABR
formula in (Hagan et al., 2014): α = 35%, β = 0.25, ρ = −10%, ν = 100% and forward f=1%.
The finite difference grid extends to `nd=4` standard deviations and is composed of `N=500` steps for the probability density dimension ϑ and 5 time steps.

```julia
using ArbitrageFreeSABR

forward = 1.0; expiry = 1.0;
α = 0.35; β = 0.25; ν = 1.0; ρ = -0.1;
N = 500; timesteps = 5; nd = 4;
maturity = SABRMaturity(α,β,ρ,ν,forward,expiry,ArbitrageFreeSABRModel())
density = makeTransformedDensityLawsonSwayne(maturity, N, timesteps, nd)
isCall = true; strike = forward;
priceTransformedDensity(density, isCall, strike, ArbitrageFreeSABR.midpoint)
```

The reference price is 0.149701955629

###Plot the implied density
Here, we reproduce the plot of the implied density, using the parameters of Hagan et al. (2014).

```julia
using ArbitrageFreeSABR
using Gadfly

forward = 1.0; expiry = 1.0;
α = 0.35; β = 0.25; ν = 1.0; ρ = -0.1;
N = 50; timesteps = 5; nd = 4;
maturity = SABRMaturity(α,β,ρ,ν,forward,expiry,ArbitrageFreeSABRModel())
density = makeTransformedDensityLawsonSwayne(maturity, N, timesteps, nd)
plot(x=density.zm, y=density.ϑ, Geom.line,  Guide.ylabel("Grid density"),Guide.xlabel("ϑ"))
ε = 1e-6; h = 2.0/1000; strikes = collect(0:1000) * h .+ ε;
impliedDensity = zeros(length(strikes));
price(strike) = priceTransformedDensity(density, true, strike, ArbitrageFreeSABR.midpoint)
@. impliedDensity = (price(strikes+ε)-2*price(strikes) +price(strikes-ε)) /ε^2
plot(x=strikes[2:end],y=impliedDensity[2:end], Geom.line, Guide.ylabel("Implied density"),Guide.xlabel("Strike"))
```

![Implied density of the arbitrage-free SABR model, using Hagan (2014) parameters](hagan_density.svg?raw=true "Implied density of the arbitrage-free SABR model, using Hagan (2014)")

###Implied volatility of the free-boundary SABR model
Here is an example of how to use the free-boundary SABR model instead of the more classic SABR model. We first plot the Bachelier (normal or basis point) volatities, and then the implied density, using the same parameters as Antonov et al.

```julia
using Gadfly
using ArbitrageFreeSABR

forward = 50*1e-4; expiry = 3.0;
β = 0.1; α = 0.5*forward^(1-β); ν = 0.30; ρ = -0.30;
N = 500; timesteps = 50; nd = 4;
maturity = SABRMaturity(α,β,ρ,ν,forward,expiry,FreeBoundarySABRModel())
density = makeTransformedDensityLawsonSwayne(maturity, N, timesteps, nd)

h = 5.0/1000; strikes = (collect(0:1000) * h .- 1.0)*forward ; prices = zeros(length(strikes));
price(strike) = priceTransformedDensity(density, true, strike, ArbitrageFreeSABR.midpoint)
@. prices = price(strikes)
vols = zeros(length(strikes));
vol(price, strike) = bachelierVolatility(price, true, strike, expiry, forward)
@. vols = vol(prices, strikes)
plot(x=strikes/forward,y=vols, Geom.line, Guide.ylabel("Implied volatility"),Guide.xlabel("Strike in forward units"))
```

![Implied Bachelier volatility of the free-boundary SABR model, using Antonov et al. parameters](antonov_density.svg?raw=true "Implied Bachelier volatility of the free-boundary SABR model, using Antonov et al. parameters")

```julia
ε = 1e-4; impliedDensity = zeros(length(strikes));
@. impliedDensity = (price(strikes+ε)-2*price(strikes) +price(strikes-ε)) /ε^2
plot(x=strikes[1:end]/forward,y=impliedDensity[1:end], Geom.line, Guide.ylabel("Implied density"),Guide.xlabel("Strike in forward units"))
```
![Implied density of the free-boundary SABR model, using Antonov et al. parameters](antonov_density.svg?raw=true "Implied density of the free-boundary SABR model, using Antonov et al. parameters")
### Convergence table

## Testing

In a Julia REPL session, enter `pkg` mode and run `test ArbitrageFreeSABR`.
