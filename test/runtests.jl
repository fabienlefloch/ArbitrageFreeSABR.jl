using ArbitrageFreeSABR, Test

@testset "HaganSet" begin
    forward = 1.0; expiry = 1.0;
    alpha = 0.35; beta = 0.25; nu = 1.0; rho = -0.1;
    N = 500; timesteps = 5; nd= 4;
    maturity = SABRMaturity(alpha,beta,rho,nu,forward,expiry,ArbitrageFreeSABRModel())
    density = makeTransformedDensityLawsonSwayne(maturity, N, timesteps, nd)
    @test isapprox(priceTransformedDensity(density, true, forward, ArbitrageFreeSABR.trapezoidal), 0.149701955629, atol=1e-8)
end

@testset "ImpliedVolatilitySet" begin
    price = 271.43234885190117
    strike = 275.0
    f = 356.73063159822254
    tte = 1.5917808219178082
    isCall = false
    tolerance = 1e-8
    vol = impliedVolatilityLiSORTS(price, isCall, strike, f, tte, 1.0, 0.0, tolerance,64)
    @test isapprox(blackScholesFormula(isCall, strike, f, vol*vol*tte, 1.0, 1.0), price, atol=tolerance)
end

@testset "BachelierVolatilitySet" begin
    price = 271.43234885190117
    strike = 275.0
    f = 356.73063159822254
    tte = 1.5917808219178082
    isCall = false
    tolerance = 1e-13
    vol = bachelierVolatility(price, isCall, strike, f, tte,1.0)
    @test isapprox(bachelierFormula(isCall, strike, f, vol*vol*tte, 1.0), price, atol=tolerance)
end
