using ArbitrageFreeSABR, Test

@testset "HaganSet" begin
    forward = 1.0; expiry = 1.0;
    alpha = 0.35; beta = 0.25; nu = 1.0; rho = -0.1;
    N = 500; timesteps = 5; nd= 4;
    maturity = SABRMaturity(alpha,beta,rho,nu,forward,expiry)
    density = makeTransformedDensityLawsonSwayne(maturity, N, timesteps, nd)
    @test isapprox(priceTransformedDensity(density, true, forward, ArbitrageFreeSABR.trapezoidal), 0.149701955629, atol=1e-8)
end
