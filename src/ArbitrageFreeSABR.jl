module ArbitrageFreeSABR
using LinearAlgebra

include("BlackScholes.jl")
include("ImpliedVolatility.jl")
include("BachelierVolatility.jl")

"Structure for the parameters of the SABR model for a given option maturity"
abstract type SABRModel end
struct ArbitrageFreeSABRModel <: SABRModel end
struct FreeBoundarySABRModel <: SABRModel end


struct SABRMaturity{T<:Real,M<:SABRModel}
  α::T
  β::T
  ρ::T
  ν::T
  forward::T
  time::T
  model::M
end


"Probability density in the finite difference grid, at a given time, expressed in the transformed coordinate theta"
struct TransformedDensity{T<:Real, M<:SABRModel}
  maturity::SABRMaturity{T,M}
  zmin::T
  zmax::T
  h::T
  zm::Vector{T}
  ϑ::Vector{T}
  PL::T
  PR::T
end

export ArbitrageFreeSABRModel, FreeBoundarySABRModel, SABRMaturity, makeTransformedDensityLawsonSwayne, priceTransformedDensity

function solveStep(ϑ::Vector{G}, PL::G, PR::G, h::G, Fm::Vector{G}, Cm::Vector{G}, Em::Vector{G}, dt::G) where {G}
  frac = dt / (2 * h)
  M = length(ϑ)
  B = zeros(M)
  A = zeros(M - 1)
  C = zeros(M - 1)
  @. B[2:M-1] = 1.0 + frac * (Cm[2:M-1] * Em[2:M-1] * (1.0 / (Fm[3:M] - Fm[2:M-1]) + 1.0 / (Fm[2:M-1] - Fm[1:M-2])))
  @. C[2:M-1] = -frac * Cm[3:M] * Em[3:M] / (Fm[3:M] - Fm[2:M-1])
  @. A[1:M-2] = -frac * Cm[1:M-2] * Em[1:M-2] / (Fm[2:M-1] - Fm[1:M-2])
  B[1] = Cm[1] / (Fm[2] - Fm[1]) * Em[1]
  C[1] = Cm[2] / (Fm[2] - Fm[1]) * Em[2]
  B[M] = Cm[M] / (Fm[M] - Fm[M-1]) * Em[M]
  A[M-1] = Cm[M-1] / (Fm[M] - Fm[M-1]) * Em[M-1]
  tri = Tridiagonal(A, B, C)
  p1 = ϑ[1]
  pm = ϑ[M]
  ϑ[1] = 0
  ϑ[M] = 0
  ϑ0 = tri \ ϑ
  PL0 = PL + dt * Cm[2] / (Fm[2] - Fm[1]) * Em[2] * ϑ0[2]
  PR0 = PR + dt * Cm[M-1] / (Fm[M] - Fm[M-1]) * Em[M-1] * ϑ0[M-1]
  ϑ[1] = p1
  ϑ[M] = pm
  return ϑ0, PL0, PR0
end

function computeBoundaries(m::SABRMaturity{G,ArbitrageFreeSABRModel}, nd) where {G}
  zmin = -nd * sqrt(m.time)
  zmax = -zmin
  if (m.β < 1)
    ybar = -m.forward^(1 - m.β) / (1 - m.β)
    zbar = -1 / m.ν * log((sqrt(1 - m.ρ^2 + (m.ρ + m.ν * ybar / m.α)^2) - m.ρ - m.ν * ybar / m.α) / (1 - m.ρ))
    if (zbar > zmin)
      zmin = zbar
    end
  end
  return zmin, zmax
end

function computeBoundaries(m::SABRMaturity{G,FreeBoundarySABRModel}, nd) where {G}
  zmin = -nd * sqrt(m.time)
  zmax = -zmin
  return zmin, zmax
end

function yofz(m::SABRMaturity{G,M}, zm) where {G,M <: SABRModel}
  return m.α / m.ν * (sinh.(m.ν * zm) + m.ρ * (cosh.(m.ν * zm) .- 1))
end
function fofy(m::SABRMaturity{G,ArbitrageFreeSABRModel}, ym) where {G}
  return (m.forward^(1 - m.β) .+ (1 - m.β) * ym).^(1 / (1 - m.β))
end
function fofy(m::SABRMaturity{G,FreeBoundarySABRModel}, ym) where {G}
  u = sign(m.forward) * abs(m.forward)^(1 - m.β) .+ (1 - m.β) * ym
  return sign.(u) .* ((abs.(u)) .^ (1 / (1 - m.β)))
end
function cofy(m::SABRMaturity{G,ArbitrageFreeSABRModel}, ym, Fm) where {G}
  return sqrt.(m.α^2 .+ 2 * m.ρ * m.α * m.ν * ym + m.ν^2 * ym.^2) .* Fm.^(m.β)
end
function cofy(m::SABRMaturity{G,FreeBoundarySABRModel}, ym, Fm) where {G}
  return sqrt.(m.α^2 .+ 2 * m.ρ * m.α * m.ν * ym + m.ν^2 * ym.^2) .* abs.(Fm) .^ (m.β)
end
function gammaOfF(m::SABRMaturity{G,ArbitrageFreeSABRModel}, Fm, j0::Int) where {G}
  Gamma = (Fm.^m.β .- m.forward^m.β) ./ (Fm .- m.forward)
  Gamma[j0+1] = m.β / m.forward.^(1 - m.β)
  return Gamma
end

function gammaOfF(m::SABRMaturity{G,FreeBoundarySABRModel}, Fm, j0::Int) where {G}
  Gamma = (abs.(Fm) .^ m.β .- abs(m.forward)^m.β) ./ (Fm .- m.forward)
  Gamma[j0+1] = sign(m.forward) * m.β / abs(m.forward).^(1 - m.β)
  return Gamma
end
function yOfStrike(m::SABRMaturity{T, ArbitrageFreeSABRModel}, strike) where {T}
  return (strike.^(1 - m.β) .- m.forward^(1 - m.β)) / (1 - m.β)
end

function yOfStrike(m::SABRMaturity{T, FreeBoundarySABRModel}, strike) where {T}
  return (sign.(strike) .* abs.(strike).^(1 - m.β) .- m.forward^(1 - m.β)) / (1 - m.β)
end

function makeForward(m::SABRMaturity{T, M}, z) where {T, M <: SABRModel}
  y = yofz(m,z)
  return fofy(m,y)
end

function makeTransformedDensityLawsonSwayne(maturity::SABRMaturity{G,M}, N::Int, timesteps::Int, nd::Real) where {G, M <: SABRModel}
  α = maturity.α
  β = maturity.β
  ν = maturity.ν
  ρ = maturity.ρ
  forward = maturity.forward
  T = maturity.time
    # initialize h,F and Q
  zmin, zmax = computeBoundaries(maturity, nd)
  J = N - 2
  h0 = (zmax - zmin) / J
  j0 = ceil((0 - zmin) / h0)
  h = (0 - zmin) / (j0 - 0.5)
  z = collect(0:(J+1)) * h .+ zmin
  zmax = z[J+2]
  zm = z .- 0.5 * h
  ym = yofz(maturity, zm)
  ym[1] = ym[2] #avoid negative
  Fm = fofy(maturity, ym)
  Cm = cofy(maturity, ym, Fm)
  Gammam = gammaOfF(maturity, Fm, Int(j0))
  ymax = yofz(maturity, zmax)
  ymin = yofz(maturity, zmin)
  Fmax = fofy(maturity, ymax)
  Fmin = fofy(maturity, ymin)
  Fm[1] = 2 * Fmin - Fm[2]
  Fm[J+2] = 2 * Fmax - Fm[J+1]
  Cm[1] = Cm[2]
  Cm[J+2] = Cm[J+1]
  dt = T / timesteps
    #Ptotal = sum(h*ϑ[2:J+1])+PL+PR
    #Ftotal = Fm[2:J+1]*ϑ[2:J+1]*h+Fmin*PL+Fmax*PR
  b = 1 - sqrt(2) / 2 #Lawson Swayne param
  dt1 = dt * b
  dt2 = dt * (1 - 2 * b)
  Em = ones(J + 2)
  Emdt1 = exp.(ρ * ν * α * Gammam * dt1)
  Emdt1[1] = Emdt1[2]
  Emdt1[J+2] = Emdt1[J+1]
  Emdt2 = exp.(ρ * ν * α * Gammam * dt2)
  Emdt2[1] = Emdt2[2]
  Emdt2[J+2] = Emdt2[J+1]
  PL = 0.0
  PR = 0.0
  ϑ = zeros(J + 2)
  ϑ[Int(j0)+1] = 1.0 / h
  for t = 1:timesteps
    @. Em *= Emdt1
    ϑ1, PL1, PR1 = solveStep(ϑ, PL, PR, h, Fm, Cm, Em, dt1)
    @. Em *= Emdt1
    ϑ2, PL2, PR2 = solveStep(ϑ1, PL1, PR1, h, Fm, Cm, Em, dt1)
    @. ϑ = (sqrt(2) + 1) * ϑ2 - sqrt(2) * ϑ1
    PL = (sqrt(2) + 1) * PL2 - sqrt(2) * PL1
    PR = (sqrt(2) + 1) * PR2 - sqrt(2) * PR1
    @. Em *= Emdt2
      #Ptotal = sum(h*ϑ[2:J+1])+PL+PR
      #Ftotal = Fm[2:J+1]*ϑ[2:J+1]+h+Fmin*PL+Fmax*PR
  end
  density = TransformedDensity(maturity, zmin, zmax, h, zm, ϑ, PL, PR)
  return density
end

@enum DensitySmoothing none = 1 trapezoidal = 2 midpoint = 3

function priceTransformedDensity(
  density::TransformedDensity,
  isCall::Bool,
  strike::G,
  smoothing::DensitySmoothing
) where {G}
  maturity = density.maturity
  α = maturity.α
  β = maturity.β
  ν = maturity.ν
  ρ = maturity.ρ
  forward = maturity.forward
  ystrike = yOfStrike(maturity, strike)
  zstrike = -1 / ν * log.((sqrt.(1 .- ρ^2 .+ (ρ + ν * ystrike / α).^2) .- ρ - ν * ystrike / α) / (1 - ρ))
  if isCall
    if (zstrike <= density.zmin)
      p = forward .- strike
    else
      if (zstrike >= density.zmax)
        p = 0
      else
        Fmax = makeForward(maturity, density.zmax)
        p = (Fmax .- strike) * density.PR
        k0 = Int.(ceil.((zstrike .- density.zmin) / density.h))
        ztilde = density.zmin .+ k0 * density.h
        ftilde = makeForward(maturity, ztilde)
        term = ftilde .- strike
        Fm = makeForward(maturity, density.zm[k0+1:end]) #zm[k0+1] = zmin+k0*h
        if smoothing != none && term > 0
          dFdz = 0.0
          if smoothing == trapezoidal
            dFdz = (ftilde - Fm[1]) / (ztilde - density.zm[k0+1])
          elseif smoothing == midpoint
            dFdz = 2 * (ftilde - Fm[1]) / density.h
          end
          p += 0.5 * term * term * density.ϑ[k0+1] / dFdz
        end
        p += sum((Fm[2:end-1] .- strike) * density.h .* density.ϑ[k0+2:end-1])
      end
    end
  else
  end

  return p
end


end # module
