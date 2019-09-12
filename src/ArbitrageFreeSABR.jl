module ArbitrageFreeSABR
using LinearAlgebra

struct SABRMaturity{T <: Real} 
    alpha::T
    beta::T
    rho::T
    nu::T
    forward::T
    time::T
end

struct TransformedDensity{T <: Real}
    maturity::SABRMaturity{T}
    zmin::T
    zmax::T
    h::T
    zm::Vector{T}
    P::Vector{T}
    PL::T
    PR::T
end

export SABRMaturity, makeTransformedDensityLawsonSwayne,priceTransformedDensity

function yOfStrike(m::SABRMaturity, strike) 
    return (strike.^(1-m.beta) .- m.forward^(1-m.beta))/(1-m.beta)
  end

  function makeForward(m::SABRMaturity, z)
    y = m.alpha/m.nu*(sinh.(m.nu*z)+m.rho*(cosh.(m.nu*z) .- 1.0))
    return (m.forward^(1-m.beta).+(1-m.beta)*y).^(1/(1-m.beta))
  end

  function solveStep(P::Vector{G}, PL::G, PR::G, h::G, Fm::Vector{G}, Cm::Vector{G}, Em::Vector{G}, dt::G) where {G}
    frac = dt/(2*h); M = length(P)
    B=zeros(M)
    A=zeros(M-1)
    C=zeros(M-1)
    @. B[2:M-1] = 1.0 + frac*(Cm[2:M-1] * Em[2:M-1]*(1.0 / (Fm[3:M]-Fm[2:M-1])+ 1.0 /(Fm[2:M-1]-Fm[1:M-2])))
    @. C[2:M-1] = -frac* Cm[3:M]*Em[3:M]/(Fm[3:M]-Fm[2:M-1])
    @. A[1:M-2] = -frac* Cm[1:M-2]*Em[1:M-2]/(Fm[2:M-1]-Fm[1:M-2])
    B[1] = Cm[1]/(Fm[2]-Fm[1])*Em[1]
    C[1] = Cm[2]/(Fm[2]-Fm[1])*Em[2]
    B[M] = Cm[M]/(Fm[M]-Fm[M-1])*Em[M]
    A[M-1] = Cm[M-1]/(Fm[M]-Fm[M-1])*Em[M-1]  
    tri = Tridiagonal(A,B,C)
    p1 = P[1]; pm = P[M]; P[1] = 0; P[M] = 0
    P0 = tri\P
    PL0 = PL + dt*Cm[2]/(Fm[2]-Fm[1])*Em[2]*P0[2]
    PR0 = PR + dt*Cm[M-1]/(Fm[M]-Fm[M-1])*Em[M-1]*P0[M-1]
    P[1] = p1; P[M] = pm
    return P0,PL0,PR0
  end

  function computeBoundaries(alpha, beta, nu, rho, forward, T, nd)
    zmin = -nd*sqrt(T); zmax = -zmin
    if (beta < 1) 
      ybar = -forward^(1-beta)/(1-beta);
      zbar = -1/nu*log((sqrt(1-rho^2+(rho+nu*ybar/alpha)^2)-rho-nu*ybar/alpha)/(1-rho))
      if (zbar > zmin)
        zmin = zbar
      end    
    end
    return zmin, zmax
  end
  
  function yofz(m::SABRMaturity{G}, zm) where {G}
    return m.alpha/m.nu*(sinh.(m.nu*zm)+m.rho*(cosh.(m.nu*zm) .- 1))
  end
  function fofy(m::SABRMaturity{G}, ym) where {G}
    return (m.forward^(1-m.beta) .+ (1-m.beta)*ym).^(1/(1-m.beta))
  end
  function cofy(m::SABRMaturity{G}, ym, Fm) where {G}
    return sqrt.(m.alpha^2 .+ 2*m.rho*m.alpha*m.nu*ym + m.nu^2*ym.^2).*Fm.^(m.beta)
  end
  function gammaOfF(m::SABRMaturity{G}, Fm, j0::Int) where {G}
    Gamma = (Fm.^m.beta .- m.forward^m.beta)./(Fm .- m.forward)
    Gamma[j0+1] = m.beta/m.forward.^(1-m.beta)
    return Gamma
  end
  
  
function makeTransformedDensityLawsonSwayne(maturity::SABRMaturity{G}, N::Int, timesteps::Int, nd::Real) where {G}
    alpha = maturity.alpha
    beta = maturity.beta
    nu = maturity.nu
    rho = maturity.rho
    forward = maturity.forward
    T = maturity.time
    # initialize h,F and Q
    zmin, zmax = computeBoundaries(alpha, beta, nu, rho, forward, T, nd)
    J = N-2; h0 = (zmax-zmin)/J; j0 = ceil((0-zmin)/h0)
    h = (0-zmin)/(j0-0.5)
    z = collect(0:(J+1))*h .+ zmin; zmax = z[J+2]; zm = z .- 0.5*h
    ym = yofz(maturity, zm)
    ym[1] = ym[2] #avoid negative
    Fm = fofy(maturity, ym)
    Cm = cofy(maturity, ym, Fm)
    Gammam = gammaOfF(maturity, Fm, Int(j0))
    ymax = yofz(maturity, zmax); ymin = yofz(maturity, zmin)
    Fmax = fofy(maturity, ymax); Fmin = fofy(maturity, ymin)
    Fm[1] = 2*Fmin-Fm[2]; Fm[J+2]= 2*Fmax - Fm[J+1]
    Cm[1] = Cm[2]; Cm[J+2] = Cm[J+1]
    dt = T/timesteps
    #Ptotal = sum(h*P[2:J+1])+PL+PR
    #Ftotal = Fm[2:J+1]*P[2:J+1]*h+Fmin*PL+Fmax*PR
    b = 1 - sqrt(2)/2 #Lawson Swayne param
    dt1 = dt*b; dt2 = dt*(1-2*b);  Em = ones(J+2)
    Emdt1 = exp.(rho*nu*alpha*Gammam*dt1); Emdt1[1]=Emdt1[2]; Emdt1[J+2]=Emdt1[J+1]
    Emdt2= exp.(rho*nu*alpha*Gammam*dt2); Emdt2[1]=Emdt2[2]; Emdt2[J+2]=Emdt2[J+1]
    PL = 0.0; PR = 0.0; P = zeros(J+2); P[Int(j0)+1]=1.0/h
    for t = 1:timesteps
        @. Em *= Emdt1
        P1,PL1,PR1 = solveStep(P,PL,PR,h,Fm,Cm,Em,dt1)
        @. Em *= Emdt1
        P2,PL2,PR2 = solveStep(P1,PL1,PR1,h,Fm,Cm,Em,dt1)
        @. P=(sqrt(2)+1)*P2-sqrt(2)*P1
        PL=(sqrt(2)+1)*PL2-sqrt(2)*PL1
        PR=(sqrt(2)+1)*PR2-sqrt(2)*PR1    
        @. Em *= Emdt2        
      #Ptotal = sum(h*P[2:J+1])+PL+PR
      #Ftotal = Fm[2:J+1]*P[2:J+1]+h+Fmin*PL+Fmax*PR
    end
    density = TransformedDensity(maturity,zmin,zmax,h,zm,P,PL,PR)
    return density
  end

@enum DensitySmoothing none=1 trapezoidal=2 midpoint=3
  function priceTransformedDensity(density::TransformedDensity{G},isCall::Bool, strike::G, smoothing::DensitySmoothing) where {G}
    maturity = density.maturity
    alpha = maturity.alpha
    beta = maturity.beta
    nu = maturity.nu
    rho = maturity.rho
    forward = maturity.forward
    ystrike = yOfStrike(maturity,strike)
    zstrike = -1/nu*log((sqrt(1-rho^2+(rho+nu*ystrike/alpha)^2)-rho-nu*ystrike/alpha)/(1-rho))
    if isCall
    if (zstrike < density.zmin)
      p = forward-strike
    else 
      if (zstrike > density.zmax)
        p = 0
      else
        Fmax = makeForward(maturity, density.zmax)
        p = (Fmax - strike) * density.PR
        k0 = Int(ceil((zstrike-density.zmin)/density.h))
        ztilde = density.zmin + k0*density.h
        ftilde = makeForward(maturity, ztilde)
        term = ftilde - strike
        Fm = makeForward(maturity, density.zm[k0+1:end]) #zm[k0+1] = zmin+k0*h
        if smoothing != none && term > 0
          dFdz = 0.0
          if smoothing == trapezoidal 
            dFdz = (ftilde-Fm[1]) / (ztilde-density.zm[k0+1])
          elseif smoothing == midpoint
            dFdz = 2*(ftilde-Fm[1]) / density.h
          end
          p += 0.5 * term * term * density.P[k0+1]/dFdz
        end
         p += sum((Fm[2:end-1] .- strike) * density.h .* density.P[k0+2:end-1])
      end
    end
else
  end

return p
end


end # module
