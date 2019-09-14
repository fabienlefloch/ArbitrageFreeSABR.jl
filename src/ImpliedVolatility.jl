using Distributions
using SpecialFunctions

export impliedVolatilityLiSORTS

function impliedVolatilitySqrtTimeRationalGuess(x::T, c::T) where {T}
    m00 = -0.00006103098165
    n00 = 1.0
    m01 = 5.33967643357688
    n01 = 22.96302109010794
    m10 = -0.40661990365427
    n10 = -0.48466536361620
    m02 = 3.25023425332360
    n02 = -0.77268824532468
    m11 = -36.19405221599028
    n11 = -1.34102279982050
    m20 = 0.08975394404851
    n20 = 0.43027619553168
    m03 = 83.84593224417796
    n03 = -5.70531500645109
    m12 = 41.21772632732834
    n12 = 2.45782574294244
    m21 = 3.83815885394565
    n21 = -0.04763802358853
    m30 = -0.21619763215668
    n30 = -0.03326944290044
    num = m00 + m10 * x + m20 * x^2 + m30 * x^3 + c * (m01 + m11 * x + m21 * x^2) + c^2 * (m02 + m12 * x) + c^3 * m03
    den = n00 + n10 * x + n20 * x^2 + n30 * x^3 + c * (n01 + n11 * x + n21 * x^2) + c^2 * (n02 + n12 * x) + c^3 * n03
    return num / den
end

function impliedVolatilityLiRationalGuess(
    price::T,
	isCall::Bool,
	strike::T,
    forward::T,
    timeToExpiry::T,
    df::T
) where {T}
    c = price / forward / df
    ex = forward / strike
    x0 = log(ex)
    x = x0
    if !isCall
        if x0 <= 0
            c = c + 1 - 1.0 / ex # put call parity
        else
			#duality + put call parity
            c = ex * c
            x = -x0
            ex = 1.0 / ex
        end
    else
        if x0 > 0
			# use c(-x0, v)
            c = (forward * (c - 1) + strike) / strike #in out duality
            x = -x0
            ex = 1.0 / ex
        end
    end
    if c >= 1 / ex || c >= 1
        error("price ", c, " higher than intrinsic value ", 1 / ex)
    elseif c <= 0
        error("price ", c, " negative")
    end
    sqrtte = sqrt(timeToExpiry)
    return impliedVolatilitySqrtTimeRationalGuess(x, c) / sqrtte
end

function impliedVolatilityLiSORTS(
    price::T,
    isCall::Bool,
    strike::T,
    forward::T,
    timeToExpiry::T,
    df::T,
    impliedVolGuess::T,
    tolerance::T,
    maxIterations::Int
) where {T}
    c = price / forward / df
    ex = forward / strike
    x0 = log(ex)
    x = x0
    if !isCall
        if x0 <= 0
            c = c + 1 - 1.0 / ex # put call parity
        else
			#duality + put call parity
            c = ex * c
            x = -x0
            ex = 1.0 / ex
        end
    else
        if x0 > 0 			# use c(-x0, v)
            c = (forward * c - forward + strike) / strike # in out duality
            x = -x0
            ex = 1.0 / ex
        end
    end
    if c > 1 / ex
        return 0
    end
    sqrtte = sqrt(timeToExpiry)
    v0 = 0.0
    if impliedVolGuess == 0.0
        v0 = impliedVolatilitySqrtTimeRationalGuess(x, c)
    else
        v0 = impliedVolGuess * sqrtte
    end
    return computeLiSORTS(c, x, ex, v0, sqrtte, tolerance, maxIterations)
end

const OneOverSqrtTwo = sqrt(0.5)

function computeLiSORTS(c::T, x::T, ex::T, v0::T, sqrttte::T, tolerance::T, maxIterations::Int) where {T}
    v = v0
    if maxIterations <= 0
        return v
    end
    exinvsqrt = 1.0 / sqrt(ex)
    h = x / v
    t = v / 2

    Np = erfcx(-OneOverSqrtTwo * (h + t))
    Nm = erfcx(-OneOverSqrtTwo * (h - t))
    norm = exinvsqrt * 0.5 * exp(-0.5 * (h^2 + t^2))
    cEstimate = norm * (Np - Nm)
    if cEstimate < 0
        cEstimate = 0.0
    end

    omega = 1.0
    iterations = 0
    normal = Normal()
    while abs(c - cEstimate) > tolerance && iterations < maxIterations
        phi = (v^2 + 2 * x) / (v^2 - 2 * x)
        alpha = (1 + omega) / (1 + phi)
        F = c + norm * (Nm + omega * Np)
        Nm1 = quantile(normal, F / (1 + omega))
        G = Nm1 + sqrt(Nm1^2 - 2 * x)
        v = alpha * G + (1 - alpha) * v
        h = x / v
        t = v / 2
        Np = erfcx(-OneOverSqrtTwo * (h + t))
        Nm = erfcx(-OneOverSqrtTwo * (h - t))
        norm = exinvsqrt * 0.5 * exp(-0.5 * (h^2 + t^2))
        cEstimate = norm * (Np - Nm)
        if cEstimate < 0
            cEstimate = 0
        end
        iterations += 1
    end
    sigma = v / sqrttte
    return sigma
end
