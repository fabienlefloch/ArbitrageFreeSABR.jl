using Statistics

export bachelierFormula

"Price of a vanilla option under the Bachelier/Normal model"
function bachelierFormula(isCall::Bool, strike::T, forward::T, variance::T, discountDf::T) where {T}
    sign = 1.0
    if !isCall
        sign = -1.0
    end
    if variance == 0
        return Max(sign * (forward - strike), 0) * discountDf
    end
    sqrtvar = sqrt(variance)
    d = sign * (forward - strike) / sqrtvar
    if forward == strike
        return sqrtvar / sqrt(2 * pi)
    end
    normal = Normal()
    Nd = cdf(normal, d)
    nd = pdf(normal, d)
    return discountDf * (sign * (forward - strike) * Nd + sqrtvar * nd)
end
