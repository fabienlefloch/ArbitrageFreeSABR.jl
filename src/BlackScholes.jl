using Statistics

export blackScholesFormula
function blackScholesFormula(isCall::Bool, strike::T, spot::T, variance::T, driftDf::T, discountDf::T) where {T}

    sign = 1.0
    if !isCall
        sign = -1.0
    end
    forward = spot / driftDf
    if variance < eps()
        price = discountDf * max(sign * (forward - strike), 0.0)
        return price
    elseif spot < eps()
        if isCall
            return 0.0
        else
            return discountDf * strike
        end
    elseif strike < eps()
        if isCall
            return discountDf * forward
        else
            return 0.0
        end
    else
        sqrtVar = sqrt(variance)
        d1 = 1.0 / sqrtVar * log(forward / strike) + 0.5 * sqrtVar
        d2 = d1 - sqrtVar
        normal = Normal()
        nd1 = cdf(normal, sign * d1)
        nd2 = cdf(normal, sign * d2)
        price = sign * discountDf * (forward * nd1 - strike * nd2)
        return price
    end
end
