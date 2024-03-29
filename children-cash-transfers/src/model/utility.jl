"""
    utility(j,p,md::model_data,kτ,kε,kω,t)

Evaluate the utility from making choice j given the parameters p, model data md, and the state, x.

In this version of the model, utility takes the form:

`` uₘ(j,kτ,kε,kω,t) = α̃₁,ₖ × log(Yₘ(j,kτ,kε,kω,t)) + α̃₂,ₖ × log(112 - 30Hⱼ) - αₛ,ₖ Sⱼ - (αₐ,ₖ + Rₘ,ₜ×αᵣ,₁(1-Hⱼ)) Aⱼ - (αₕ,ₖ + Rₘ,ₜ×αᵣ,₂Aⱼ) Hⱼ ``

where:
- `j` indexes the discrete choice and the function [`j_inv`](@ref) returns the choices ``Sⱼ, Aⱼ, Hⱼ`` that are associated with choice ``j``
- `kτ` indexes type
- `kε` indexes the value of the wage shock
- `kω` indicates how many years of welfare have been used so far
- ``md::model_data`` holds the fixed and observable components of the state space (state of residence, family composition)
- Following the paper, we have the below formula for `α̃₁` and `α̃₂`:
    - ` α̃₁,ₖ = 1 + p.αθ[kτ] × sum([Γₐ(p.βΓx,md.age_kid[f,t]) for f in 1:md.TotKids]) `
    - ` α̃₂,ₖ = p.αl + p.αθ[kτ] × sum([Γₐ(p.βΓτ,md.age_kid[f,t]) for f in 1:md.TotKids]) `
where we call the function [`Γₐ`](@ref) to approximate the terms ``Γxₐ`` and ``Γτₐ`` from the paper
- The function [`budget`](@ref) returns net income ``Yₘ(j,kτ,kε,kω,t)`` given the state and observed choice `j`.
- The vectors ``αₛ,αₐ,αₕ,αᵣ`` are stored in `p` as `αS`,`αA`,`αH`,`αR`
- ``Rₘ,ₜ`` indicates whether a work requirement is in effect (stored in `md`)
"""
function utility(S,A,H,p,md::model_data,kτ,kε,kω,t)

    sum_term1 = 0.0
    sum_term2 = 0.0

    for f in 1:md.TotKids
        age_kid_f_t = md.age_kid[f, t]
        sum_term1 += Γₐ(p.βΓx, age_kid_f_t)
        sum_term2 += Γₐ(p.βΓτ, age_kid_f_t)
    end

    a = (1 + p.αθ[kτ] * sum_term1) * log(budget(S, A, H, p, md, kτ, kε, kω, t))
    b = (p.αl + p.αθ[kτ] * sum_term2) * log(112 - 30 * H)
    c = p.αS[kτ] * S + ( p.αA[kτ] + md.R[t] * p.αR[1] * (1 - H) ) * A #<- missing brackets? also αA used twice. what about αS?
    d = (p.αH[kτ] + md.R[t] * p.αR[2] * A) * H #<- αR used twice. want to use αH?

    return a + b - c - d

end

"""
    Γₐ(β,a)

Evaluate the quadratic approximation for the expressions Γxₐ and Γτₐ from the paper.
"""
Γₐ(β,a) = (a>=0 && a<=17) * exp(β[1]*a + β[2]a)

"""
    budget(S,A,H,p,md::model_data,kτ,kε,kω,t)
    
    This function determines whether eligible for food stamps and, if working, calculates the wage. Then it calls the function `budget` found in the module `Transfers.jl`
"""
function budget(S,A,H,p,md::model_data,kτ,kε,kω,t)
    year = md.y0 + t - 1 
    if H==1
        W = exp(logwage(p,md::model_data,kτ,kε,t))
    else
        W = 0.
    end

    if md.TL[t] && kω==md.Kω #<- ineligible due to time limit in effect
        return Transfers.budget(W,0.,md.SOI[t],year,md.num_kids[t],md.cpi[t],S)
    else
        return Transfers.budget(W,0.,md.SOI[t],year,md.num_kids[t],md.cpi[t],S+A)
    end
end

"""
    logwage(p,md::model_data,kτ,kε,t)

Returns log wages given observable data `md::model_data` and the state: `kτ,kε` according to the formula:

`` log(Wₘ,ₖ,ₜ) = βw[kτ,1] + βw[kτ,2]×Ageₘ,ₜ + p.βw[kτ,3]×Age²ₘ,ₜ + ση * η[kε] ``

where 
- ``βw₀,βw₁,βw₂``` are stored sequentially in the vector `p.βw`
- ``Ageₘ,ₜ`` is stored at `md.age0[t]`
-  ``η`` is stored in `p.ηgrid` and ση is also stored in p
"""
function logwage(p,md::model_data,kτ,kε,t)
    age = md.age0 + t
    return p.βw[1,kτ] + p.βw[2,kτ]*age + p.βw[3,kτ]*age^2 + p.ση*p.εgrid[kε]

end
# note that age0 in model data is an integer, not a vector, so age must be calculated as above
# note the change of indexing for βw