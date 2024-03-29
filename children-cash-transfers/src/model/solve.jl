"""
    solve!(logP,V,vj,p,md::model_data)

A function to solve the model by backward induction and store the choice probabilities in logP, a 3D array.
    - `V` is a 2D array that stores in one column the value of arriving next period in state `k` while the other column gets filled with values for the current period. 
    In the next iteration, the columns switch.
    - `md` contains the information necessary to determine `T`, the finite horizon of the problem, It is also used as an argument in `iterate!`, `calc_vj`, and `utility`

This version of the model contains three state variables: `kε`, which indexes the wage shock, `kω` which tracks welfare usage once time limits are applied 
(which is given by `md.TL[t]` in period `t`), and `kτ` which indexes type. The tuple `(p.Kε,md.Kω,p.Kτ)` is used an argument for `LinearIndices` and 
`CartesianIndices` for convenient indexation of the state.

"""
# the function will eventually go here

"""
    iterate!(logP,V,vj,p,md,t,tnow,state_idx)

A function to calculate the choice probabilities `logP[j,k]` for each choice `j` and state `k` where:
- `logP[j,k]` is a 2D array of choice probabilties
- `V` is a 2D array of values where `V[:,3-tnow]` holds the value of arriving next period in each state
- `vj` is a buffer for storing the value of each choice `j` when filling in `logP[:,k]`
- `p` is a named tuple containing the parameters as well as `B` and `C`, the arrays that summarize the nested logit structure
- `state_idx` is a named tuple containing:
    - `K`: the size of the state space    (;K,k_idx,k_inv) = state_idx
    - `k_idx`: a Linear index that converts the tuple of individual states to a single index
    - `k_inv`: a cartesian index that inverts the single index back to a tuple

This function calls [`calc_vj`](@ref) to evaluate choice-specific values and [`nested_logit`](@ref) to calculate choice probabilities given vj.
Note that the `logP` and `V` are buffers that may be re-used to solve different sized models, and so `state_idx.K` should always be used for iteration, not 
`axes(V,1)` or `axes(logP,2)`, each of which may be larger than `K`.

"""
# function eventually here

"""
    calc_vj(j,V,md::model_data,state,pars,t)
A function to calculate the choice-specific value ``vⱼ`` given by:
    `` vⱼ(k) = uⱼ(k) + β * ∑ fⱼ(k'|k) × V(k') `` where
- `j` indexes the discrete choice
- `state` holds the individual state variables `kω`,`kε`, `kτ`
- `state` also holds `k_idx` a linear index for the triple above
- `V[k′]` holds the value of arriving tomorrow in state `k′`
- ``fⱼ(k'|k)`` is the probability of being in state k' tomorrow given choice `j` and current state `k`.

Note that the implementation of the sum and `fⱼ` will be model specific to maximize efficiency of iteration.
This version of the model contains three state variables: `kε`, which indexes the wage shock, `kω` which tracks welfare usage once time limits are applied (which is given by `md.TL[t]` in period `t`), and `kτ` which indexes type. Of these, `kω` evolves deterministically as `` kω' = max(Kω,kω + A(j)*TL) `` where `TL` is a binary indicating whether time limits apply, `A(j)` indicates whether welfare use is associated with choice `j` and `Kω` is 1 + the state adopted time limit. `kϵ` evolves stochastically as described in `states_transitions.jl`.

"""
# function eventually here

""" 
    get_model(p,md)
    
    A function that returns an appropriately sized triple (V,vj,logP)
    where:
    - V stores inclusive values
    - vj stores choice-specific values
    - logP stores choice probabilities (J x K x T)
"""
function get_model(p,md)
    T = md.T
    J = 6
    K = p.Kε * p.Kτ * md.Kω
    logP = zeros(J,K,T)
    V = zeros(K,2)
    vj = zeros(J) 
    return (;logP,vj,V)
end


# a function that takes initial values vj and calculates nested logit probabilities
# Some notes:
# - log-choice probabilities are written to the vector logP
# - vj is also used as a buffer to store inclusive values at each layer
# - **note**: this algorithm works as long as the partitions are written with nodes in ascending order. For example:
#   - e.g. the partition B[l] = [[1,2],[3,4,5]] is ok because the value for nest [1,2] will be written to vj[1] and vj[2] for nest [3,4,5]
#   - e.g. the partition B[l] = [[2,3],[1,4,5]] will be incorrect because the nest [2,3] will write to vj[1], which still needs to be used for the nest [1,4,5]
function nested_logit(logP,vj;B,C,σ)
    fill!(logP,0.)
    for l ∈ eachindex(B)
        Cₗ = C[l] #<- each Cₗ is a K(l)-vector of vectors containing the choices that are descendents of each node k.
        for k ∈ eachindex(B[l])
            bₖ = B[l][k] #<- indicates which nodes are in the kth partition of layer l
            vmax = -Inf
            # find the maximum
            for j ∈ bₖ
                vj[j]>vmax ? vmax=vj[j] : nothing
            end
            norm = 0.
            for j ∈ bₖ
                norm += exp((vj[j] - vmax) / σ[l])
            end
            norm = log(norm)
            # then: 
            for jₗ ∈ bₖ
                for j ∈ Cₗ[jₗ]
                    logP[j] += (vj[jₗ] - vmax) / σ[l] - norm
                end
            end
            vj[k] = vmax + σ[l] * norm #<- write the inclusive value of node k to vj
        end
    end
    return vj[1]
end

function plain_logit(logP,vj;B,σ)
    norm = 0.
    vmax = -Inf
    for j ∈ B
        vj[j]>vmax ? vmax=vj[j] : nothing
    end
    for j ∈ B
        norm += exp((vj[j] - vmax) / σ)
    end
    norm = log(norm)
    for j ∈ B
        logP[j] = (vj[j] - vmax) / σ - norm
    end
    vj[1] = vmax + σ * norm
end