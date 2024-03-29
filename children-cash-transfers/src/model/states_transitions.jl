"""
    # TODO: describe income process
    fε(kε′ | kε) = πW * (kε==kε′) + 0.5(1-πW)*(kϵ′=max(1,kε-1)) + 0.5(1-πW)*(kϵ′=min(kε+1,Kε))

"""

# TODO: write a function to evalute these probabilities only over states in which f(kε′|kε) is positive
# e.g. function that returns (max(1,kϵ-1),kϵ,min(p.Kϵ,kϵ+1)) with (0.5(1-πW),πW,0.5(1-πW))
# then use this in calc_vj

function fε(kε, Kε, πW)
    if kε == 1
        return (1, 2), (0.5 * (1 - πW), πW + 0.5 * (1 - πW))
    elseif kε == Kε
        return (Kε - 1, Kε), (0.5 * (1 - πW) + πW, 0.5 * (1 - πW))
    else
        return (kε - 1, kε, kε + 1), (0.5 * (1 - πW), πW, 0.5 * (1 - πW))
    end
end
# <- missing end statement.

# note here that if kϵ=1, the function will return: (1,2),(πW,0.5*(1-πW)) which doesn't sum to 1.
# given how calc_vj is written, it would be OK to return (1,1,2),(0.5*(1-πW),πW,0.5*(1-πW))
# returning (1,2), (πW + 0.5*(1-πW),0.5*(1-πW)) would also work.
# I think it would be better to return the latter, which will make more sense when calculating the likelihood of transitions later on
# I think the function should also return the probabilties as a tuple (immutable with no allocations) instead of as an array (mutable with allocations)