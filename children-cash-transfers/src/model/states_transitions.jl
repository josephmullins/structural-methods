"""
    # TODO: describe income process
    fε(kε′ | kε) = πW * (kε==kε′) + 0.5(1-πW)*(kϵ′=max(1,kε-1)) + 0.5(1-πW)*(kϵ′=min(kε+1,Kε))

"""

function fε(kε, Kε, πW)
    if kε == 1
        return (1, 2), (πW + 0.5 * (1 - πW), 0.5 * (1 - πW))
    elseif kε == Kε
        return (Kε - 1, Kε), (0.5 * (1 - πW), 0.5 * (1 - πW) + πW)
    else
        return (kε - 1, kε, kε + 1), (0.5 * (1 - πW), πW, 0.5 * (1 - πW))
    end
end
