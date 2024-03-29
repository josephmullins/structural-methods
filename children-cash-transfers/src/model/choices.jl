# - here we define choices and their mapping to a single index
# - we also define the structure of the nested logit

# choices are: 
# S ∈(0,1) food stamps
# A ∈(0,1) AFDC/TANF (A=1 only if S=1)
# P = S + A indexes program participation (sometimes useful)
# H ∈(0,1) labor supply (restricted 30 hours)
# F ∈(0,1) formal care (F=1 only if H=1)

j_idx(S,A,H) = (S + A)*2 + H + 1

J = 6 #<- 6 choices, total
function j_inv(j)
    p = fld(j-1,3)
    H::Int64 = mod(j-1,3)
    S::Int64 = p>0
    A::Int64 = p==2
    return S,A,p,H
end

# this is the nested logit structure for this model
function get_nests()

    B₁ = [[1,2],[3,4],[5,6]]
    C₁ = [[1,],[2,],[3,],[4,],[5,],[6,]]
    
    B₂ = [[1,2,3]]
    C₂ = [[1,2],[3,4],[5,6]]

    B = (B₁,B₂)
    C = (C₁,C₂)
    return (;B,C)
end