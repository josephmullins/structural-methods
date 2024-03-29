include("../src/Transfers.jl")


function pars(Kτ::Int64,Kε::Int64)
    β = 0.98
    αC = 1.
    αθ = 0.1 * ones(Kτ)
    αH = zeros(Kτ)
    αA = zeros(Kτ)
    αS = zeros(Kτ)
    αR = zeros(2)
    σ = fill(2.,2)


    αl=1 #added this in for the utility function

    σ_W = 2.

    πW = 0.8

    σε = 2.
    ση = 2.

    βΓx = zeros(2)
    βΓτ = zeros(2)
    βw = [LinRange(6,7.5,Kτ)';zeros(2,Kτ)]
    
    βτ = zeros(23,Kτ-1) #<- to be determined
    πε = ones(Kε,Kτ) / Kε
    
    εgrid = LinRange(-1,1,Kε) 
    return (;Kτ,Kε,β,αC,αθ,αH,αA,αS,αR,σ,σ_W,σε,βΓx,βΓτ,βw,βτ,πε,εgrid,αl,πW,ση)
end

struct model_data #
    T::Int64 #<- length of problem

    y0::Int64 #<- year to begin problem
    age0::Int64 # <- mother's age at start of problem
    SOI::Vector{Int64} #<- state SOI in each year
    num_kids::Vector{Int64} #<- number of kids in household that are between age 0 and 17
    TotKids::Int64 #<- indicares the total number of children that the mother will have over the available panel
    age_kid::Matrix{Int64} #< age_kid[f,t] is the age of child f at time t. Will be negative if child not born yet.
    cpi::Vector{Float64} #<- cpi

    R::Vector{Int64} #<- indicates if work requirement in time t
    Kω::Int64 #<- indicates length of time limit once introduced
    TL::Vector{Bool} #<- indicating that time limit is effective
end

function test_model()
    T = 12
    y0 = 1990
    age0 = 23
    SOI = fill(12,T)
    num_kids = fill(1,T)
    TotKids = 1
    age_kid = [a for j in 1:1, a in 1:T] 
    cpi = fill(1.,T)

    R = fill(0,T)
    Kω = 1
    TL = fill(false,T)
    return model_data(T,y0,age0,SOI,num_kids,TotKids,age_kid,cpi,R,Kω,TL)
end

include("model/choices.jl")
include("model/utility.jl")
include("model/states_transitions.jl")
include("model/solve.jl")