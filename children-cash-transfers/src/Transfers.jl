module Transfers 
export budget, SNAP, TANF
using CSV,DataFrames
using DelimitedFiles
# some notes on this budget function:
# - compute everything in monthly terms
# - assumes that earnings and nonlabor income have been deflated by cpi
# - assumes that the tax function takes annual values
# - assumes the state is SOI code
# - assumes nonlabor income is not reported in benefit calculations
function budget(E,N,state,year,num_kids,cpi,p)
    E *= cpi #<- re-inflate income variables
    N *= cpi

    tax = 0.
    if E>0
        tax = TAX(E*12,state,year,num_kids)
        tax /= 12
        # in case you see this later and think there should be a dtax /= 12...
        # remember: d[f(α*x)/α]/dx = f'(α*x)
    end
    snap = 0.
    tanf = 0.
    if p>0
        if p>1 && num_kids>0
            tanf = TANF(E,0.,state,year,num_kids)
        end
        snap = SNAP(E,0.,year,num_kids,tanf)
    end
    y = (E + N + snap + tanf - tax)/cpi
    return y
end

# - alternative version assumes everything is in nominal terms
function budget(E,N,state,year,num_kids,p)
    return budget(E,N,state,year,num_kids,1.,p)
end

# - this version is for when TANF entitlements run out (i.e. food stamps only)
#  - why no cpi here??
function budget_inel(E,N,state,year,num_kids,p)
    tax = 0
    if E>0
        tax = TAX(E*12,state,year,num_kids)[1]
        tax /= 12
    end
    snap = 0
    tanf = 0
    if p>0
        snap = SNAP(E,0.,year,num_kids,0.)
    end
    y = E + N + snap + tanf - tax
    return y
end


# ------------------ SNAP ---------------------------------------------------------------- #
#snap_rules_csv = CSV.read("data/WelfareRules/SNAPRules.csv",DataFrame)
fp = joinpath(@__DIR__, "..", "data", "WelfareRules", "SNAPRules.csv")
snap_rules_csv = CSV.read(fp,DataFrame)

const PG::Vector{Float64} = convert(Array{Float64,1},snap_rules_csv.PG) 
const MA::Vector{Float64} = convert(Array{Float64,1},snap_rules_csv.MA)

function SNAP(E,N,year,numchild,tanf)
    ii = (year-1979)*5 + numchild+1
    return SNAP(E,N,PG[ii],tanf,0.,MA[ii])
end

function SNAP(E,N,PG,WB,SD,M)
    # ARGUMENTS
    # E - earnings
    # N - non-labor income
    # PG - poverty guideline
    # WB - welfare benefit received
    # SD - shelter deduction
    # M - maximum benefit
    #  --------------------------------------- #
    # - 1. Eligibility test
    EG = (E+N)<1.3*PG
    NE = max(0.8*E + N + WB - SD - 134,0.)
    EN = NE<PG
    # - 2. payment
    snap = max(M-0.3*NE,0)
    return EG*EN*snap
end



# -------------------- TANF ------------------------------------------ #
#tanf_rules_csv = CSV.read("data/WelfareRules/WelfareRules.csv",DataFrame)
fp = joinpath(@__DIR__, "..", "data", "WelfareRules", "WelfareRules.csv")
tanf_rules_csv = CSV.read(fp,DataFrame)

const NS::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.NS)
const GR::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.GR)
const NR::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.NR)
const EDD::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.EDD)
const ERD::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.ERD)
const BDD::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.BDD)
const BRD::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.BRD)
const BS::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.BS)
const MB::Vector{Float64} = convert(Array{Float64,1},tanf_rules_csv.MB)

function TANF(E,N::Float64,SOI::Int64,year::Int64,numchild::Int64)
    ii = (year-1975)*51*4 + (SOI-1)*4 + numchild
    return TANF(E,N,NS[ii],GR[ii],NR[ii],EDD[ii],ERD[ii],BDD[ii],BRD[ii],BS[ii],MB[ii])
end

function TANF(E,N,NS,rg,rn,dd0,rd0,dd1,rd1,PS,M)
    # ARGUMENTS:
    # E - earnings
    # N - non-labor income
    # NS - need standard
    # rg - gross adjustment
    # rn - net adjustment
    # dd0 - dollar disregard on eligibility
    # rd0 - rate of disregard on eligibility
    # PS - payment standard
    # M - max payment
    # dd1 - dollar disregard on payment (same as eligibility)
    # rd1 - rate of disregard on payment (same as eligibility)
    # ------------------------------------------ #
    # - 1. Eligibility tests
    if rg>0
        EG = (E+N)<rg*NS # gross test
    else
        EG = 1
    end

    net_income = max(E-dd0,0.)*(1-rd0) + N
    if rn>0
        EN = net_income<rn*NS # net test
    else
        EN = 1
    end

    # - 2. payment
    net_income = max(E-dd1,0.)*(1-rd1) + N #<- also an error here? if dd1>0, there is a flat region
    #println(1-rd1)
    tanf = max(min(M,PS-net_income),0)
    return EG*EN*tanf
end


# --------------------- TAXES ------------------------------------ #
#const TDAT = readdlm("data/WelfareRules/TaxsimVector.csv",',')[:] # great
fp = joinpath(@__DIR__, "..", "data", "WelfareRules", "TaxsimVector.csv")
const TDAT = readdlm(fp,',')[:] # great

# -- TDAT is data ordered as follows:  year,state,depx,earnings
# years: 1970-1976 (federal only), 1977-1978
# states: 0-51 (SOI codes, not FIPS!)
# depx - 0,1,2,3,4
# earnings in $500 grids, up to 99,000


# this function needs to take annual thousands of dollars!!!
function TAX(E,state,year,num_kids)
    # step 1: locate the bracket in which the mother's earnings lies (this can be mechanical)
    Ns = 52
    Ne = 100
    Nd = 5
    N_ = (1976-1970+1)*Nd*Ne
    # first, get indices
    if E>99000
        ie = 98
    else
        ie = Int(fld(E,1000))
    end
    if year<=1976
        iy = year-1970
        ii = iy*Ne*Nd + num_kids*Ne + ie+1
    else
        iy = year-1977
        ii = N_+iy*Ne*Ns*Nd + state*Ne*Nd + num_kids*Ne + ie+1
    end
    T0 = TDAT[ii]
    T1 = TDAT[ii+1]
    E0 = ie*1000
    tax = T0 + (T1-T0)/1000*(E-E0)
    return tax
end
end