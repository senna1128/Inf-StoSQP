
## Load packages
using NLPModels
using JuMP
using LinearOperators
using OptimizationProblems
using MathProgBase
using ForwardDiff
using CUTEst
using NLPModelsJuMP
using LinearAlgebra
using Distributed
using Ipopt
using DataFrames
using PyPlot
using MATLAB
using Glob
using DelimitedFiles
using Random
using Distributions
using NLPModelsIpopt

cd("/.../StoSQPFrame")
Prob = readdlm(string(pwd(),"/../Parameter/problems.txt"))

module Parameter
    # Parameters of StoSQP with augmented Lagrangian
    struct StoSQP
        verbose                            # Do we create dump dir?
        # stopping parameters
        MaxIter::Int                       # Maximum Iteration
        # fixed parameters
        Rep::Int                           # Number of Independent runs
        tau::Int                           # Number of iterations for inexact solver
        DBeta::Array{Float64}              # Decay stepsize 2/(t^(c_2)) with 0.5<c_2<=1
        # test parameters
        Sigma::Array{Float64}              # variance of gradient
    end
end


using Main.Parameter
include("StoSQPMain.jl")


#######################################
#########  run main file    ###########
#######################################
function main()
    ## run StoSQP framework
    include("../Parameter/Param.jl")
    StoSQPR = StoSQPMain(StoSQPSet, Prob)
end

main()
