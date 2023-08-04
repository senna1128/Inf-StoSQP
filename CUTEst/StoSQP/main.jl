
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

cd("/.../StoSQP")
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
        c_1::Float64                       # beta_t = c_1/t^{c_2}
        c_2::Array{Float64}
        c_3::Float64                       # chi_t = beta_t^{c_3}
        # test parameters
        Sigma::Array{Float64}              # variance of gradient
    end
end
using Main.Parameter
include("../Parameter/Param.jl")
include("StoSQPMain.jl")


#######################################
#########  run main file    ###########
#######################################
function main()
    ## run StoSQP framework
    StoSQPR = StoSQPMain(StoSQPSet, Prob)
end

main()

