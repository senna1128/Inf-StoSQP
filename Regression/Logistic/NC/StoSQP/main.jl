## Load packages
using NLPModels
using JuMP
using LinearOperators
using OptimizationProblems
using MathProgBase
using ForwardDiff
using NLPModelsJuMP
using LinearAlgebra
using Distributed
using Ipopt
using DataFrames
using Glob
using DelimitedFiles
using Random
using Distributions
using NLPModelsIpopt
using MATLAB

cd("/.../StoSQP")

module Parameter
    # Parameters of StoSQP with augmented Lagrangian
    struct StoSQP
        verbose                            # Do we create dump dir?
        # stopping parameters
        MaxIter::Int                       # Maximum Iteration
        # fixed parameters
        Rep::Int                           # Number of Independent runs
        tau::Array{Int64}                  # Number of iterations for inexact solver
        c_1::Float64                       # beta_t = c_1/t^{c_2}
        c_2::Array{Float64}
        c_3::Float64                       # chi_t = beta_t^{c_3}
        SigToe::Array{Float64}             # Toeplitz covariance
        SigEqui::Array{Float64}            # Equi Corr covariance
        # Data parameters
        D::Array{Int64}                    # problem dimension
        M::Array{Float64}                  # number of constraints
    end
end
using Main.Parameter
include("../Parameter/Param.jl")

include("StoSQPMain.jl")
#######################################
#########  run main file    ###########
#######################################
function main()
    StoSQPMain(StoSQPSet)
end

main()
