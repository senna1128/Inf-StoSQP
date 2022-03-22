include("StoSQP.jl")
struct StoSQPResult
	X_true::Array
	Lam_true::Array
	X1::Array
    Lam1::Array
	Radius::Array
	ErrXLam::Array
	ErrK::Array
    KKT::Array
    Time::Array
end


## Implement StoSQP framework for whole problem set
# StoSQPSet: parameters of Berahas algorithm
# Prob: problem name set

function StoSQPMain(StoSQPSet, Prob)
	Max_Iter = StoSQPSet.MaxIter
	Rep = StoSQPSet.Rep
	tau = StoSQPSet.tau
	DBeta = StoSQPSet.DBeta
	Sigma = StoSQPSet.Sigma
    LenDStep = length(DBeta)
    LenSigma = length(Sigma)

    StoSQPR = Array{StoSQPResult}(undef, length(Prob))

    # Go over all Problems
	Idprob = 1
	while Idprob <= length(Prob)
		# load problem
        nlp = CUTEstModel(Prob[Idprob])
        # solve exactly
        stats = ipopt(nlp, print_level=0)
		# check if the nlp model is solved
		if stats.solver_specific[:internal_msg] != :Solve_Succeeded
			StoSQPR[Idprob] = StoSQPResult([0],[0],[0],[0],[0],[0],[0],[0],[0])
			path = string("../Solution/StoSQP", Idprob, ".mat")
		    Result = StoSQPR[Idprob]
		    write_matfile(path; Result)
		else
			X_true, Lam_true = stats.solution, stats.multipliers
			# Implement StoSQP framework
			# define results vector for decay stepsize
	        X1Step = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
	        Lam1Step = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
			RadiusStep = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
			ErrXLamStep = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
			ErrKStep = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
	        KKTStep = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
	        TimeStep = reshape([[] for i=1:LenDStep for j=1:LenSigma],(LenDStep,LenSigma))
			# Loop over decay stepsize, sigma, replicate
			i = 1
	        while i <= LenDStep
	            j = 1
	            while j <= LenSigma
	                rep = 1
	                while rep <= Rep
						println("StoSQP DecayStep","-",Idprob,"-",i,"-",j,"-",rep)
						X1, Lam1, Radius, ErrXLam, ErrK, KKT, Time, IdSing = StoSQP(nlp,DBeta[i],Sigma[j],Max_Iter,tau,X_true,Lam_true)
						if IdSing == 1
							break
						else
							push!(X1Step[i,j],X1)
							push!(Lam1Step[i,j],Lam1)
							push!(RadiusStep[i,j],Radius)
							push!(ErrXLamStep[i,j],ErrXLam)
							push!(ErrKStep[i,j],ErrK)
							push!(KKTStep[i,j],KKT)
							push!(TimeStep[i,j],Time)
							rep += 1
						end
					end
					j += 1
				end
				i += 1
			end
			StoSQPR[Idprob] = StoSQPResult(X_true,Lam_true,X1Step,Lam1Step,RadiusStep,ErrXLamStep,ErrKStep,KKTStep,TimeStep)
			path = string("../Solution/StoSQP", Idprob, ".mat")
		    Result = StoSQPR[Idprob]
		    write_matfile(path; Result)
		end
		finalize(nlp)
		Idprob += 1
	end

	return StoSQPR
end
