
include("StoSQP.jl")
struct StoSQPResult
	X_true::Array
	Lam_true::Array
	X1::Array
    Lam1::Array
	CoverRate::Int64
	Radius::Array
	ErrXLam::Array
	ErrK::Array
    KKT::Array
    Time::Float64
end


## Implement StoSQP framework for whole problem set
# StoSQPSet: parameters of Berahas algorithm
# Prob: problem name set

function StoSQPMain(StoSQPSet,Prob)
	Max_Iter = StoSQPSet.MaxIter
	Rep = StoSQPSet.Rep
	tau = StoSQPSet.tau
	c_1 = StoSQPSet.c_1
	c_2 = StoSQPSet.c_2
	c_3 = StoSQPSet.c_3
	Sigma = StoSQPSet.Sigma
	LenDStep = length(c_2)
	LenSigma = length(Sigma)
	StoSQPR = Array{StoSQPResult}(undef,length(Prob),LenDStep,LenSigma,Rep)

	# Go over all Problems
	for Idprob = 1:length(Prob)
		# load problem
		nlp = CUTEstModel(Prob[Idprob])
		# solve exactly
        stats = ipopt(nlp,print_level=0)
		# check if the nlp model is solved
		if stats.solver_specific[:internal_msg] == :Solve_Succeeded
			X_true, Lam_true = stats.solution, stats.multipliers
			# Implement StoSQP framework
			# Loop over decay stepsize, sigma, replicate
			i = 1
	        while i <= LenDStep
	            j = 1
	            while j <= LenSigma
	                rep = 1
	                while rep <= Rep
						X1,Lam1,Radius,ErrXLam,ErrK,KKT,Time,IdCov,IdSing = StoSQP(nlp,c_1,c_2[i],c_3,Sigma[j],Max_Iter,tau,X_true,Lam_true)
						if IdSing == 1
							break
						else
							StoSQPR[Idprob,i,j,rep] = StoSQPResult(X_true,Lam_true,X1,Lam1,Int.(IdCov),Radius,ErrXLam,ErrK,KKT,Time)
							path1 = string("../Solution/StoSQP",Idprob,"/D",i,"S",j)
							if !isdir(path1)
								mkpath(path1)
							end
							path = string(path1,"/rep",rep,".mat")
							Result = StoSQPR[Idprob,i,j,rep]
						    write_matfile(path; Result)
							rep += 1
						end
					end
					j += 1
				end
				i += 1
			end
		end
		finalize(nlp)
	end


	return StoSQPR
end
