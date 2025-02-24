
include("StoSQP.jl")

struct StoSQPResult
	ErrX::Array
	ErrXLam::Array
	ErrK::Array
    KKT::Array
	Time::Float64
	Per_Flops::Float64
	Radius::Array
	Cov_Id::Array
end

struct CaseResult
	Case_MSE_X::Array
	Case_MSE_XLam::Array
	Case_MSE_K::Array
	Case_MSE_KKT::Array
	Case_Run_Time::Array
	Case_Ave_Flops::Array
	Case_Radius_Entry::Array
	Case_Radius_Ave::Array
	Case_Cov_Id_Entry::Array
	Case_Cov_Id_Ave::Array
	Case_Cov_Id_Ave_XLam::Array
end




## Implement StoSQP framework for whole problem set
# StoSQPSet: parameters of Berahas algorithm
# Prob: problem name set

function StoSQPMain(StoSQPSet,Prob)
	Max_Iter,Rep,tau = StoSQPSet.MaxIter,StoSQPSet.Rep,StoSQPSet.tau
	c_1,c_2,c_3,Sigma = StoSQPSet.c_1,StoSQPSet.c_2,StoSQPSet.c_3,StoSQPSet.Sigma
	LenDStep, LenSigma = length(c_2), length(Sigma)

	# Go over all Problems
	for Idprob = 5:5
		Random.seed!(2023)
		# load problem
		nlp = CUTEstModel(Prob[Idprob])
		# solve exactly
        stats = ipopt(nlp,print_level=0)
		# check if the nlp model is solved
		if stats.solver_specific[:internal_msg] == :Solve_Succeeded
			X_true,Lam_true,nx = stats.solution,stats.multipliers,nlp.meta.nvar
			G_star, B_star = jac(nlp,X_true), hess(nlp,X_true,Lam_true)
			path = string("../Solution/Prob",Idprob)
			if !isdir(path)
				mkpath(path)
			end
			write_matfile(string(path,"/SaveSolution.mat");X_true,Lam_true,G_star,B_star)
			# Implement StoSQP framework
			# Loop over decay stepsize, sigma, replicate
			XSeq = [[] for i = 1:LenDStep, j = 1:LenSigma]
			XLamSeq = [[] for i = 1:LenDStep, j = 1:LenSigma]
			KSeq = [[] for i = 1:LenDStep, j = 1:LenSigma]
			KKTSeq = [[] for i = 1:LenDStep, j = 1:LenSigma]
			TimeSeq = [[] for i = 1:LenDStep, j = 1:LenSigma]
			FlopsSeq = [[] for i = 1:LenDStep, j = 1:LenSigma]
			RadiusSeq = [Array{Float64}(undef,nx,0) for i = 1:LenDStep, j = 1:LenSigma]
			Cov_IdSeq = [Array{Float64}(undef,nx,0) for i = 1:LenDStep, j = 1:LenSigma]

			MSE_X = zeros(LenDStep,LenSigma,2)
			MSE_XLam = zeros(LenDStep,LenSigma,2)
			MSE_K = zeros(LenDStep,LenSigma,2)
			MSE_KKT = zeros(LenDStep,LenSigma,2)
			Run_Time = zeros(LenDStep,LenSigma,2)
			Ave_Flops = zeros(LenDStep,LenSigma,2)
			Radius_Entry = Array{Any}(undef,LenDStep,LenSigma)
			Radius_Ave = zeros(LenDStep,LenSigma,2)
			Cov_Id_Entry = Array{Any}(undef,LenDStep,LenSigma)
			Cov_Id_Ave = zeros(LenDStep,LenSigma,2)

			# Go over all cases
			for i = 1:LenDStep, j = 1:LenSigma
				for rep = 1:Rep
					println("StoSQP Prob-",Idprob,"-Step-",i,"-Sig-",j,"-Rep-",rep)
					ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,IdSing =
					StoSQP(nlp,c_1,c_2[i],c_3,Sigma[j],Max_Iter,tau,X_true,Lam_true)
					if IdSing == 1
						break
					else
						## Save Result
						Result = StoSQPResult(ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,Int.(Cov_Id_XLam))
						path1 = string("../Solution/Prob",Idprob,"/D",i,"S",j)
						if !isdir(path1)
							mkpath(path1)
						end
						path = string(path1,"/rep",rep,".mat")
						write_matfile(path; Result)

						push!(XSeq[i,j],ErrX[end])
						push!(XLamSeq[i,j],ErrXLam[end])
						push!(KSeq[i,j],ErrK[end])
						push!(KKTSeq[i,j],KKT[end])
						push!(TimeSeq[i,j],Time)
						push!(FlopsSeq[i,j],Per_Flops)
						RadiusSeq[i,j] = hcat(RadiusSeq[i,j],Radius)
						Cov_IdSeq[i,j] = hcat(Cov_IdSeq[i,j],Cov_Id)
					end
				end
				## Compute average result
				MSE_X[i,j,:] = [mean(XSeq[i,j]); std(XSeq[i,j])]
				MSE_XLam[i,j,:] = [mean(XLamSeq[i,j]); std(XLamSeq[i,j])]
				MSE_K[i,j,:] = [mean(KSeq[i,j]); std(KSeq[i,j])]
				MSE_KKT[i,j,:] = [mean(KKTSeq[i,j]); std(KKTSeq[i,j])]
				Run_Time[i,j,:] = [mean(TimeSeq[i,j]); std(TimeSeq[i,j])]
				Ave_Flops[i,j,:] = [mean(FlopsSeq[i,j]); std(FlopsSeq[i,j])]
				Radius_Entry[i,j] = [mean(RadiusSeq[i,j],dims=2) std(RadiusSeq[i,j],dims=2)]
				Radius_Ave[i,j,:] = [mean(RadiusSeq[i,j]); std(RadiusSeq[i,j])]
				Cov_Id_Entry[i,j] = [mean(Cov_IdSeq[i,j],dims=2) std(Cov_IdSeq[i,j],dims=2)]
				Cov_Id_Ave[i,j,:] = [mean(Cov_IdSeq[i,j]); std(Cov_IdSeq[i,j])]

				Case_Result = CaseResult(MSE_X[i,j,:],MSE_XLam[i,j,:],
				MSE_K[i,j,:],MSE_KKT[i,j,:],
				Run_Time[i,j,:],Ave_Flops[i,j,:],
				Radius_Entry[i,j],Radius_Ave[i,j,:],
				Cov_Id_Entry[i,j],Cov_Id_Ave[i,j,:])
				path1 = string("../Solution/Prob",Idprob,"/D",i,"S",j)
				path = string(path1,"/CaseSummary.mat")
				write_matfile(path;Case_Result)
			end
			path = string("../Solution/Prob",Idprob,"/Summary.mat")
			write_matfile(path;MSE_X,MSE_XLam,MSE_K,MSE_KKT,Run_Time,Ave_Flops,Radius_Entry,Radius_Ave,
			Cov_Id_Entry,Cov_Id_Ave)
		end
		finalize(nlp)
		
	end


end
