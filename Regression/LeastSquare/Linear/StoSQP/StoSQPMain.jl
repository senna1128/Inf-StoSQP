
include("StoSQP.jl")
include("ConstM.jl")

struct StoSQPResult
	ErrX::Float64
	ErrXLam::Float64
	ErrK::Float64
    KKT::Float64
	Time::Float64
	Per_Flops::Float64
	Radius::Array
	Cov_Id::Array
	Cov_Id_XLam::Int64
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

struct TrueSolution
	X_true::Array
	AA::Array
end

## Implement StoSQP framework for whole problem set
# StoSQPSet: parameters of Berahas algorithm

function StoSQPMain(StoSQPSet, ConstMSet)
	## Extract parameters
	Max_Iter,tau = StoSQPSet.MaxIter,StoSQPSet.tau
	c_1,c_2,c_3 = StoSQPSet.c_1,StoSQPSet.c_2,StoSQPSet.c_3
	SigToe,SigEqui = StoSQPSet.SigToe,StoSQPSet.SigEqui
	D,M = StoSQPSet.D,StoSQPSet.M
	SampleSize = ConstMSet.SampleSize
	Rep = max(StoSQPSet.Rep,ConstMSet.Rep)
	LenTau,LenStep,LenDim,LenM = length(tau),length(c_2),length(D),length(M)
	LenSigToe,LenSigEqui = length(SigToe),length(SigEqui)

	# Go over all cases
	for IdStep=1:LenStep, IdDim=1:4, IdM=1:LenM
		Random.seed!(2023)
		nx, nlam = D[IdDim],Int(floor(D[IdDim]^M[IdM]))
		X_true = range(0,stop=1,length=nx)|>collect
		AA = rand(Normal(0,0.5),nlam,nx)
		# Save True Solution
		path = string("../Solution/S",IdStep,"D",IdDim,"M",IdM)
		if !isdir(path)
			mkpath(path)
		end
		True_Solution = TrueSolution(X_true,AA)
		write_matfile(string(path,"/SaveSolution.mat");True_Solution)
		# Result Array (the second column is the start of StoSQP with diff tau)
		MSE_X = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		MSE_XLam = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		MSE_K = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		MSE_KKT = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		Run_Time = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		Ave_Flops = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		Radius_Entry = Array{Any}(undef,LenSigToe+LenSigEqui,1+LenTau)
		Radius_Ave = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		Cov_Id_Entry = Array{Any}(undef,LenSigToe+LenSigEqui,1+LenTau)
		Cov_Id_Ave = zeros(LenSigToe+LenSigEqui,1+LenTau,2)
		Cov_Id_Ave_XLam = zeros(LenSigToe+LenSigEqui,1+LenTau,2)

		## Toeplitz Covariance
		for IdSigToe = 1:LenSigToe
			## Constrained M Solver
			# Result Array
			XSeq, XLamSeq, KSeq, KKTSeq = [], [], [], []
			TimeSeq, FlopsSeq, Cov_IdSeq_XLam = [], [], []
			RadiusSeq, Cov_IdSeq = Array{Float64}(undef,nx,0), Array{Float64}(undef,nx,0)
			for IdRep = 1:Rep
				println("ConstM-Step-",IdStep,"-M-",IdM,"-D-",IdDim,"-CovT-",IdSigToe,"-Rep-",IdRep)
				ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,IdSing,Cov_Id_XLam =
				ConstM(SampleSize,nx,nlam,X_true,AA,SigToe[IdSigToe],1)
				if IdSing == 0
					## Save Result
					Result = StoSQPResult(ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,Int.(Cov_Id_XLam))
					path1 = string("../Solution/S",IdStep,"D",IdDim,
					"M",IdM,"/ConstM/Toe",IdSigToe)
					if !isdir(path1)
						mkpath(path1)
					end
					path = string(path1,"/rep",IdRep,".mat")
					write_matfile(path; Result)

					push!(XSeq,ErrX)
					push!(XLamSeq,ErrXLam)
					push!(KSeq,ErrK)
					push!(KKTSeq,KKT)
					push!(TimeSeq,Time)
					push!(FlopsSeq,Per_Flops)
					RadiusSeq = hcat(RadiusSeq,Radius)
					Cov_IdSeq = hcat(Cov_IdSeq,Cov_Id)
					push!(Cov_IdSeq_XLam,Int.(Cov_Id_XLam))
				end
			end
			## Compute average case
			MSE_X[IdSigToe,1,:] = [mean(XSeq); std(XSeq)]
			MSE_XLam[IdSigToe,1,:] = [mean(XLamSeq); std(XLamSeq)]
			MSE_K[IdSigToe,1,:] = [mean(KSeq); std(KSeq)]
			MSE_KKT[IdSigToe,1,:] = [mean(KKTSeq); std(KKTSeq)]
			Run_Time[IdSigToe,1,:] = [mean(TimeSeq); std(TimeSeq)]
			Ave_Flops[IdSigToe,1,:] = [mean(FlopsSeq); std(FlopsSeq)]
			Radius_Entry[IdSigToe,1] = [mean(RadiusSeq,dims=2) std(RadiusSeq,dims=2)]
			Radius_Ave[IdSigToe,1,:] = [mean(RadiusSeq); std(RadiusSeq)]
			Cov_Id_Entry[IdSigToe,1] = [mean(Cov_IdSeq,dims=2) std(Cov_IdSeq,dims=2)]
			Cov_Id_Ave[IdSigToe,1,:] = [mean(Cov_IdSeq); std(Cov_IdSeq)]
			Cov_Id_Ave_XLam[IdSigToe,1,:] = [mean(Cov_IdSeq_XLam); std(Cov_IdSeq_XLam)]
			## Save Case Result
			Case_Result = CaseResult(MSE_X[IdSigToe,1,:],MSE_XLam[IdSigToe,1,:],
			MSE_K[IdSigToe,1,:],MSE_KKT[IdSigToe,1,:],
			Run_Time[IdSigToe,1,:],Ave_Flops[IdSigToe,1,:],
			Radius_Entry[IdSigToe,1],Radius_Ave[IdSigToe,1,:],
			Cov_Id_Entry[IdSigToe,1],Cov_Id_Ave[IdSigToe,1,:],
			Cov_Id_Ave_XLam[IdSigToe,1,:])
			path1 = string("../Solution/S",IdStep,"D",IdDim,
			"M",IdM,"/ConstM/Toe",IdSigToe)
			path = string(path1,"/CaseSummary.mat")
			write_matfile(path;Case_Result)

			## StoSQP Methods
			for IdTau=1:LenTau
				# Result Array
				XSeq, XLamSeq, KSeq, KKTSeq = [], [], [], []
				TimeSeq, FlopsSeq, Cov_IdSeq_XLam = [], [], []
				RadiusSeq, Cov_IdSeq = Array{Float64}(undef,nx,0), Array{Float64}(undef,nx,0)
				for IdRep = 1:Rep
					println("StoSQP-Step-",IdStep,"-M-",IdM,"-D-",IdDim,"-Tau-",IdTau,"-CovT-",IdSigToe,"-Rep-",IdRep)
					ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,IdSing,Cov_Id_XLam =
					StoSQP(c_1,c_2[IdStep],c_3,Max_Iter,tau[IdTau],
					nx,nlam,X_true,AA,SigToe[IdSigToe],1)
					if IdSing == 0
						## Save Result
						Result = StoSQPResult(ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,Int.(Cov_Id_XLam))
						path1 = string("../Solution/S",IdStep,"D",IdDim,
						"M",IdM,"/StoSQPTau",IdTau,"/Toe",IdSigToe)
						if !isdir(path1)
							mkpath(path1)
						end
						path = string(path1,"/rep",IdRep,".mat")
						write_matfile(path; Result)

						push!(XSeq,ErrX)
						push!(XLamSeq,ErrXLam)
						push!(KSeq,ErrK)
						push!(KKTSeq,KKT)
						push!(TimeSeq,Time)
						push!(FlopsSeq,Per_Flops)
						RadiusSeq = hcat(RadiusSeq,Radius)
						Cov_IdSeq = hcat(Cov_IdSeq,Cov_Id)
						push!(Cov_IdSeq_XLam,Int.(Cov_Id_XLam))
					end
				end
				## Compute average case
				MSE_X[IdSigToe,1+IdTau,:] = [mean(XSeq); std(XSeq)]
				MSE_XLam[IdSigToe,1+IdTau,:] = [mean(XLamSeq); std(XLamSeq)]
				MSE_K[IdSigToe,1+IdTau,:] = [mean(KSeq); std(KSeq)]
				MSE_KKT[IdSigToe,1+IdTau,:] = [mean(KKTSeq); std(KKTSeq)]
				Run_Time[IdSigToe,1+IdTau,:] = [mean(TimeSeq); std(TimeSeq)]
				Ave_Flops[IdSigToe,1+IdTau,:] = [mean(FlopsSeq); std(FlopsSeq)]
				Radius_Entry[IdSigToe,1+IdTau] = [mean(RadiusSeq,dims=2) std(RadiusSeq,dims=2)]
				Radius_Ave[IdSigToe,1+IdTau,:] = [mean(RadiusSeq); std(RadiusSeq)]
				Cov_Id_Entry[IdSigToe,1+IdTau] = [mean(Cov_IdSeq,dims=2) std(Cov_IdSeq,dims=2)]
				Cov_Id_Ave[IdSigToe,1+IdTau,:] = [mean(Cov_IdSeq); std(Cov_IdSeq)]
				Cov_Id_Ave_XLam[IdSigToe,1+IdTau,:] = [mean(Cov_IdSeq_XLam); std(Cov_IdSeq_XLam)]
				## Save Case Result
				Case_Result = CaseResult(MSE_X[IdSigToe,1+IdTau,:],MSE_XLam[IdSigToe,1+IdTau,:],
				MSE_K[IdSigToe,1+IdTau,:],MSE_KKT[IdSigToe,1+IdTau,:],
				Run_Time[IdSigToe,1+IdTau,:],Ave_Flops[IdSigToe,1+IdTau,:],
				Radius_Entry[IdSigToe,1+IdTau],Radius_Ave[IdSigToe,1+IdTau,:],
				Cov_Id_Entry[IdSigToe,1+IdTau],Cov_Id_Ave[IdSigToe,1+IdTau,:],
				Cov_Id_Ave_XLam[IdSigToe,1+IdTau,:])
				path1 = string("../Solution/S",IdStep,"D",IdDim,
				"M",IdM,"/StoSQPTau",IdTau,"/Toe",IdSigToe)
				path = string(path1,"/CaseSummary.mat")
				write_matfile(path;Case_Result)
			end
		end

		## Equi Covariance
		for IdSigEqui = 1:LenSigEqui
			## Constrained M Solver
			# Result Array
			XSeq, XLamSeq, KSeq, KKTSeq = [], [], [], []
			TimeSeq, FlopsSeq, Cov_IdSeq_XLam = [], [], []
			RadiusSeq, Cov_IdSeq = Array{Float64}(undef,nx,0), Array{Float64}(undef,nx,0)
			for IdRep = 1:Rep
				println("ConstM-Step-",IdStep,"-M-",IdM,"-D-",IdDim,"-CovE-",IdSigEqui,"-Rep-",IdRep)
				ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,IdSing,Cov_Id_XLam =
				ConstM(SampleSize,nx,nlam,X_true,AA,SigEqui[IdSigEqui],2)
				if IdSing == 0
					## Save Result
					Result = StoSQPResult(ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,Int.(Cov_Id_XLam))
					path1 = string("../Solution/S",IdStep,"D",IdDim,
					"M",IdM,"/ConstM/Equ",IdSigEqui)
					if !isdir(path1)
						mkpath(path1)
					end
					path = string(path1,"/rep",IdRep,".mat")
					write_matfile(path; Result)

					push!(XSeq,ErrX)
					push!(XLamSeq,ErrXLam)
					push!(KSeq,ErrK)
					push!(KKTSeq,KKT)
					push!(TimeSeq,Time)
					push!(FlopsSeq,Per_Flops)
					RadiusSeq = hcat(RadiusSeq,Radius)
					Cov_IdSeq = hcat(Cov_IdSeq,Cov_Id)
					push!(Cov_IdSeq_XLam,Int.(Cov_Id_XLam))
				end
			end
			## Compute average case
			MSE_X[LenSigToe+IdSigEqui,1,:] = [mean(XSeq); std(XSeq)]
			MSE_XLam[LenSigToe+IdSigEqui,1,:] = [mean(XLamSeq); std(XLamSeq)]
			MSE_K[LenSigToe+IdSigEqui,1,:] = [mean(KSeq); std(KSeq)]
			MSE_KKT[LenSigToe+IdSigEqui,1,:] = [mean(KKTSeq); std(KKTSeq)]
			Run_Time[LenSigToe+IdSigEqui,1,:] = [mean(TimeSeq); std(TimeSeq)]
			Ave_Flops[LenSigToe+IdSigEqui,1,:] = [mean(FlopsSeq); std(FlopsSeq)]
			Radius_Entry[LenSigToe+IdSigEqui,1] = [mean(RadiusSeq,dims=2) std(RadiusSeq,dims=2)]
			Radius_Ave[LenSigToe+IdSigEqui,1,:] = [mean(RadiusSeq); std(RadiusSeq)]
			Cov_Id_Entry[LenSigToe+IdSigEqui,1] = [mean(Cov_IdSeq,dims=2) std(Cov_IdSeq,dims=2)]
			Cov_Id_Ave[LenSigToe+IdSigEqui,1,:] = [mean(Cov_IdSeq); std(Cov_IdSeq)]
			Cov_Id_Ave_XLam[LenSigToe+IdSigEqui,1,:] = [mean(Cov_IdSeq_XLam); std(Cov_IdSeq_XLam)]
			## Save Case Result
			Case_Result = CaseResult(MSE_X[LenSigToe+IdSigEqui,1,:],MSE_XLam[LenSigToe+IdSigEqui,1,:],
			MSE_K[LenSigToe+IdSigEqui,1,:],MSE_KKT[LenSigToe+IdSigEqui,1,:],
			Run_Time[LenSigToe+IdSigEqui,1,:],Ave_Flops[LenSigToe+IdSigEqui,1,:],
			Radius_Entry[LenSigToe+IdSigEqui,1],Radius_Ave[LenSigToe+IdSigEqui,1,:],
			Cov_Id_Entry[LenSigToe+IdSigEqui,1],Cov_Id_Ave[LenSigToe+IdSigEqui,1,:],
			Cov_Id_Ave_XLam[LenSigToe+IdSigEqui,1,:])
			path1 = string("../Solution/S",IdStep,"D",IdDim,
			"M",IdM,"/ConstM/Equ",IdSigEqui)
			path = string(path1,"/CaseSummary.mat")
			write_matfile(path;Case_Result)

			## StoSQP Methods
			for IdTau=1:LenTau
				# Result Array
				XSeq, XLamSeq, KSeq, KKTSeq = [], [], [], []
				TimeSeq, FlopsSeq, Cov_IdSeq_XLam = [], [], []
				RadiusSeq, Cov_IdSeq = Array{Float64}(undef,nx,0), Array{Float64}(undef,nx,0)
				for IdRep = 1:Rep
					println("StoSQP-Step-",IdStep,"-M-",IdM,"-D-",IdDim,"-Tau-",IdTau,"-CovE-",IdSigEqui,"-Rep-",IdRep)
					ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,IdSing,Cov_Id_XLam =
					StoSQP(c_1,c_2[IdStep],c_3,Max_Iter,tau[IdTau],
					nx,nlam,X_true,AA,SigEqui[IdSigEqui],2)
					if IdSing == 0
						## Save Result
						Result = StoSQPResult(ErrX,ErrXLam,ErrK,KKT,Time,Per_Flops,Radius,Cov_Id,Int.(Cov_Id_XLam))
						path1 = string("../Solution/S",IdStep,"D",IdDim,
						"M",IdM,"/StoSQPTau",IdTau,"/Equ",IdSigEqui)
						if !isdir(path1)
							mkpath(path1)
						end
						path = string(path1,"/rep",IdRep,".mat")
						write_matfile(path; Result)

						push!(XSeq,ErrX)
						push!(XLamSeq,ErrXLam)
						push!(KSeq,ErrK)
						push!(KKTSeq,KKT)
						push!(TimeSeq,Time)
						push!(FlopsSeq,Per_Flops)
						RadiusSeq = hcat(RadiusSeq,Radius)
						Cov_IdSeq = hcat(Cov_IdSeq,Cov_Id)
						push!(Cov_IdSeq_XLam,Int.(Cov_Id_XLam))
					end
				end
				## Compute average case
				MSE_X[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(XSeq); std(XSeq)]
				MSE_XLam[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(XLamSeq); std(XLamSeq)]
				MSE_K[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(KSeq); std(KSeq)]
				MSE_KKT[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(KKTSeq); std(KKTSeq)]
				Run_Time[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(TimeSeq); std(TimeSeq)]
				Ave_Flops[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(FlopsSeq); std(FlopsSeq)]
				Radius_Entry[LenSigToe+IdSigEqui,1+IdTau] = [mean(RadiusSeq,dims=2) std(RadiusSeq,dims=2)]
				Radius_Ave[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(RadiusSeq); std(RadiusSeq)]
				Cov_Id_Entry[LenSigToe+IdSigEqui,1+IdTau] = [mean(Cov_IdSeq,dims=2) std(Cov_IdSeq,dims=2)]
				Cov_Id_Ave[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(Cov_IdSeq); std(Cov_IdSeq)]
				Cov_Id_Ave_XLam[LenSigToe+IdSigEqui,1+IdTau,:] = [mean(Cov_IdSeq_XLam); std(Cov_IdSeq_XLam)]
				## Save Case Result
				Case_Result = CaseResult(MSE_X[LenSigToe+IdSigEqui,1+IdTau,:],MSE_XLam[LenSigToe+IdSigEqui,1+IdTau,:],
				MSE_K[LenSigToe+IdSigEqui,1+IdTau,:],MSE_KKT[LenSigToe+IdSigEqui,1+IdTau,:],
				Run_Time[LenSigToe+IdSigEqui,1+IdTau,:],Ave_Flops[LenSigToe+IdSigEqui,1+IdTau,:],
				Radius_Entry[LenSigToe+IdSigEqui,1+IdTau],Radius_Ave[LenSigToe+IdSigEqui,1+IdTau,:],
				Cov_Id_Entry[LenSigToe+IdSigEqui,1+IdTau],Cov_Id_Ave[LenSigToe+IdSigEqui,1+IdTau,:],
				Cov_Id_Ave_XLam[LenSigToe+IdSigEqui,1+IdTau,:])
				path1 = string("../Solution/S",IdStep,"D",IdDim,
				"M",IdM,"/StoSQPTau",IdTau,"/Equ",IdSigEqui)
				path = string(path1,"/CaseSummary.mat")
				write_matfile(path;Case_Result)
			end
		end
		path = string("../Solution/S",IdStep,"D",IdDim,"M",IdM,"/Summary.mat")
		write_matfile(path;MSE_X,MSE_XLam,MSE_K,MSE_KKT,Run_Time,Ave_Flops,Radius_Entry,Radius_Ave,
		Cov_Id_Entry,Cov_Id_Ave,Cov_Id_Ave_XLam)
	end

end
