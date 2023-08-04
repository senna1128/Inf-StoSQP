
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

function StoSQPMain(StoSQPSet)

	Max_Iter,Rep,tau = StoSQPSet.MaxIter,StoSQPSet.Rep,StoSQPSet.tau
	c_1,c_2,c_3 = StoSQPSet.c_1,StoSQPSet.c_2,StoSQPSet.c_3
	SigToe,SigEqui = StoSQPSet.SigToe,StoSQPSet.SigEqui
	D,M = StoSQPSet.D,StoSQPSet.M
	LenTau,LenStep,LenDim,LenM = length(tau),length(c_2),length(D),length(M)
	LenSigToe,LenSigEqui = length(SigToe),length(SigEqui)

	# Go over all cases
	for IdStep=1:LenStep, IdDim=1:LenDim, IdM=1:LenM, IdTau=1:LenTau
		nx, nlam = D[IdDim],Int(floor(D[IdDim]^M[IdM]))
		X_true = range(0,stop=1,length=nx)|>collect
		AA = rand(Normal(0,1),nlam,nx)
		# Toeplitz Covariance
		CountConvT, CountRadiusT = [[] for i = 1:LenSigToe], [[] for i = 1:LenSigToe]
		for IdSigToe = 1:LenSigToe, IdRep = 1:Rep
			println("Step","-",IdStep,"-D-",IdDim,"-M-",IdM,"-Tau-",IdTau,"-CovT-",IdSigToe,"-Rep-",IdRep)
			X1,Lam1,Radius,ErrXLam,ErrK,KKT,Time,IdCov,IdSing=
			StoSQP(c_1,c_2[IdStep],c_3,Max_Iter,tau[IdTau],
			nx,nlam,X_true,AA,SigToe[IdSigToe],1)
			if IdSing == 0
				push!(CountConvT[IdSigToe],IdCov)
				push!(CountRadiusT[IdSigToe],Radius[end])
				path1 = string("../Solution/Step",IdStep,
				"Dim",IdDim,"M",IdM,"Tau",IdTau,"/Toe",IdSigToe)
				if !isdir(path1)
					mkpath(path1)
				end
				path = string(path1,"/rep",IdRep,".mat")
				Result = StoSQPResult(X_true,zeros(nlam),X1,Lam1,Int.(IdCov),Radius,ErrXLam,ErrK,KKT,Time)
				write_matfile(path; Result)
			end
		end
		write_matfile(string("../Solution/Step",IdStep,
				"Dim",IdDim,"M",IdM,"Tau",IdTau,"/ToeCR.mat");CountConvT,CountRadiusT)
		# Equi Covariance
		CountConvE, CountRadiusE = [[] for i = 1:LenSigEqui], [[] for i = 1:LenSigEqui]
		for IdSigEqui = 1:LenSigEqui, IdRep = 1:Rep
			println("Step","-",IdStep,"-D-",IdDim,"-M-",IdM,"-Tau-",IdTau,"-CovE-",IdSigEqui,"-Rep-",IdRep)
			X1,Lam1,Radius,ErrXLam,ErrK,KKT,Time,IdCov,IdSing=
			StoSQP(c_1,c_2[IdStep],c_3,Max_Iter,tau[IdTau],
			nx,nlam,X_true,AA,SigEqui[IdSigEqui],2)
			if IdSing == 0
				push!(CountConvE[IdSigEqui],IdCov)
				push!(CountRadiusE[IdSigEqui],Radius[end])
				path1 = string("../Solution/Step",IdStep,
				"Dim",IdDim,"M",IdM,"Tau",IdTau,"/Equ",IdSigEqui)
				if !isdir(path1)
					mkpath(path1)
				end
				path = string(path1,"/rep",IdRep,".mat")
				Result = StoSQPResult(X_true,zeros(nlam),X1,Lam1,Int.(IdCov),Radius,ErrXLam,ErrK,KKT,Time)
				write_matfile(path; Result)
			end
		end
		write_matfile(string("../Solution/Step",IdStep,
				"Dim",IdDim,"M",IdM,"Tau",IdTau,"/EquCR.mat");CountConvE,CountRadiusE)
	end

end
