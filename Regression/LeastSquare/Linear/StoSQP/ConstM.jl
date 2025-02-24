# Implement Constrained M estimator
# SampleSize: sample size
# nx, nlam: dimensions of variable and constraints
# X_true, AA: true solution, Jacobian
# RR, IdR: covariance setup and indicator of covaraince type

function ConstM(SampleSize,nx,nlam,X_true,AA,RR,IdR,EPS=1e-6,Max_Iter=1e3,rho=10)
    ## Define Covariance Matrix
    if IdR == 1
        Sigma = [RR^abs(i-j) for i=1:nx, j=1:nx] + 5*Matrix(1.0I,nx,nx)
    else
        Sigma = RR*ones(nx,nx)+(1-RR)*Matrix(1.0I,nx,nx) + 5*Matrix(1.0I,nx,nx)
    end
    ## Compute true Hessian and Constraint RHS
    K_star = hcat(vcat(Sigma,AA),vcat(AA',zeros(nlam,nlam)))
    BB = AA*X_true
    ## True value
    Lam_true = zeros(nlam)
    ## Initialization variables
    t, CompFlops, X_t, Lam_t = 1, 0, zeros(nx), ones(nlam)
    Step_size, KKT, NewDir_t, Hxx = 0, 1, zeros(nx+nlam), Matrix(1.0I,nx,nx)
    ## Offline Data Generation
    a_t = rand(MvNormal(zeros(nx),Sigma),SampleSize)'
    eps_t = rand(Normal(0,1),SampleSize)
    aa_t, aeps_t = a_t'*a_t/SampleSize, a_t'*eps_t/SampleSize
    ## Start the SQP loop
    Time = time()
    while t <= Max_Iter && KKT >= EPS
        # Objective quantities
        f_val = (X_t-X_true)'*aa_t*(X_t-X_true)/2-aeps_t'*(X_t-X_true)
        c_val = AA*X_t-BB
        grad_f = aa_t*(X_t-X_true)-aeps_t
        jac_c = AA
        CompFlops += 2*SampleSize*(nx+1)+nx*nlam
        # Solve for Newton System
        K_t = [Hxx AA';AA zeros(nlam,nlam)]
        nabL_t = [grad_f+jac_c'*Lam_t;c_val]
        KKT, NewDir_t = norm(nabL_t), lu(K_t)\-nabL_t
        CompFlops += (nx+nlam)^3+nx*nlam
        # Line search on the merit function
        Step_size = 1
        Quant1 = f_val+rho*norm(c_val,1)
        Quant2 = grad_f'*NewDir_t[1:nx]+rho*(norm(AA*X_t-BB+AA*NewDir_t[1:nx],1)-norm(c_val,1))
        CompFlops += nx*nlam
        while Step_size >= EPS
            x_new = X_t + Step_size*NewDir_t[1:nx]
            f_val_new = (x_new-X_true)'*aa_t*(x_new-X_true)/2-aeps_t'*(x_new-X_true)
            c_val_new = AA*x_new-BB
            CompFlops += SampleSize*(nx+1)+nx*nlam
            if f_val_new+rho*norm(c_val_new,1) <= Quant1+0.5*Step_size*Quant2
                break
            end
            Step_size *= 0.9
        end
        # Update primal and dual
        X_t = X_t + Step_size*NewDir_t[1:nx]
        Lam_t = Lam_t + Step_size*NewDir_t[nx+1:end]
        t = t + 1
    end
    Time = time() - Time
    if KKT >= 1e-3
        return 0,0,0,0,0,0,0,0,1,0
    else
        # Compute Final Error
        ErrX = norm(X_t-X_true)
        ErrXLam = norm([X_t-X_true,Lam_t-Lam_true])
        K_t = [aa_t AA';AA zeros(nlam, nlam)]
        ErrK = opnorm(K_t-K_star,1)
        # Compute the 95% CI radius over primal individuals
        grad_f = aa_t*(X_t-X_true)-aeps_t
        SamCov = a_t'*Diagonal((a_t*(X_t-X_true)-eps_t).^2)*a_t/SampleSize - grad_f*grad_f'
        Omega_t = inv(Matrix(K_t))*hcat(vcat(SamCov,zeros(nlam,nx)),zeros(nx+nlam,nlam))*inv(Matrix(K_t))
        CompFlops += SampleSize*nx^2+(nx+nlam)^3
        ## vectorize inference
        Radius = 1.96*sqrt.(diag(Omega_t)[1:nx]/SampleSize)
        Cov_Id = (X_true.>= X_t-Radius).*(X_true.<=X_t+Radius)

        ## Test Inf
        COV_value = sum(Omega_t)/(nx+nlam)^2
        Radius_XLam = 1.96*sqrt(COV_value/SampleSize)
        Cov_Id_XLam = mean([X_true;Lam_true])>=(mean([X_t;Lam_t])-Radius_XLam) && mean([X_true;Lam_true])<=(mean([X_t;Lam_t])+Radius_XLam)

        return ErrX,ErrXLam,ErrK,KKT,Time,CompFlops/t,Radius,Cov_Id,0,Cov_Id_XLam
    end
end
