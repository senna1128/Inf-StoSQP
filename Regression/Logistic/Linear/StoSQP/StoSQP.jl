# Impmement Stochastic SQP framework
# c1,c2,c3: step decay ratio
# Max_Iter: maximum iterations
# ttau: number of inner iterations of the solver
# nx,nlam: dimension of variable and constraint
# X_true, AA: optimal solution, constraints
# RR, IdR: covariance setup and indicator of covaraince type
# Output:
# X1, Lam1, Radius, ErrXLam, ErrK, KKT, Time, IdCoverage, IdSing


function StoSQP(c1,c2,c3,Max_Iter,ttau,nx,nlam,X_true,AA,RR,IdR,EPS=1e-8,Pen=2)
    ## Set Buffer Size
    Buff= Int(floor(Max_Iter/1.1))
    ## Define Covariance Matrix
    if IdR == 1
        Sigma = [RR^abs(i-j) for i=1:nx, j=1:nx] + 5*Matrix(1.0I,nx,nx)
    else
        Sigma = RR*ones(nx,nx)+(1-RR)*Matrix(1.0I,nx,nx) + 5*Matrix(1.0I,nx,nx)
    end
    ## Compute K_star
    BB = AA*X_true
    ## True value
    Lam_true = zeros(nlam)
    ## Initialize variables
    t,CompFlops,X_t,Lam_t,Step_size = 1,0,zeros(nx),ones(nlam),0
    cum_barg_t,NewDir_t,K_t = zeros(nx),zeros(nx+nlam),zeros(nx+nlam,nx+nlam)
    cum_barggT_t,cum_bar2x_f_t = zeros(nx,nx),Matrix(1.0I,nx,nx)
    ## Start the loop
    Time = time()
    while t <= Max_Iter+1
        # Step 1: generate random gradient and Hessian
        a_t = rand(MvNormal(zeros(nx),Sigma))
        b_t = 2*rand(Bernoulli(1/(1+exp(-a_t'*X_t))))-1
        barg_t = -b_t/(1+exp(b_t*a_t'*X_t))*a_t
        if t > Buff
            cum_barggT_t = (t-Buff)/(t-Buff+1)*cum_barggT_t + 1/(t-Buff+1)*barg_t*barg_t'
            cum_barg_t = (t-Buff)/(t-Buff+1)*cum_barg_t + 1/(t-Buff+1)*barg_t
        end
        CompFlops += nx+nx^2
        # define stochastic Lagrangian gradient and Hessian
        c_t = AA*X_t - BB
        G_tTLam_t = AA'*Lam_t
        bar_nab_L_t = [barg_t+G_tTLam_t+Pen*(X_t-X_true); c_t]
        cum_bar2x_L_t = cum_bar2x_f_t + Pen*Matrix(1.0I,nx,nx)
        K_t = [cum_bar2x_L_t AA'; AA zeros(nlam,nlam)]
        # Update Hessian
        bar_nab_x2f_t = a_t*a_t'/((1+exp(a_t'*X_t))*(1+exp(-a_t'*X_t)))
        cum_bar2x_f_t = t/(t+1)*cum_bar2x_f_t + 1/(t+1)*bar_nab_x2f_t
        CompFlops += 2*nx*nlam + 2*nx^2
        # Step 2: solve Newton system via Randomized Solve
        if ttau == 0
            NewDir_t = lu(K_t)\-bar_nab_L_t
            CompFlops += (nx+nlam)^3
        else
            NewDir_t = zeros(nx+nlam)
            for inner_iter = 1:ttau
                # random pick a index
                j = sample(1:(nx+nlam))
                NewDir_t = NewDir_t - (K_t[j,:]'*NewDir_t+bar_nab_L_t[j])/norm(K_t[j,:])^2*K_t[:,j]
            end
            CompFlops += ttau*(nx+nlam)
        end
        # Step 3: update the iterate
        beta_t = c1/t^c2
        chi_t = beta_t^c3
        Step_size = rand(Uniform(beta_t,beta_t+chi_t))
        X_t = X_t + Step_size*NewDir_t[1:nx]
        Lam_t = Lam_t + Step_size*NewDir_t[nx+1:end]
        t = t + 1
    end
    Time = time() - Time
    # Compute Final Error
    ErrX = norm(X_t-X_true)
    ErrXLam = norm([X_t-X_true,Lam_t-Lam_true])
    ErrK = 0
    a_t = rand(MvNormal(zeros(nx),Sigma),Int(1e5))'
    Exp1, Exp2 = exp.(a_t*X_t), exp.(a_t*X_true)
    g_t = mean(a_t.*(Exp1-Exp2)./((1 .+ Exp1).*(1 .+ Exp2)),dims=1)'[:,1] + Pen*(X_t-X_true)
    G_tTLam_t, c_t = AA'*Lam_t, AA*X_t-BB
    KKT = norm([g_t+G_tTLam_t;c_t])
    # Compute the 95% CI radius over primal individuals
    Omega_t = inv(Matrix(K_t))*hcat(vcat(cum_barggT_t-cum_barg_t*cum_barg_t',zeros(nlam,nx)),zeros(nx+nlam,nlam))*inv(Matrix(K_t))/(2*(c2<1)+(2-1/c1)*(c2==1))
    CompFlops += (nx+nlam)^3
    ## vectorize inference
    Radius = 1.96*sqrt.(Step_size*diag(Omega_t)[1:nx])
    Cov_Id = (X_true.>= X_t-Radius).*(X_true.<=X_t+Radius)
    ## Test Inf
    COV_value = sum(Omega_t)/(nx+nlam)^2
    Radius_XLam = 1.96*sqrt(Step_size)*sqrt(COV_value)
    Cov_Id_XLam = mean([X_true;Lam_true])>=(mean([X_t;Lam_t])-Radius_XLam) && mean([X_true;Lam_true])<=(mean([X_t;Lam_t])+Radius_XLam)

    return ErrX,ErrXLam,ErrK,KKT,Time,CompFlops/t,Radius,Cov_Id,0,Cov_Id_XLam
end
