# Impmement Stochastic SQP framework
# nlp: model
# c1,c2,c3: step decay ratio
# sigma2: estimate variance
# Max_Iter: maximum iterations
# tau: number of inner iterations of the solver
# X_true, Lam_true: optimal solution
# Output:
# X1, Lam1, Radius, ErrXLam, ErrK, KKT, Time,IdCov,IdSing


function StoSQP(nlp,c1,c2,c3,sigma2,Max_Iter,tau,X_true,Lam_true,EPS=1e-8)
    Buff= Int(floor(Max_Iter/1.1))
    # Problem dimension
    nx, nlam = nlp.meta.nvar, nlp.meta.ncon
    # Compute K^star
    G_star, B_star = jac(nlp,X_true), hess(nlp,X_true,Lam_true)
    K_star = hcat(vcat(B_star,G_star),vcat(G_star',zeros(nlam, nlam)))
    # Initialize variables
    t, CompFlops, X_t, Lam_t, Step_size = 0, 0, nlp.meta.x0, nlp.meta.y0, 0
    cum_barg_t, NewDir_t, K_t = zeros(nx), zeros(nx+nlam), zeros(nx+nlam, nx+nlam)
    cum_barggT_t, cum_bar2x_f_t = zeros(nx,nx), Matrix(I,nx,nx)
    CovM = sigma2*(Matrix(I,nx,nx)+ones(nx,nx))
    ErrX, ErrXLam, ErrK, KKT = [], [], [], []
    # Start the loop
    Time = time()
    while t <= Max_Iter
        # Step 0: compute deterministic quantities
        push!(ErrX, norm(X_t-X_true))
        push!(ErrXLam, norm([X_t-X_true;Lam_t-Lam_true]))
        # Constraints value and Jacobian
        c_t, G_t = consjac(nlp, X_t)
        cum_bar2x_L_t = cum_bar2x_f_t + hess(nlp,X_t,Lam_t,obj_weight=0.0)
        CompFlops += (nx+nlam)^2
        # Compute the reduced Hessian and do Hessian modification
        Q_t, R_t = qr(Matrix(G_t'))
        if abs(R_t[end,end])< EPS || ErrXLam[end]>1e10
            return [],[],[],[],0,0,[],[],1,0
        end
        delta_tt = eigmin(Hermitian(Q_t[:,nlam+1:end]'cum_bar2x_L_t*Q_t[:,nlam+1:end],:L))
        if nlam<nx && delta_tt <= EPS
            delta_t = abs(delta_tt) + 0.1
        else
            delta_t = 0
        end
        K_t = hcat(vcat(cum_bar2x_L_t+delta_t*Matrix(I,nx,nx),G_t),vcat(G_t',zeros(nlam,nlam)))
        push!(ErrK, opnorm(K_t-K_star,1))

        # Step 1: generate random gradient and Hessian
        g_t,G_tTLam_t,nab_x2f_t = grad(nlp,X_t),jtprod(nlp,X_t,Lam_t),hess(nlp,X_t,zeros(nlam))
        barg_t = rand(MvNormal(g_t,CovM))
        if t > Buff
            cum_barggT_t = (t-Buff)/(t-Buff+1)*cum_barggT_t + 1/(t-Buff+1)*barg_t*barg_t'
            cum_barg_t = (t-Buff)/(t-Buff+1)*cum_barg_t + 1/(t-Buff+1)*barg_t
        end
        push!(KKT,norm([g_t+G_tTLam_t; c_t]))
        # define stochastic Lagrangian gradient
        bar_nab_L_t = [barg_t+G_tTLam_t; c_t]
        Rand_t = rand(Normal(0,sigma2^(1/2)),nx,nx)
        bar_nab_x2f_t = Hermitian(nab_x2f_t,:L) + (Rand_t+Rand_t')/2
        cum_bar2x_f_t = t/(t+1)*cum_bar2x_f_t + 1/(t+1)*bar_nab_x2f_t
        CompFlops += 2*nx^2
        # Step 2: solve Newton system via Randomized Solve
        NewDir_t = zeros(nx+nlam)
        for inner_iter = 1:tau
            # random pick a index
            j = sample(1:(nx+nlam))
            NewDir_t = NewDir_t - (K_t[j,:]'*NewDir_t+bar_nab_L_t[j])/norm(K_t[j,:])^2*K_t[:,j]
        end
#        NewDir_t = lu(K_t)\-bar_nab_L_t
        CompFlops += tau*(nx+nlam)
        # Step 3: update the iterate
        beta_t = c1/(t+1)^c2
        chi_t = beta_t^c3
        Step_size = rand(Uniform(beta_t,beta_t+chi_t))
        X_t = X_t + Step_size*NewDir_t[1:nx]
        Lam_t = Lam_t + Step_size*NewDir_t[nx+1:end]
        t = t + 1
    end
    Time = time() - Time
    # Compute the 95% CI radius over primal individuals
    Omega_t = inv(Matrix(K_t))*hcat(vcat(cum_barggT_t-cum_barg_t*cum_barg_t',zeros(nlam,nx)),zeros(nx+nlam,nlam))*inv(Matrix(K_t))/(2*(c2<1)+(2-1/c1)*(c2==1))
    CompFlops += (nx+nlam)^3
    ## vectorize inference
    Radius = 1.96*sqrt.(Step_size*diag(Omega_t)[1:nx])
    Cov_Id = (X_true.>= X_t-Radius).*(X_true.<=X_t+Radius)

    return ErrX,ErrXLam,ErrK,KKT,Time,CompFlops/Max_Iter,Radius,Cov_Id,0

end
