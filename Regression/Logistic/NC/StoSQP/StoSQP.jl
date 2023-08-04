# Impmement Stochastic SQP framework
# c1,c2,c3: step decay ratio
# Max_Iter: maximum iterations
# ttau: number of inner iterations of the solver
# nx,nlam: dimension of variable and constraint
# X_true, AA: optimal solution, constraints
# RR, IdR: covariance setup and indicator of covaraince type
# Output:
# X1, Lam1, Radius, ErrXLam, ErrK, KKT, Time, IdCoverage, IdSing

function StoSQP(c1,c2,c3,Max_Iter,ttau,nx,nlam,X_true,AA,RR,IdR,EPS=1e-8,Pen=1)
    ## Set Buffer Size
    Buff= Int(floor(Max_Iter/1.5))
    ## Define Covariance Matrix
    if IdR == 1
        Sigma = [RR^abs(i-j) for i=1:nx, j=1:nx]
    else
        Sigma = RR*ones(nx,nx)+(1-RR)*Matrix(1.0I,nx,nx)
    end
    ## True value
    Lam_true = zeros(nlam)

    ## Initialize variables
    t,X_t,Lam_t,NewDir_t = 1,ones(nx),ones(nlam),zeros(nx+nlam)
    cum_barg_t,cum_barggT_t,cum_bar2x_f_t = zeros(nx),zeros(nx,nx),Matrix(1.0I,nx,nx)
    X1, Lam1, Radius, ErrXLam, ErrK, KKT = [], [], [], [], [], []
    ## Start the loop
    Time = time()
    while t <= Max_Iter+1
        # Step 0: compute deterministic quantities
        push!(X1, mean(X_t)), push!(Lam1, mean(Lam_t))
        push!(ErrXLam, norm([X_t-X_true;Lam_t-Lam_true]))
        # Constraints value and Jacobian
        c_t = 0.5*norm(X_t)^2 - AA
        cum_bar2x_L_t = cum_bar2x_f_t + Pen*Matrix(1.0I,nx,nx)
        K_t = hcat(vcat(cum_bar2x_L_t,X_t'),vcat(X_t,zeros(nlam,nlam)))
        push!(ErrK, 0)
        # Compute the 95% CI radius
        if t > Buff
            Omega_t = inv(Matrix(K_t))*hcat(vcat(cum_barggT_t-cum_barg_t*cum_barg_t',zeros(nlam,nx)),zeros(nx+nlam,nlam))*inv(Matrix(K_t))
        else
            Omega_t = zeros(nx+nlam,nx+nlam)
        end
        beta_t = c1/t^c2
        chi_t = beta_t^c3
        COV_value = sum(Omega_t)/(nx+nlam)^2
        if COV_value < 0
            return [],[],[],[],[],[],0,0,1
        end
        radius = 1.96*sqrt(beta_t)*sqrt(COV_value/(2*(c2<1)+(2-1/c1)*(c2==1)))
        push!(Radius,radius)
        # Step 1: generate random gradient and Hessian
        a_t = rand(MvNormal(zeros(nx),Sigma))
        b_t = 2*rand(Bernoulli(1/(1+exp(-a_t'X_t))))-1
        barg_t = -b_t/(1+exp(b_t*a_t'X_t))*a_t
        G_tTLam_t = X_t.*Lam_t
        if t > Buff
            cum_barggT_t = (t-Buff)/(t-Buff+1)*cum_barggT_t + 1/(t-Buff+1)*barg_t*barg_t'
            cum_barg_t = (t-Buff)/(t-Buff+1)*cum_barg_t + 1/(t-Buff+1)*barg_t
        end
#        cum_barggT_t = t^log(t+1)/(t+1)^log(t+2)*cum_barggT_t + (1-t^log(t+1)/(t+1)^log(t+2))*barg_t*barg_t'
#        cum_barg_t = t^log(t+1)/(t+1)^log(t+2)*cum_barg_t + (1-t^log(t+1)/(t+1)^log(t+2))*barg_t
        bar_nab_L_t = vcat(barg_t+Pen*(X_t-X_true)+G_tTLam_t,c_t)
        push!(KKT,norm(bar_nab_L_t))
        # define stochastic Lagrangian gradient
        bar_nab_x2f_t = a_t*a_t'/((1+exp(a_t'X_t))*(1+exp(-a_t'X_t)))
        cum_bar2x_f_t = t/(t+1)*cum_bar2x_f_t + 1/(t+1)*bar_nab_x2f_t
        # Step 2: solve Newton system via Randomized Solve
        if ttau == 0
            NewDir_t = lu(K_t)\-bar_nab_L_t
        else
            NewDir_t = zeros(nx+nlam)
            for inner_iter = 1:ttau
                # random pick a index
                j = sample(1:(nx+nlam))
                NewDir_t = NewDir_t - (K_t[j,:]'*NewDir_t+bar_nab_L_t[j])/(K_t^2)[j,j]*K_t[:,j]
            end
        end
        # Step 3: update the iterate
        Step_size = rand(Uniform(beta_t,beta_t+chi_t))
        X_t = X_t + Step_size*NewDir_t[1:nx]
        Lam_t = Lam_t + Step_size*NewDir_t[nx+1:end]
        t = t + 1
    end
    Time = time() - Time
    IdCov = (sum(X_true)+sum(Lam_true))/(nx+nlam)>=((X1[end]*nx+Lam1[end]*nlam)/(nx+nlam)-Radius[end]) && (sum(X_true)+sum(Lam_true))/(nx+nlam)<=((X1[end]*nx+Lam1[end]*nlam)/(nx+nlam)+Radius[end])
    return X1,Lam1,Radius,ErrXLam,ErrK,KKT,Time,IdCov,0
end
