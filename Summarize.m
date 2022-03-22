%% This fucntion draws plot for simulation
clear all; close all
%delete('./Figure/*')
%rmdir('./Figure/*','s')

DBeta = [0.5,0.6,0.7,0.8,0.9,1.0];
Sigma = [1e-8, 1e-4, 1e-2, 1e-1, 1];

%% Summarize the results
% Idprob = [2,32,22]

load(['./Solution/StoSQP', num2str(Idprob),'.mat'])
if length(Result.Time)==1 || sum(cellfun(@isempty,Result.Time),'all') >= 1
    continue
end
% Make Output Directory
mkdir(['./Figure/Prop' num2str(Idprob)])
%% Show Convergence Plot
% Show Convergence of Iterate
L_Tail = 1;
R_Tail = 1e5+1;
t = L_Tail:R_Tail;
Sigma_Map = jet(5);
Sigma_Legend = {'$\sigma^2=10^{-8}$','$\sigma^2=10^{-4}$','$\sigma^2=10^{-2}$', ...
        '$\sigma^2=10^{-1}$','$\sigma^2=1$','Theory'};

for step = 1:6
    fig = figure(step);
    for sigma = 1:5
        plot(L_Tail:R_Tail,[Result.ErrXLam{step,sigma}{1}{L_Tail:R_Tail}],'Color',Sigma_Map(sigma,:),'LineWidth',1.5)
        hold on
    end
    y = sqrt(log(t)./t.^DBeta(step));
    plot(L_Tail:R_Tail,y,'Color','red','LineWidth',1.5)
    hold off
    set(gca, 'YScale', 'log')
    set(gca,'fontsize',18)
    xlabel('Iteration $t$','Interpreter','latex','Fontsize',35) 
    ylabel('$\|(x_t-x^\star,\lambda_t-\lambda^\star)\|$','Interpreter','latex','Fontsize',35)         
    xlim([L_Tail R_Tail])
%    ylim([1e-2 1e2])
%    legend((Sigma_Legend),'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',7)
    filename = ['./Figure/Prop' num2str(Idprob) '/ErrXLamS' num2str(step) '.png'];
    print('-dpng', filename)
end
% Show Convergence of KKT Matrix
for step = 1:6
    fig = figure(step);
    for sigma = 1:5
        plot(L_Tail:R_Tail,[Result.ErrK{step,sigma}{1}{L_Tail:R_Tail}],'Color',Sigma_Map(sigma,:),'LineWidth',1.5)
        hold on
    end
    y = sqrt(log(t)./t.^DBeta(step));
    plot(L_Tail:R_Tail,y,'Color','red','LineWidth',1.5)
    hold off
    set(gca, 'YScale', 'log')
    set(gca,'fontsize',18)
    xlabel('Iteration $t$','Interpreter','latex','Fontsize',35) 
    ylabel('$\|K_t-K^\star\|$','Interpreter','latex','Fontsize',35)         
    xlim([L_Tail R_Tail])
%        ylim([1e-2 1e2])
%        legend((Sigma_Legend),'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',7)
    filename = ['./Figure/Prop' num2str(Idprob) '/ErrKS' num2str(step) '.png'];
    print('-dpng', filename)
end
% Show Convergence of KKT Residual
for step = 1:6
    fig = figure(step);
    for sigma = 1:5
        plot(L_Tail:R_Tail,[Result.KKT{step,sigma}{1}{L_Tail:R_Tail}],'Color',Sigma_Map(sigma,:),'LineWidth',1.5)
        hold on
    end
    y = sqrt(log(t)./t.^DBeta(step));
    plot(L_Tail:R_Tail,y,'Color','red','LineWidth',1.5)
    hold off
    set(gca, 'YScale', 'log')
    set(gca,'fontsize',18)
    xlabel('Iteration $t$','Interpreter','latex','Fontsize',35) 
    ylabel('$\|\nabla L_t\|$','Interpreter','latex','Fontsize',35)         
    xlim([L_Tail R_Tail])
%        ylim([1e-2 1e2])
%        legend((Sigma_Legend),'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',7)
    filename = ['./Figure/Prop' num2str(Idprob) '/ErrKKTS' num2str(step) '.png'];
    print('-dpng', filename)
end
    
%% Asymptotic normality of X1 and Lam1
% Asymptotic normality of X1    
t = 1:R_Tail;
L_Tail = 0.5e5;
for step = 1:6
    for sigma = 1:5
        fig = figure(step*10+sigma);
        Beta_t = 2./t.^(DBeta(step));
        X1Err = sqrt(1./Beta_t).*([Result.X1{step,sigma}{1}{:}] - Result.X_true(1));
        histogram(X1Err(L_Tail:end),'Normalization','pdf')
        hold on
        pd = fitdist(X1Err(L_Tail:end)','normal');
        x_pdf = linspace(min(X1Err(L_Tail:end)),max(X1Err(L_Tail:end)));
        y_pdf = pdf(pd,x_pdf);
        line(x_pdf,y_pdf,'LineWidth',2,'Color','red')
        hold off
        ylabel('Probability','FontSize',45)         
        xlabel('$[x_t]_1-x_1^\star$','Interpreter','latex','FontSize',55)
        filename = ['./Figure/Prop' num2str(Idprob) '/HistX1S' num2str(step) 'Sig' num2str(sigma) '.png'];
        print('-dpng', filename)
    end
end
% Asymptotic normality of Lam1    
t = 1:R_Tail;
L_Tail = 0.5e5;
for step = 1:6
    for sigma = 1:5
        fig = figure(step*10+sigma);
        Beta_t = 2./t.^(DBeta(step));
        Lam1Err = sqrt(1./Beta_t).*([Result.Lam1{step,sigma}{1}{:}] - Result.Lam_true(1));
        histogram(Lam1Err(L_Tail:end),'Normalization','pdf')
        hold on
        pd = fitdist(Lam1Err(L_Tail:end)','normal');
        x_pdf = linspace(min(Lam1Err(L_Tail:end)),max(Lam1Err(L_Tail:end)));
        y_pdf = pdf(pd,x_pdf);
        line(x_pdf,y_pdf,'LineWidth',2,'Color','red')
        hold off
        ylabel('Probability','FontSize',45)         
        xlabel('$[\lambda_t]_1-\lambda_1^\star$','Interpreter','latex','FontSize',55)
        filename = ['./Figure/Prop' num2str(Idprob) '/HistLam1S' num2str(step) 'Sig' num2str(sigma) '.png'];
        print('-dpng', filename)
    end
end
%% Show Confidence Interval of X1+Lam1
t = 1:R_Tail;
L_Tail = R_Tail-100;
for step = 1:6
    for sigma = 1:5
        fig = figure(step*10+sigma);
        Beta_t = 2./t.^(DBeta(step));
        XLam1 = [Result.X1{step,sigma}{1}{:}] + [Result.Lam1{step,sigma}{1}{:}];
        QA = [Result.Radius{step,sigma}{1}{:}];
        Rad = 1.96*sqrt(Beta_t).*sqrt(QA(1,:)+2*QA(2,:)+QA(3,:));
        XLam1_true = Result.X_true(1)+ Result.Lam_true(1);
        plot(L_Tail:R_Tail,XLam1(L_Tail:end)+Rad(L_Tail:end), 'Color','green', 'LineWidth', 1.5);
        hold on;
        plot(L_Tail:R_Tail,XLam1(L_Tail:end)-Rad(L_Tail:end),'Color','blue', 'LineWidth', 1.5);
        hold on
        yline(XLam1_true,'LineWidth',2,'Color','red')
        hold off
%            ylim([min([XLam1(L_Tail:end)-Rad(L_Tail:end) XLam1_true]-0.1), max([XLam1(L_Tail:end)+Rad(L_Tail:end) XLam1_true]+0.1)])
        xlim([L_Tail,R_Tail])
        ylabel('$95\%$ CI of $x_1^\star+\lambda_1^\star$','Interpreter','latex','FontSize',35)
        xlabel('Iteration $t$','Interpreter','latex','FontSize',35)
        filename = ['./Figure/Prop' num2str(Idprob) '/ConfS' num2str(step) 'Sig' num2str(sigma) '.png'];
        print('-dpng', filename)            
    end
end




