%% This fucntion draws plot for simulation
clear all; close all
DBeta = [0.501];
Sigma = [1e-4, 1e-2, 1e-1, 1];
Rep = 200;
IdPProb = [1,2,3,4,5];

%% Summarize the results
for ii = 1:length(IdPProb)
    Idprob = IdPProb(ii);
    for IdDBeta = 1:length(DBeta)
        Id_comp = ones(length(Sigma),1);
        for IdSigma = 1:length(Sigma)
            load_file = ['./Solution/StoSQP',num2str(Idprob),'/D',num2str(IdDBeta),'S',num2str(IdSigma),'/rep1.mat'];
            if isfile(load_file)
                load(load_file);
                Res{IdSigma} = Result;
            else
                Id_comp(IdSigma) = 0;
            end
        end
        if sum(Id_comp) < length(Sigma)
            continue
        end
        % Make Output Directory
%        mkdir(['./Figure/Prop' num2str(Idprob)]) 
        Sigma_Map = jet(length(Sigma)+1);
        Sigma_Map = Sigma_Map(2:end,:);        

        % Show Convergence of Iterate
        L_Tail = 1; R_Tail = 1e5+1; t = L_Tail:R_Tail;
        Sigma_Legend = {'$\sigma^2=10^{-4}$','$\sigma^2=10^{-2}$', ...
            '$\sigma^2=10^{-1}$','$\sigma^2=1$','Theory'};
        fig = figure(IdDBeta);
        for IdSigma = 1:length(Sigma)
            plot(L_Tail:R_Tail,[Res{IdSigma}.ErrXLam{L_Tail:R_Tail}], ...
                'Color',Sigma_Map(IdSigma,:),'LineWidth',1.5)
            hold on
        end
        y = sqrt(log(t)./t.^DBeta(IdDBeta));
        plot(L_Tail:R_Tail,y,'Color','red','LineWidth',1.5)
        hold off
        set(gca, 'YScale', 'log')
        set(gca,'fontsize',18)
        xlabel('Iteration $t$','Interpreter','latex','Fontsize',35)
        ylabel('$\|(x_t-x^\star,\lambda_t-\lambda^\star)\|$','Interpreter','latex','Fontsize',35)
        xlim([L_Tail R_Tail])

        % Show Convergence of KKT Matrix
        fig = figure(10 + IdDBeta);
        for IdSigma = 1:length(Sigma)
            plot(L_Tail:R_Tail,[Res{IdSigma}.ErrK{L_Tail:R_Tail}], ...
                'Color',Sigma_Map(IdSigma,:),'LineWidth',1.5)
            hold on
        end
        y = sqrt(log(t)./t.^DBeta(IdDBeta));
        plot(L_Tail:R_Tail,y,'Color','red','LineWidth',1.5)
        hold off
        set(gca, 'YScale', 'log')
        set(gca,'fontsize',18)
        xlabel('Iteration $t$','Interpreter','latex','Fontsize',35) 
        ylabel('$\|K_t-K^\star\|$','Interpreter','latex','Fontsize',35)
        xlim([L_Tail R_Tail])
        
        % Show Convergence of KKT Residual
        fig = figure(100 + IdDBeta);
        for IdSigma = 1:length(Sigma)
            plot(L_Tail:R_Tail,[Res{IdSigma}.KKT{L_Tail:R_Tail}], ...
                'Color',Sigma_Map(IdSigma,:),'LineWidth',1.5)
            hold on
        end
        y = sqrt(log(t)./t.^DBeta(IdDBeta));
        plot(L_Tail:R_Tail,y,'Color','red','LineWidth',1.5)
        hold off
        set(gca, 'YScale', 'log')
        set(gca,'fontsize',18)
        xlabel('Iteration $t$','Interpreter','latex','Fontsize',35)
        ylabel('$\|\nabla L_t\|$','Interpreter','latex','Fontsize',35)
        xlim([L_Tail R_Tail])

        % Show Confidence Interval of X1+Lam1
        t = 1:R_Tail; L_Tail = R_Tail-100;
        for IdSigma = 1:length(Sigma)
            fig = figure(200+IdSigma);
            Beta_t = 1./t.^(DBeta(IdDBeta));
            XLam1 = [Res{IdSigma}.X1{:}] + [Res{IdSigma}.Lam1{:}];
            Rad = [Res{IdSigma}.Radius{:}];
            XLam1_true = Res{IdSigma}.X_true(1)+ Res{IdSigma}.Lam_true(1);
            plot(L_Tail:R_Tail,XLam1(L_Tail:R_Tail)+Rad(L_Tail:R_Tail), 'Color','green', 'LineWidth', 1.5);
            hold on;
            plot(L_Tail:R_Tail,XLam1(L_Tail:R_Tail)-Rad(L_Tail:R_Tail),'Color','blue', 'LineWidth', 1.5);
            hold on
            yline(XLam1_true,'LineWidth',2,'Color','red')
            hold off
            xlim([L_Tail,R_Tail])
            ylabel('$95\%$ CI of $x_1^\star+\lambda_1^\star$','Interpreter','latex','FontSize',35)
            xlabel('Iteration $t$','Interpreter','latex','FontSize',35)
        end

        % save cover rate to the file
        file_path_1 = './Figure/CoverRate.txt';
        file_path_2 = './Figure/Radius.txt';
        rate_cover = ones(length(Sigma),1);
        radius = ones(length(Sigma),2);
        for IdSigma = 1:length(Sigma)
            folder_path = ['./Solution/StoSQP',num2str(Idprob),'/D',num2str(IdDBeta),'S',num2str(IdSigma)];
            directory = dir([folder_path,'/*.mat']);
            NumCov = numel(directory);
            rrate = zeros(NumCov,1);
            rradius = zeros(NumCov,1);
            iii = 1; 
            for dat = directory'
                load([folder_path,'/',dat.name]);
                rrate(iii) = double(Result.CoverRate);
                rradius(iii) = Result.Radius{end};
                iii = iii + 1;
            end
            rate_cover(IdSigma,1) = mean(rrate);
            radius(IdSigma,1) = mean(rradius);
            radius(IdSigma,2) = std(rradius);
        end

        if isfile(file_path_1)
            fileID = fopen(file_path_1,'a');
            fprintf(fileID,'%1.0f %8.4f %11.4f %14.4f %17.4f\r\n',Idprob, ...
                rate_cover(1,1),rate_cover(2,1),rate_cover(3,1),rate_cover(4,1));
        else
            fileID = fopen(file_path_1,'w');
            fprintf(fileID,'%1.0f %8.4f %11.4f %14.4f %17.4f\r\n',Idprob, ...
                rate_cover(1,1),rate_cover(2,1),rate_cover(3,1),rate_cover(4,1));
        end

        if isfile(file_path_2)
            fileID = fopen(file_path_2,'a');
            fprintf(fileID,'%1.0f %8.4f %11.4f %14.4f %17.4f\r\n',Idprob, ...
                radius(1,1),radius(2,1),radius(3,1),radius(4,1));
            fprintf(fileID,'%1.0f %8.4f %11.4f %14.4f %17.4f\r\n',Idprob, ...
                radius(1,2),radius(2,2),radius(3,2),radius(4,2));
        else
            fileID = fopen(file_path_2,'w');
            fprintf(fileID,'%1.0f %8.4f %11.4f %14.4f %17.4f\r\n',Idprob, ...
                radius(1,1),radius(2,1),radius(3,1),radius(4,1));
            fprintf(fileID,'%1.0f %8.4f %11.4f %14.4f %17.4f\r\n',Idprob, ...
                radius(1,2),radius(2,2),radius(3,2),radius(4,2));
        end
  
    end
end






