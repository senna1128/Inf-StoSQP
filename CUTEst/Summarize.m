%% This fucntion draws plot for simulation
clear all; close all
delete('./Figure/*')
rmdir('./Figure/*','s')

DBeta = [0.501];
Sigma = [1e-4, 1e-2, 1e-1, 1];
Rep = 200;
%% Go over Problem Folder
d = dir('./Solution');
isub = [d(:).isdir]; 
nameFolds = {d(isub).name}; % Names of all directories
nameFolds = nameFolds(~ismember(nameFolds, {'.', '..'}))';
numFolds = sort(cellfun(@(fold) str2double(fold(5:end)),nameFolds));
Problems = readlines('./Parameter/problems.txt');


%% Plot the results
for ii = 1:length(numFolds)
    Idprob = numFolds(ii);
    for IdDBeta = 1:length(DBeta)
        Id_comp = ones(length(Sigma),1);
        for IdSigma = 1:length(Sigma)
            load_file = ['./Solution/Prob',num2str(Idprob),'/D',num2str(IdDBeta),'S',num2str(IdSigma),'/rep2.mat'];
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
        %% Make Output Directory
        mkdir(['./Figure/Prop' num2str(Idprob) '_']+Problems(Idprob))
        %% Show Convergence of Iterate
        L_Tail = 1; R_Tail = 1e5+1; t = L_Tail:R_Tail;
        Sigma_Map = jet(length(Sigma)+1);
        Sigma_Legend = {'$\sigma^2=10^{-4}$','$\sigma^2=10^{-2}$', ...
            '$\sigma^2=10^{-1}$','$\sigma^2=1$','Theory'};
        fig = figure(IdDBeta);
        for IdSigma = 1:length(Sigma)
            plot(L_Tail:R_Tail,[Res{IdSigma}.ErrXLam{L_Tail:R_Tail}], ...
                'Color',Sigma_Map(IdSigma+1,:),'LineWidth',1.5)
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
%        ylim([1e-2 1e2])
%        legend((Sigma_Legend),'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',29)
        filename = ['./Figure/Prop' num2str(Idprob) '_']+Problems(Idprob)+['/ErrXLamS' num2str(IdDBeta) '.png'];
        print('-dpng', filename)

        %% Show Convergence of KKT Matrix
        fig = figure(10 + IdDBeta);
        for IdSigma = 1:length(Sigma)
            plot(L_Tail:R_Tail,[Res{IdSigma}.ErrK{L_Tail:R_Tail}], ...
                'Color',Sigma_Map(IdSigma+1,:),'LineWidth',1.5)
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
%        ylim([1e-2 1e2])
%        legend((Sigma_Legend),'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',7)
        filename = ['./Figure/Prop' num2str(Idprob) '_']+Problems(Idprob)+['/ErrKS' num2str(IdDBeta) '.png'];
        print('-dpng', filename)

        %% Show Convergence of KKT Residual
        fig = figure(100 + IdDBeta);
        for IdSigma = 1:length(Sigma)
            plot(L_Tail:R_Tail,[Res{IdSigma}.KKT{L_Tail:R_Tail}], ...
                'Color',Sigma_Map(IdSigma+1,:),'LineWidth',1.5)
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
%        ylim([1e-2 1e2])
%        legend((Sigma_Legend),'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',7)
        filename = ['./Figure/Prop' num2str(Idprob) '_']+Problems(Idprob)+['/ErrKKTS' num2str(IdDBeta) '.png'];
        print('-dpng', filename)
    end
end

%% Save results to txt file
clear all; close all
delete('./Figure/SummaryS*')

DBeta = [0.501];
Sigma = [1e-4, 1e-2, 1e-1, 1];
Rep = 200;
%% Go over Problem Folder
d = dir('./Solution');
isub = [d(:).isdir]; 
nameFolds = {d(isub).name}; % Names of all directories
nameFolds = nameFolds(~ismember(nameFolds, {'.', '..'}))';
numFolds = sort(cellfun(@(fold) str2double(fold(5:end)),nameFolds));
Problems = readlines('./Parameter/problems.txt');

for IdDBeta = 1:length(DBeta)
    file_path = ['./Figure/SummaryS' num2str(IdDBeta) '.txt'];
    fileID = fopen(file_path,'w');
    fprintf(fileID,'%17s %23s %31s %27s %20s\r\n','Sigma','MSE (std)', ...
        'Rate (std)','Length (std)','FLOPs');
    %% Go over problems
    for ii = 1:length(numFolds)
        Idprob = numFolds(ii);
        load(['./Solution/Prob',num2str(Idprob),'/SaveSolution.mat'])
        if rank(full(G_star*G_star')) == size(G_star,1)
            % Index for inference
            IdIndex = find(vecnorm(G_star'*inv(G_star*G_star')*G_star*eye(size(G_star,2))- eye(size(G_star,2)),2)>=0.05);
            fileID = fopen(file_path,'a');
            fprintf(fileID,'\n%1s',Problems(Idprob));            
            % load summary results
            load(['./Solution/Prob',num2str(Idprob),'/Summary.mat'])
            % Go over Sigma
            for IdSigma = 1:length(Sigma)
                % Save MSE 
                fileID = fopen(file_path,'a');
                fprintf(fileID,'%12.4e %16.4e(%0.4e)', Sigma(IdSigma), ...
                    MSE_X(IdDBeta,IdSigma,1),MSE_X(IdDBeta,IdSigma,2));
                % Save Coverage Rate
                Rate = mean(Cov_Id_Entry{IdDBeta,IdSigma}(IdIndex,1));
                fileID = fopen(file_path,'a');
                fprintf(fileID,'%18.4f(%0.4f)',Rate,sqrt(Rate*(1-Rate)));
                % Save Radius
                Radius_mean = mean(Radius_Entry{IdDBeta,IdSigma}(IdIndex,1));
                diff_mean = Radius_Entry{IdDBeta,IdSigma}(IdIndex,1) - Radius_mean;
                Radius_std = sqrt(mean(Radius_Entry{IdDBeta,IdSigma}(IdIndex,2).^2)+mean(diff_mean.^2));
                fileID = fopen(file_path,'a');
                fprintf(fileID,'%20.4f(%0.4f)',Radius_mean,Radius_std);
                % Save FLOPs
                fileID = fopen(file_path,'a');
                fprintf(fileID,'%22.4f\r\n',Ave_Flops(IdDBeta,IdSigma,1));
            end
        end
    end
end







