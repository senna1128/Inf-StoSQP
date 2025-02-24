%% This fucntion draws plot for simulation
clear all; close all
delete('./Figure/*')
%rmdir('./Figure/*','s')

%% Extract parameter
Rep=200; tau=[0,20,40,60]; c_1=1; c_2=[0.501]; c_3=2; R_Tail=1e5+1;
SigToe=[0,0.4,0.5,0.6]; SigEqui=[0.1,0.2,0.3];
D=[5,20,40,60]; M=[0.5];
LenTau=length(tau); LenStep=length(c_2); LenDim=4;
LenM=length(M); LenSigToe=length(SigToe); LenSigEqui=length(SigEqui);

%% Save Results to txt file
for IdDim = 1:LenDim
    for IdStep = 1:LenStep
        for IdM = 1: LenM
            file_path = ['./Figure/SummaryS' num2str(IdStep)...
                'D' num2str(IdDim) 'M' num2str(IdM) '.txt'];
            fileID = fopen(file_path,'w');
            fprintf(fileID,'%11s %19s %30s %27s %20s %20s\r\n','Tau','MSE (std)', ...
                'Rate (std)','Length (std)','FLOPs');
            % load results summary
            load(['./Solution/S',num2str(IdStep),'D',num2str(IdDim), ...
                'M',num2str(IdM),'/Summary.mat'])
            % Save Results
            for IdVar = 1:LenSigToe+LenSigEqui
                if IdVar == 1
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'\n%1s','Iden');
                elseif 2<= IdVar & IdVar<=LenSigToe
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'\n%1s',['Top',num2str(SigToe(IdVar))]);
                else
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'\n%1s',['Equ',num2str(SigEqui(IdVar-LenSigToe))]);
                end
                for IdMethod = 1:LenTau+1
                    if IdMethod == 1
                        fileID = fopen(file_path,'a');
                        fprintf(fileID,'%7s','ConstM');
                    else
                        fileID = fopen(file_path,'a');
                        fprintf(fileID,'%10s',num2str(tau(IdMethod-1)));
                    end
                    % Save MSE
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'%16.4e(%0.4e)', MSE_X(IdVar,IdMethod,1), ...
                        MSE_X(IdVar,IdMethod,2));
                    % Save Coverage Rate
                    Rate = mean(Cov_Id_Entry{IdVar,IdMethod}(:,1));
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'%18.4f(%0.4f)',Rate,sqrt(Rate*(1-Rate)));
                    % Save Radius
                    Radius_mean = mean(Radius_Entry{IdVar,IdMethod}(:,1));
                    diff_mean = Radius_Entry{IdVar,IdMethod}(:,1) - Radius_mean;
                    Radius_std = sqrt(mean(Radius_Entry{IdVar,IdMethod}(:,2).^2)+mean(diff_mean.^2));
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'%20.4f(%0.4f)',Radius_mean,Radius_std);
                    % Save FLOPs
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'%22.4f',Ave_Flops(IdVar,IdMethod,1));
                end
            end
        end
    end
end


