%% This fucntion draws plot for simulation
clear all; close all
%delete('./Figure/*')
%rmdir('./Figure/*','s')

%% Extract parameter
R_Tail=1e5+1; Rep=200; tau=[0,20,40,60]; c_1=1; c_2=[0.501]; c_3=2;
SigToe=[0,0.4,0.5,0.6]; SigEqui=[0.1,0.2,0.3];
D=[5,20,40,60]; M=[0.5];
LenTau=length(tau); LenStep=length(c_2); LenDim=length(D);
LenM=length(M); LenSigToe=length(SigToe); LenSigEqui=length(SigEqui);

%% Save Coverage Rate and Length
Rate_Cov_T = zeros(LenStep,LenDim,LenM,LenTau,LenSigToe,Rep);
Rate_Cov_E = zeros(LenStep,LenDim,LenM,LenTau,LenSigEqui,Rep);
Radius_T = zeros(LenStep,LenDim,LenM,LenTau,LenSigToe,Rep);
Radius_E = zeros(LenStep,LenDim,LenM,LenTau,LenSigEqui,Rep);

for IdStep = 1:LenStep 
    for IdDim = 1:LenDim
        for IdM = 1: LenM
            for IdTau = 1:LenTau
                load_file = ['Step',num2str(IdStep),'Dim',num2str(IdDim),'M',num2str(IdM),'Tau',num2str(IdTau)];
                load(['./Solution/',load_file,'/ToeCR.mat']);
                load(['./Solution/',load_file,'/EquCR.mat']);
                for IdSigToe = 1:LenSigToe
                    for IdRep = 1:Rep
                        Rate_Cov_T(IdStep,IdDim,IdM,IdTau,IdSigToe,IdRep) = double(CountConvT{IdSigToe}{IdRep});
                        Radius_T(IdStep,IdDim,IdM,IdTau,IdSigToe,IdRep) = CountRadiusT{IdSigToe}{IdRep};
                    end
                end
                for IdSigEqui = 1:LenSigEqui
                    for IdRep = 1:Rep
                        Rate_Cov_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,IdRep) = double(CountConvE{IdSigEqui}{IdRep});
                        Radius_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,IdRep) = CountRadiusE{IdSigEqui}{IdRep};
                    end
                end
            end
        end
    end
end

file_path = './Figure/CoverRate.txt';
for IdStep = 1:LenStep 
    for IdDim = 1:LenDim
        for IdM = 1: LenM
            for IdTau = 1:LenTau
                if isfile(file_path)
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'\r\n\r\n Step %1s D %3s M %5s Tau %7s\r\n',num2str(c_2(IdStep)),num2str(D(IdDim)),num2str(M(IdM)),num2str(tau(IdTau)));
                    for IdSigToe = 1:LenSigToe
                        Format = [' Toe',num2str(IdSigToe),' %',num2str(IdSigToe*2-1),'.4f %' num2str(IdSigToe*2),'.4f'];
                        fprintf(fileID,Format,mean(Rate_Cov_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)),std(Rate_Cov_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)));
                    end
                    fprintf(fileID,'\r\n');
                    for IdSigEqui = 1:LenSigEqui
                        Format = [' Equ',num2str(IdSigEqui),' %',num2str(IdSigEqui*2-1),'.4f %' num2str(IdSigEqui*2),'.4f'];
                        fprintf(fileID,Format,mean(Rate_Cov_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)),std(Rate_Cov_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)));
                    end                    
                else
                    fileID = fopen(file_path,'w');
                    fprintf(fileID,'\r\n Step %1s D %3s M %5s Tau %7s\r\n',num2str(c_2(IdStep)),num2str(D(IdDim)),num2str(M(IdM)),num2str(tau(IdTau)));
                    for IdSigToe = 1:LenSigToe
                        Format = [' Toe',num2str(IdSigToe),' %',num2str(IdSigToe*2-1),'.4f %' num2str(IdSigToe*2),'.4f'];
                        fprintf(fileID,Format,mean(Rate_Cov_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)),std(Rate_Cov_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)));
                    end
                    fprintf(fileID,'\r\n');
                    for IdSigEqui = 1:LenSigEqui
                        Format = [' Equ',num2str(IdSigEqui),' %',num2str(IdSigEqui*2-1),'.4f %' num2str(IdSigEqui*2),'.4f'];
                        fprintf(fileID,Format,mean(Rate_Cov_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)),std(Rate_Cov_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)));
                    end                    
                end
            end
        end
    end
end

file_path = './Figure/Radius.txt';
for IdStep = 1:LenStep 
    for IdDim = 1:LenDim
        for IdM = 1: LenM
            for IdTau = 1:LenTau
                if isfile(file_path)
                    fileID = fopen(file_path,'a');
                    fprintf(fileID,'\r\n\r\n Step %1s D %3s M %5s Tau %7s\r\n',num2str(c_2(IdStep)),num2str(D(IdDim)),num2str(M(IdM)),num2str(tau(IdTau)));
                    for IdSigToe = 1:LenSigToe
                        Format = [' Toe',num2str(IdSigToe),' %',num2str(IdSigToe*2-1),'.4f %' num2str(IdSigToe*2),'.4f'];
                        fprintf(fileID,Format,mean(Radius_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)*2),std(Radius_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)*2));
                    end
                    fprintf(fileID,'\r\n');
                    for IdSigEqui = 1:LenSigEqui
                        Format = [' Equ',num2str(IdSigEqui),' %',num2str(IdSigEqui*2-1),'.4f %' num2str(IdSigEqui*2),'.4f'];
                        fprintf(fileID,Format,mean(Radius_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)*2),std(Radius_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)*2));
                    end                    
                else
                    fileID = fopen(file_path,'w');
                    fprintf(fileID,'\r\n Step %1s D %3s M %5s Tau %7s\r\n',num2str(c_2(IdStep)),num2str(D(IdDim)),num2str(M(IdM)),num2str(tau(IdTau)));
                    for IdSigToe = 1:LenSigToe
                        Format = [' Toe',num2str(IdSigToe),' %',num2str(IdSigToe*2-1),'.4f %' num2str(IdSigToe*2),'.4f'];
                        fprintf(fileID,Format,mean(Radius_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)*2),std(Radius_T(IdStep,IdDim,IdM,IdTau,IdSigToe,:)*2));
                    end
                    fprintf(fileID,'\r\n');
                    for IdSigEqui = 1:LenSigEqui
                        Format = [' Equ',num2str(IdSigEqui),' %',num2str(IdSigEqui*2-1),'.4f %' num2str(IdSigEqui*2),'.4f'];
                        fprintf(fileID,Format,mean(Radius_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)*2),std(Radius_E(IdStep,IdDim,IdM,IdTau,IdSigEqui,:)*2));
                    end                    
                end
            end
        end
    end
end




