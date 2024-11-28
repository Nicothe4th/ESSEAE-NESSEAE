%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Analysis of Urban Dataset with 4 End-members with NESSEAE
%
% Daniel Ulises Campos-Delgado & Nicolas Mendoza Chavarria
% October/2024
% UASLP-ULPGC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('EBEAE');

load('Urban/Urban_F210.mat');
N=4;
load(['Urban/end' num2str(N) '_groundTruth.mat']);

M=M./sum(M,1);
A=A./sum(A,1);
Yo=Y(SlectBands,:);
Z=Yo./sum(Yo,1);
[L,K]=size(Z);
F=2;
U=2;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Gorund-Truths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure;
Thr=0.6;
labelT={'(a)','(b)','(c)','(d)'};
channels=1:162;
for i=1:N
    subplot(N,1,i)
    data=Z(:,A(i,:)>Thr);
    dataMax=max(data,[],2);
    dataMin=min(data,[],2);
    x2=[channels, fliplr(channels)];
    inBetween = [dataMax' fliplr(dataMin')];
    patch(x2, inBetween, 'r','FaceAlpha',0.2, 'EdgeColor','none');
    %plot(Z(:,A(i,:)>Thr),':'); grid on;
    hold on;
    plot(M(:,i),'r','LineWidth',4); axis tight; grid on;
    ylabel('Normalized intensity')
    title([labelT(i), cood(i)], 'FontSize',14,'FontWeight','normal');
end
xlabel('Spectral channel','FontSize',14);

%%
figure;
subplot(4,2,1);
plot(M,'LineWidth',2);axis tight; grid on;
xlabel('spectral channel','FontSize',12);
ylabel('Normalized intensity','FontSize',12); 
title('(a)', 'FontSize',14,'FontWeight','normal');
legend(cood)
for i=1:N
    subplot(4,8,4+i);
    imagesc(reshape(A(i,:),nRow,nCol),[0,1]); axis off;
    title(cood(i), 'FontSize',12,'FontWeight','normal')
end
annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'string', '(b)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of NEBEAE and NESSEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=5;             % Initial condition of end-members matrix
rho=1;               % Similarity weight in end-members estimation
lambda=0.1;             % Entropy weight for abundance estimation
lm=0.1;
epsilon=1e-3;
maxiter=50;
parallel=1;
normalization=1;
downsampling=0.0;      
display_iter=0;            

paramvecSN=[initcond,rho,lambda,lm, epsilon,maxiter,downsampling,parallel,display_iter];

rho2=1;
lambda2=0.25;
paramvec=[initcond,rho2,lambda2,epsilon,maxiter,downsampling,parallel,display_iter];
paramvecSN2=[initcond,rho2,lambda2,lm,epsilon,maxiter,downsampling,parallel,display_iter];


ErrorP=zeros(2+N*(N-1)/2,1);
ErrorSAM=zeros(2+N*(N-1)/2,1);
ErrorA=zeros(2+N*(N-1)/2,1);
ErrorZ=zeros(2+N*(N-1)/2,1);
CompTime=zeros(2+N*(N-1)/2,1);
EnergyV=zeros(2+N*(N-1)/2,1);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Urban Dataset');
disp(['Number of fixed end-members=' num2str(F)]);
disp(['Number of unknown end-members=' num2str(U)]);

tic;
[P2,A2,D2,S2,Zh2,J2]=NEBEAE3(Z,N,paramvec,M,1);
Tnebeae=toc; 

disp('%%%%%%%%%%%%%%%%%')
disp('NEBEAE-Supervised');
ErrorZ(1)=norm(Zh2-Z,'fro')/norm(Z,'fro');
ErrorP(1)=0;
ErrorSAM(1)=0;
ErrorA(1)=errorabundances(A,A2);
CompTime(1)=Tnebeae;
EnergyV(1)=0;
%%

subplot(4,2,3);
plot(P2,'LineWidth',2); axis tight; grid on;
xlabel('spectral channel','FontSize',12);
ylabel('Normalized intensity','FontSize',12); 
title('(c)', 'FontSize',14,'FontWeight','normal');
legend('Fixed-1','Fixed-2','Fixed-3','Fixed-4')
for i=1:N
    subplot(4,8,12+i);
    imagesc(reshape(A2(i,:),nRow,nCol),[0,1]); axis off;
    title(['Fixed-' num2str(i)], 'FontSize',12,'FontWeight','normal')
end
annotation('textbox', [0.7, 0.65, 0.1, 0.1], 'string', '(d)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');
%%

tic;
[P3,A3,D3,S3,Zh3,V3,J3]=NEBEAESN(Z,N,paramvecSN2);
Tnebeaesn=toc; 

disp('%%%%%%%%%%%%%%%%%')
disp('NEBEAE');
ErrorZ(2+N*(N-1)/2)=norm(Zh3-Z,'fro')/norm(Z,'fro');
ErrorP(2+N*(N-1)/2)=errorendmembers(M,P3);
ErrorSAM(2+N*(N-1)/2)=errorSAM(M,P3);
ErrorA(2+N*(N-1)/2)=errorabundances(A,A3);
CompTime(2+N*(N-1)/2)=Tnebeaesn;
EnergyV(2+N*(N-1)/2)=norm(V3,'fro')/norm(Z,'fro'); 
%%

subplot(4,2,7);
plot(P3,'LineWidth',2); axis tight; grid on;
xlabel('spectral channel','FontSize',12);
ylabel('Normalized intensity','FontSize',12); 
title('(g)', 'FontSize',14,'FontWeight','normal');
legend('Unknown-1','Unknown-2','Unknown-3','Unknown-4')
for i=1:N
    subplot(4,8,28+i);
    imagesc(reshape(A3(i,:),nRow,nCol),[0,1]); axis off;
    title(['Unknown-' num2str(i)], 'FontSize',12,'FontWeight','normal')
end
annotation('textbox', [0.7, 0.25, 0.1, 0.1], 'string', '(h)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');
%%

disp('%%%%%%%%%%%%%%%%%')
disp('NESSEAE');
index=2;
IndexUnknown=zeros(2,2+N*(N-1)/2);
%%
for ii=1:(N-1)
    for jj=ii+1:N
        disp(['Fixed end-members :' num2str(ii) ' and ' num2str(jj)]);
        FixedEM=[ii,jj];
        IndexUnknown(:,index)=FixedEM';
        UnknownEM=setdiff(1:N,FixedEM);

        AoF=A(FixedEM,:);
        AoU=A(UnknownEM,:);
        PF=M(:,FixedEM);
        PU=M(:,UnknownEM);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Execute NESSEAE Methodology
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tic;
        [P1,A1,D1,S1,Zh1,V1,J1]=NESSEAE(Z,N,paramvecSN,PF);
        Tnesseae=toc;  

        ErrorZ(index)=norm(Zh1-Z,'fro')/norm(Z,'fro');
        ErrorP(index)=errorendmembers(M,P1);
        ErrorSAM(index)=errorSAM(M,P1);
        ErrorA(index)=errorabundances(A,A1);
        CompTime(index)=Tnesseae;
        EnergyV(index)=norm(V1,'fro')/norm(Z,'fro'); 
        index=index+1;
    
        if (ii==2 && jj==3)

            subplot(4,2,5);
            plot(P1,'LineWidth',2); axis tight; grid on;
            xlabel('spectral channel','FontSize',12);
            ylabel('Normalized intensity','FontSize',12); 
            title('(e)', 'FontSize',14,'FontWeight','normal');
            legend('Fixed-1','Fixed-2','Unknown-1','Unknown-2')
            for i=1:N
                subplot(4,8,20+i);
                imagesc(reshape(A1(i,:),nRow,nCol),[0,1]); axis off;
                if i<3
                    title(['Fixed-' num2str(i)], 'FontSize',12,'FontWeight','normal');
                else
                    title(['Unknown-' num2str(i-2)], 'FontSize',12,'FontWeight','normal')
                end
            end
            annotation('textbox', [0.7, 0.45, 0.1, 0.1], 'string', '(f)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');

        end
    end
end

disp('%%%%%%%%%%%%%%%%%%%')
disp('Peformance Metrics')
disp('NEBEAE-Sup            NESSEAE               NEBEAE')
disp(IndexUnknown)
disp('Measurements Errors E_z')
disp(ErrorZ')
disp('Abundance Errors E_a')
disp(ErrorA')
disp('End-members Errors E_p')
disp(ErrorP')
disp('End-members SAM E_SAM')
disp(ErrorSAM')
disp('Energy Spectral Noise')
disp(EnergyV')
disp('Computational Time')
disp(CompTime')
