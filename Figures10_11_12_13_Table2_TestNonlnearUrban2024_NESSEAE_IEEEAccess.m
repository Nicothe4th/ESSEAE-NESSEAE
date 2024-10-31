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

figure;
Thr=0.6;
labelT={'(a)','(b)','(c)','(d)'};
for i=1:N
    subplot(N,1,i)
    plot(Z(:,A(i,:)>Thr),':'); grid on;
    hold on;
    plot(M(:,i),'r','LineWidth',4); axis tight;
    ylabel('Normalized intensity')
    title([labelT(i), cood(i)]);
end
xlabel('Spectral channel');


figure;
for i=1:N
    subplot(1,N,i);
    imagesc(reshape(A(i,:),nRow,nCol),[0,1]);
    title([labelT(i), cood(i)])
end
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

disp('%%%%%%%%%%%%%%%%%')
disp('NESSEAE');
index=2;
IndexUnknown=zeros(2,2+N*(N-1)/2);

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
    
        if (ii==1 && jj==2) || (ii==2 && jj==3)

            figure;
            subplot(2,1,1)
            plot(M(:,UnknownEM),'linewidth',2); grid on;
            title(['(a) Ground-truth' cood(UnknownEM(1)) cood(UnknownEM(2))])
            axis([1 size(M,1) min([min(M(:)), min(P1(:))]) max([max(M(:)), max(P1(:))])])
            subplot(2,1,2);  
            plot(P1(:,F+1:N),'linewidth',2); grid on;
            title('(b) Estimated End-members')
            axis([1 size(M,1) min([min(M(:)) min(P1(:))]) max([max(M(:)) max(P1(:))])])
    
            figure;
            Au=A(UnknownEM,:);
            A1u=A1(F+1:N,:);
            for u=1:U
                subplot(2,U,u);
                imagesc(reshape(Au(u,:),nRow,nCol),[0,1]);
                title(['(a) Ground-truth' cood(UnknownEM(u))])
                subplot(2,U,u+U);
                imagesc(reshape(A1u(u,:),nRow,nCol),[0,1]);
                title(['(b) Estimation' num2str(u) '-th'])
            end
        end
    end
end
%%
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
