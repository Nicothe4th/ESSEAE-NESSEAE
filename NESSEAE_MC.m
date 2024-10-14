%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Experimental Evaluation of NESSEAE with Synthetic Dataset
% 
%
% SEP/2024
% JNMC-DUCD-UASLP
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('EBEAE');
sSNR=[40 35 30];%[30 35 40];    
pDensity=[0.005 0.0075 0.01];%[0.01 0.0075 0.005];

N=4;                % Number of End-members
Nsamples=64;
nCol=Nsamples;
nRow=Nsamples;
ModelType=5;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 
Rep=10;
%%
rho=0.1;                   % Weight on Regularization Term  >= 0
lambda=0.1;                                                                                                                                                            % Weight on Entropy of Abundances \in [0,1)
epsilon=1e-3;            % Threshold for convergence
maxiter=50;              % Maximum number of iteration in alternated least squares approach
downsampling=0.0;       % Downsampling factor to estimate end-members \in [0,1)
parallel=0;              % Parallelization in abundance estimation process
disp_iter=0;          % Display results of iterative optimization
initcond=5;
lm=0.01;
%%
ResultsYh=zeros(length(sSNR),3,Rep);
ResultsAh=zeros(length(sSNR),3,Rep);
ResultsPh=zeros(length(sSNR),3,Rep);
ResultsTh=zeros(length(sSNR),3,Rep);
ResultsPh2=zeros(length(sSNR),3,Rep);
ResultsVh=zeros(length(sSNR),3,Rep);
ResultsDh=zeros(length(sSNR),3,Rep);


for index=1:length(sSNR)

    SNR=sSNR(index);
    density=pDensity(index);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('Synthetic Datasets');
    disp('Absorbance HSI');
    disp('Semi-Supervised Unmixing Estimation');
    disp(['SNR =' num2str(SNR) ' dB']);
    disp(['density =' num2str(density) ]);
    disp(['Number of end-members=' num2str(N)]);
    for j=1:Rep 
        [Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
        titleL='Synthetic Dataset with Absorbance Spectral Response of Hb, HbO2, Fat and Water';
        rNs =  sort(randperm(4, 2));
        P0=P0./sum(P0);
        [L,K]=size(Z);
        disp(['Iteration=' num2str(j)])   

        tic;
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('NESSEAE Analysis');
        Pu=P0(:,rNs);
        
        
        paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,disp_iter];
        [P1,A1,D1,S1,Zh1,V1,J1]=HybridNEBEAESN(Z,N,paramvec,Pu);
        tesseae=toc;
        ResultsYh(index,1,j)=norm(Zh1-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,1,j)=norm(V1-V0,'fro')/norm(V0,'fro');
        ResultsDh(index,1,j)=norm(D1-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,1,j)=errorabundances(A0,A1);
        ResultsPh(index,1,j)=errorendmembers(P0,P1);
        ResultsPh2(index,1,j)=errorSAM(P0,P1);
        ResultsTh(index,1,j)=tesseae;
        
        tic;
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('NEBEAE Sup Analysis');
        [P2,A2,D2,S2,Zh2,V2,J2]=NEBEAESN(Z,N,paramvec,P0,1);
        tsup=toc;
        ResultsYh(index,2,j)=norm(Zh2-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,2,j)=norm(V2-V0,'fro')/norm(V0,'fro');
        ResultsDh(index,2,j)=norm(D2-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,2,j)=errorabundances(A0,A2);
        ResultsPh(index,2,j)=errorendmembers(P0,P2);
        ResultsPh2(index,2,j)=errorSAM(P0,P2);
        ResultsTh(index,2,j)=tsup;
        
        tic;
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('NEBEAE Blind Analysis');
        [P3,A3,D3,S3,Zh3,V3,J3]=NEBEAESN(Z,N,paramvec);
        tblind=toc;
        ResultsYh(index,3,j)=norm(Zh3-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,3,j)=norm(V3-V0,'fro')/norm(V0,'fro');
        ResultsDh(index,3,j)=norm(D3-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,3,j)=errorabundances(A0,A3);
        ResultsPh(index,3,j)=errorendmembers(P0,P3);
        ResultsPh2(index,3,j)=errorSAM(P0,P3);
        ResultsTh(index,3,j)=tblind;
        
        
        
    end
end
%%
NP=4;
AAnd=[]; PM=[]; EEnd=[]; slash=[]; pes=[];
for jj=1:length(sSNR)
        AAnd=[AAnd; ' & '];
        PM=[PM; ' \pm '];
        EEnd=[EEnd; ' \\'];
        slash=[slash; '/'];
        pes=[pes; '$'];
end
%%
clc
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Mean Responses in Performance Metrics')
disp('SNR/Density      Sup       Semi-supervised         Blind     ');
disp('%%%%%%%%%%%%%%%');
disp('Error in Output Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) ...
    AAnd pes num2str(mean(ResultsYh(:,2,:),3),NP) PM num2str(std(ResultsYh(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsYh(:,1,:),3),NP) PM num2str(std(ResultsYh(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsYh(:,3,:),3),NP) PM num2str(std(ResultsYh(:,3,:),[],3),NP) pes...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Abundance Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsAh(:,2,:),3),NP) PM num2str(std(ResultsAh(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsAh(:,1,:),3),NP) PM num2str(std(ResultsAh(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsAh(:,3,:),3),NP) PM num2str(std(ResultsAh(:,3,:),[],3),NP) pes  ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Non-linear interaction levels Estimation');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsDh(:,2,:),3),NP) PM num2str(std(ResultsDh(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsDh(:,1,:),3),NP) PM num2str(std(ResultsDh(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsDh(:,3,:),3),NP) PM num2str(std(ResultsDh(:,3,:),[],3),NP) pes ...
    EEnd]);

disp('Error in End-member Estimation');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsPh(:,2,:),3),NP) PM num2str(std(ResultsPh(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsPh(:,1,:),3),NP) PM num2str(std(ResultsPh(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsPh(:,3,:),3),NP) PM num2str(std(ResultsPh(:,3,:),[],3),NP) pes ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation (SAM)');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsPh2(:,2,:),3),NP) PM num2str(std(ResultsPh2(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsPh2(:,1,:),3),NP) PM num2str(std(ResultsPh2(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsPh2(:,3,:),3),NP) PM num2str(std(ResultsPh2(:,3,:),[],3),NP) pes ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Computation Time');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsTh(:,2,:),3),NP) PM num2str(std(ResultsTh(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsTh(:,1,:),3),NP) PM num2str(std(ResultsTh(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsTh(:,3,:),3),NP) PM num2str(std(ResultsTh(:,3,:),[],3),NP) pes ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Sparse Noise Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsVh(:,2,:),3),NP) PM num2str(std(ResultsVh(:,2,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsVh(:,1,:),3),NP) PM num2str(std(ResultsVh(:,1,:),[],3),NP) pes ...
    AAnd pes num2str(mean(ResultsVh(:,3,:),3),NP) PM num2str(std(ResultsVh(:,3,:),[],3),NP) pes ...
    EEnd]);

%%


