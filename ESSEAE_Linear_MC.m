%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Experimental Evaluation of ESSEAE with Synthetic Dataset
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
ModelType=0;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 
Rep=10;
%%
initcond=6;             % Initial condition of end-members matrix: 6 (VCA) and 8 (SISAL).
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.1;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=50;
parallel=0;
downsampling=0.0;       % Downsampling in end-members estimation
display_iter=0;            % Display partial performance in BEAE
lm=0.01;
%%
ResultsYh=zeros(length(sSNR),3,Rep);
ResultsAh=zeros(length(sSNR),3,Rep);
ResultsPh=zeros(length(sSNR),3,Rep);
ResultsTh=zeros(length(sSNR),3,Rep);
ResultsPh2=zeros(length(sSNR),3,Rep);
ResultsVh=zeros(length(sSNR),3,Rep);
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
        disp('ESSEAE Analysis');
        Pu=P0(:,rNs);
        
        
        paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,display_iter];
        [P1,A1,S1,Zh1,V1,J1]=HybridEBEAESN(Z,N,paramvec,Pu);
        tesseae=toc;
        ResultsYh(index,1,j)=norm(Zh1-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,1,j)=norm(V1-V0,'fro')/norm(V0,'fro');
        ResultsAh(index,1,j)=errorabundances(A0,A1);
        ResultsPh(index,1,j)=errorendmembers(P0,P1);
        ResultsPh2(index,1,j)=errorSAM(P0,P1);
        ResultsTh(index,1,j)=tesseae;
        
        tic;
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('EBEAE Sup Analysis');
        [P2,A2,S2,Zh2,V2,J2]=EBEAESN(Z,N,paramvec,P0,1);
        tsup=toc;
        ResultsYh(index,2,j)=norm(Zh2-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,2,j)=norm(V2-V0,'fro')/norm(V0,'fro');
        ResultsAh(index,2,j)=errorabundances(A0,A2);
        ResultsPh(index,2,j)=errorendmembers(P0,P2);
        ResultsPh2(index,2,j)=errorSAM(P0,P2);
        ResultsTh(index,2,j)=tsup;
        
        tic;
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('EBEAE Blind Analysis');
        [P3,A3,S3,Zh3,V3,J3]=EBEAESN(Z,N,paramvec);
        tblind=toc;
        ResultsYh(index,3,j)=norm(Zh3-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,3,j)=norm(V3-V0,'fro')/norm(V0,'fro');
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
    AAnd pes num2str(mean(ResultsYh(:,2,:),3),NP) PM num2str(std(ResultsYh(:,2,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsYh(:,1,:),3),NP) PM num2str(std(ResultsYh(:,1,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsYh(:,3,:),3),NP) PM num2str(std(ResultsYh(:,3,:),[],3),NP) pes...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Abundance Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsAh(:,2,:),3),NP) PM num2str(std(ResultsAh(:,2,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsAh(:,1,:),3),NP) PM num2str(std(ResultsAh(:,1,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsAh(:,3,:),3),NP) PM num2str(std(ResultsAh(:,3,:),[],3),NP) pes...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation');
disp([num2str(int8(sSNR')) ...
    AAnd pes num2str(mean(ResultsPh(:,2,:),3),NP) PM num2str(std(ResultsPh(:,2,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsPh(:,1,:),3),NP) PM num2str(std(ResultsPh(:,1,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsPh(:,3,:),3),NP) PM num2str(std(ResultsPh(:,3,:),[],3),NP) pes...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation (SAM)');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsPh2(:,2,:),3),NP) PM num2str(std(ResultsPh2(:,2,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsPh2(:,1,:),3),NP) PM num2str(std(ResultsPh2(:,1,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsPh2(:,3,:),3),NP) PM num2str(std(ResultsPh2(:,3,:),[],3),NP) pes...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Computation Time');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsTh(:,2,:),3),NP) PM num2str(std(ResultsTh(:,2,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsTh(:,1,:),3),NP) PM num2str(std(ResultsTh(:,1,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsTh(:,3,:),3),NP) PM num2str(std(ResultsTh(:,3,:),[],3),NP) pes...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Sparse Noise Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity'))...
    AAnd pes num2str(mean(ResultsVh(:,2,:),3),NP) PM num2str(std(ResultsVh(:,2,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsVh(:,1,:),3),NP) PM num2str(std(ResultsVh(:,1,:),[],3),NP) pes...
    AAnd pes num2str(mean(ResultsVh(:,3,:),3),NP) PM num2str(std(ResultsVh(:,3,:),[],3),NP) pes...
    EEnd]);

%%
disp('%%%%%%%%%%%%%%%%')
disp('ANOVA OUTPUT STIMATION')
[n, algs, reps] = size(ResultsYh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

% Initialize comparison names for only the first algorithm
comparisonNames = {' vs_Sup', '     vs_Blind'};

% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data = squeeze(ResultsYh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data, [], 'off');
    
    % Perform post-hoc comparisons
    results = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova(index, k) = results(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable = array2table(TabAnova, 'VariableNames', comparisonNames);
pValuesTable.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);
% Display the table
disp(pValuesTable);
%%
disp('ANOVA Abundances')
[n, algs, reps] = size(ResultsAh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_A = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_A = squeeze(ResultsAh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_A, [], 'off');
    
    % Perform post-hoc comparisons
    results_A = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_A(index, k) = results_A(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_A = array2table(TabAnova_A, 'VariableNames', comparisonNames);
pValuesTable_A.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_A);
%%
disp('ANOVA EndMembers error')
[n, algs, reps] = size(ResultsPh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Ph = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Ph = squeeze(ResultsPh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Ph, [], 'off');
    
    % Perform post-hoc comparisons
    results_Ph = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Ph(index, k) = results_Ph(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Ph = array2table(TabAnova_Ph, 'VariableNames', comparisonNames);
pValuesTable_Ph.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Ph);
%%
disp('ANOVA EndMembers SAM rror')
[n, algs, reps] = size(ResultsPh2);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Ph2 = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Ph2 = squeeze(ResultsPh2(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Ph2, [], 'off');
    
    % Perform post-hoc comparisons
    results_Ph2 = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Ph2(index, k) = results_Ph2(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Ph2 = array2table(TabAnova_Ph2, 'VariableNames', comparisonNames);
pValuesTable_Ph2.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Ph2);
%%
disp('ANOVA Computational Time')
[n, algs, reps] = size(ResultsTh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Th = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values
% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Th = squeeze(ResultsTh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Th, [], 'off');
    
    % Perform post-hoc comparisons
    results_Th = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Th(index, k) = results_Th(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Th = array2table(TabAnova_Th, 'VariableNames', comparisonNames);
pValuesTable_Th.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Th);
%% 
disp('ANOVA Sparse Error')
[n, algs, reps] = size(ResultsVh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Vh = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Vh = squeeze(ResultsVh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Vh, [], 'off');
    
    % Perform post-hoc comparisons
    results_Vh = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Vh(index, k) = results_Vh(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Vh = array2table(TabAnova_Vh, 'VariableNames', comparisonNames);
pValuesTable_Vh.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Vh);
%%
